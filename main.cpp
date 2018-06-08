#include <iostream>
#include <unordered_map>
#include <map>
#include "minimizer.h"
#include "filter.h"
#include <algorithm>
#include <vector>
#include <cassert>
#include "opal.h"
#include "ScoreMatrix.hpp"
#include <ctime>
#include <unordered_set>
#include <cmath>
#include <cstring>

#define W 13
#define K 4
#define PERC 1

int OUTPUT_RESULT = 10;


using namespace std;

extern "C" {
    #include "fasta.h"
}

unordered_set<int> freqHashes;

void reduceIndexTable(minimizer::IndexTable& indexTable, double p) {
    assert(p <= 100 && p >= 0);
    int indexCnt = 0;
    for (auto it = indexTable.begin(); it != indexTable.end(); it++) {
        indexCnt += it->size();
    }
    int limit = (p/100)*indexCnt;

    int hash = 0;
    for (auto it = indexTable.begin(); it != indexTable.end(); it++) {
        if (it->size() > limit) {
            it->clear();
            freqHashes.insert(hash);
        }
        hash++;
    }

}

double toSeconds(clock_t t) {
    return ((double)t) / CLOCKS_PER_SEC;
}

int main(int argc, char **argv) {
        clock_t start = clock();
        FASTAFILE *ffp;
        ffp = OpenFASTA(argv[1]);
        FASTAFILE *ffq;
        ffq = OpenFASTA(argv[2]);
        freqHashes.clear();

        if (argc > 3) {
            OUTPUT_RESULT = atoi(argv[3]);
        }
        long long finalResult = 0;


        char *seq, *name;
        int len;
        int dbSeqCnt = 0;
        vector<char *> db;
        vector<string> dbNames;
        vector<int> dbLengths;

        int queryCnt = 0;
        vector<char *> queries;
        vector<string> queryName;
        vector<int> queryLength;

        minimizer::IndexTable indexTable;

        clock_t t = clock();
        long long size = 0;
        while (ReadFASTA(ffp, &seq, &name, &len)) {
            db.push_back(seq);
            string n = "";
            for (char *c = name;*c && *c != ' ' && *c != '\n' && *c != '\r'; c++) {
                if (*c == ' ') {
                    *c = 0;
                    break;
                }
                n += *c;

            }
            size += len;
            dbNames.push_back(n);
            dbLengths.push_back(len);
            dbSeqCnt++;
        }
        cout << dbNames.size() <<  "proteina u bazi " << size * 1. / dbNames.size() <<  endl;

        CloseFASTA(ffp);

        fprintf(stderr, "citanje baze od %d proteina iz faste: %lf s\n", dbSeqCnt, toSeconds(clock() - t));
        t = clock();
        for (int i = 0; i < dbSeqCnt; i++) {
            minimizer::addMinimizers(db[i], dbLengths[i], i, W, K, indexTable);
        }
        fprintf(stderr, "indeksacija baze: %lf s\n", toSeconds(clock() - t));
        t = clock();

        int cnt = 0;
        while (ReadFASTA(ffq, &seq, &name, &len)) {
            cnt++;
            queries.push_back(seq);

            string n = "";
            for (char *c = name;*c && *c != ' ' && *c != '\n' && *c != '\r'; c++) {
                if (*c == ' ') {
                    *c = 0;
                    break;
                }
                n += *c;

            }
            queryName.push_back(n);
            queryLength.push_back(len);
        }
        CloseFASTA(ffq);

        fprintf(stderr, "citanje %d queryija iz faste: %lf s\n", queries.size(), toSeconds(clock() - t));
        t = clock();

        reduceIndexTable(indexTable, (PERC));
        for (int i = 0; i < indexTable.size(); i++) {
            sort(indexTable[i].begin(), indexTable[i].end());
        }
        fprintf(stderr, "reduciranje cestih i sorrtiranje minimizera: %lf s\n", toSeconds(clock() - t));
        t = clock();


        //     OPAL CONFIG
        int alphabetLength = 23;
        int gapOpen = 11;
        int gapExt = 1;


        unsigned char blosumAlphabet[23] = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
                                            'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X'};
        unsigned char toBlosumIndex[26];
        for (int i = 0; i < 26; i++)
            toBlosumIndex[i] = 22;
        for (int i = 0; i < 23; i++) {
            toBlosumIndex[blosumAlphabet[i] - 'A'] = i;

        }

        int *scoreMatrix = ScoreMatrix::getBlosum62().getMatrix();




        // za svaki query napravi redukciju baze i opal nad kandidatima

        clock_t qs = clock();
        for (int q = 0; q < queries.size(); q++) {
            int reduceTo = std::min((int)dbNames.size() / 100, 1500000 / queryLength[q]);

            t = clock();
            vector<minimizer::Minimizer> mins = minimizer::computeForSequence(queries[q], queryLength[q], W, K);
            vector<minimizer::Minimizer> seqMin;
            for (minimizer::Minimizer mini: mins) {
                if (freqHashes.find(mini.h) == freqHashes.end())
                    seqMin.push_back(mini);
            }
            vector<int> similar = filter::getSimilar(seqMin, indexTable, reduceTo);
            fprintf(stderr, "prostor reduciran na %d proteina u: %lf s\n", similar.size(), toSeconds(clock() - t));
            t = clock();
            int opStart = clock();
            fprintf(stderr, "minimizera ima: %d\n", seqMin.size());


            //sortiranje similar sekvence tako da prvo dolaze manji proteini
            std::sort(similar.begin(), similar.end(), [&](int a, int b) {return dbLengths[a] < dbLengths[b];});


// Query

            unsigned char *opalQuery = new unsigned char[sizeof(unsigned char) * queryLength[q]];
            for (int i = 0; i < queryLength[q]; i++) {
                opalQuery[i] = toBlosumIndex[queries[q][i] - 'A'];
            }

// Database
            int dbOpalLength = similar.size();

            unsigned char **opalDb = new unsigned char *[sizeof(unsigned char *) * dbOpalLength];
            int *opalDbSeqsLen = new int[sizeof(int) * dbOpalLength];
            int opdb = 0;
            for (int candidate: similar) {
                unsigned char *opalCandidate = new unsigned char[sizeof(unsigned char) * dbLengths[candidate]];
                for (int i = 0; i < dbLengths[candidate]; i++) {
                    opalCandidate[i] = toBlosumIndex[db[candidate][i] - 'A'];
                }
                opalDb[opdb] = opalCandidate;
                opalDbSeqsLen[opdb] = dbLengths[candidate];
                opdb++;
            }

// Results for each sequence in database
            OpalSearchResult **results = new OpalSearchResult *[sizeof(OpalSearchResult *) * dbOpalLength];
            for (int i = 0; i < dbOpalLength; i++) {
                results[i] = new OpalSearchResult();
                opalInitSearchResult(results[i]);
            }


// Do calculation!
            t = clock();
            int resultCode = opalSearchDatabase(opalQuery, queryLength[q], opalDb, dbOpalLength, opalDbSeqsLen,
                                                gapOpen, gapExt, ScoreMatrix::getBlosum62().getMatrix(),
                                                alphabetLength,
                                                results,
                                                OPAL_SEARCH_ALIGNMENT, OPAL_MODE_SW, OPAL_OVERFLOW_SIMPLE);

            fprintf(stderr, "kraj opala u: %lf s\n", toSeconds(clock() - t));

            t = clock();
            delete opalQuery;
            for (int i = 0; i < dbOpalLength; i++)
                delete opalDb[i];
            delete opalDbSeqsLen;
            delete opalDb;

// Print scores
            vector<pair<int, int>> topScoring;
            for (int i = 0; i < similar.size(); i++) {
                topScoring.push_back({results[i]->score, i});
            }

            sort(topScoring.begin(), topScoring.end());
            fprintf(stderr, "\n\n");
            int cnt = 0;
            for (int i = (int) topScoring.size() - 1; i >= 0 && cnt < OUTPUT_RESULT; i--) {
                cnt++;
                int matchCnt = 0;
                int mismatchCnt = 0;
                int gapCnt = 0;
                int id = topScoring[i].second;
                for (int j = 0; j < results[id]->alignmentLength; j++) {
                    switch (results[id]->alignment[j]) {
                        case OPAL_ALIGN_MATCH:
                            matchCnt++;
                            break;
                        case OPAL_ALIGN_DEL:
                            if (j == 0 || results[id]->alignment[j-1] != OPAL_ALIGN_DEL)
                                gapCnt++;
                            break;
                        case OPAL_ALIGN_INS:
                            if (j == 0 || results[id]->alignment[j-1] != OPAL_ALIGN_INS)
                                gapCnt++;
                            break;
                        case OPAL_ALIGN_MISMATCH:
                            mismatchCnt++;
                            break;
                    }
                }

                printf("%s\t%s\t%.2lf\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                       queryName[q].c_str(),
                       dbNames[similar[topScoring[i].second]].c_str(),
                       matchCnt * 100. / results[id]->alignmentLength,
                       results[id]->alignmentLength,
                       mismatchCnt,
                       gapCnt,
                       results[id]->startLocationQuery,
                       results[id]->endLocationQuery,
                       results[id]->startLocationTarget,
                       results[id]->endLocationTarget,
                       topScoring[i].first);
                finalResult += topScoring[i].first;
            }
            fprintf(stderr, "\n\n");

            for (int i = 0; i < dbOpalLength; i++) {
                free(results[i]->alignment);
                delete results[i];
            }
            delete results;
        }

        fprintf(stderr, "queryije odgovorio u: %lf s \n", toSeconds(clock() - qs));
        fprintf(stderr, "total kraj: %lf s \n", toSeconds(clock() - start));
        fprintf(stderr, "FINAL: %lld\n", finalResult);
    return 0;
}
