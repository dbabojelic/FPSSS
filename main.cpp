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

#define PERC 1

int W, K;
int OUTPUT_RESULT = 10;
int BANDS;
int REDUCE_TO_PERC, REDUCE_TO_PER_QUERY_LEN;


using namespace std;

extern "C" {
#include "fasta.h"
}

vector<unordered_set<int>> freqHashes;

void reduceIndexTable(minimizer::IndexTable* indexTable, double p) {
    assert(p <= 100 && p >= 0);
    for (int band = 0; band < BANDS; band++) {
        int indexCnt = 0;
        for (auto it = indexTable[band].begin(); it != indexTable[band].end(); it++) {
            indexCnt += it->size();
        }
        int limit = (p / 100) * indexCnt;

        int hash = 0;
        for (auto it = indexTable[band].begin(); it != indexTable[band].end(); it++) {
            if (it->size() > limit) {
                it->clear();
                freqHashes[band].insert(hash);
            }
            hash++;
        }
    }

}

double toSeconds(clock_t t) {
    return ((double) t) / CLOCKS_PER_SEC;
}

int main(int argc, char **argv) {
    if (argc < 9) {
        printf("Trebate unijeti 9 argumenata: db_file, query_file, best_cnt, window_size, k-mer_size, bands, reduce_to_perc, reduce_to_per_query_len\n");
        return 0;
    }

    clock_t start = clock();
    FASTAFILE *ffp;
    ffp = OpenFASTA(argv[1]);
    FASTAFILE *ffq;
    ffq = OpenFASTA(argv[2]);
    freqHashes.clear();
    OUTPUT_RESULT = atoi(argv[3]);
    W = atoi(argv[4]);
    K = atoi(argv[5]);
    BANDS = atoi(argv[6]);
    REDUCE_TO_PERC = atoi(argv[7]);
    REDUCE_TO_PER_QUERY_LEN = atoi(argv[8]);

    fprintf(stderr, "Ucitao sve cli argumente!\n");

    freqHashes.resize(BANDS);

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

    minimizer::IndexTable* indexTable = new minimizer::IndexTable[BANDS];

    clock_t t = clock();
    long long size = 0;
    while (ReadFASTA(ffp, &seq, &name, &len)) {
        db.push_back(seq);
        string n = "";
        for (char *c = name; *c && *c != ' ' && *c != '\n' && *c != '\r'; c++) {
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
    CloseFASTA(ffp);

    fprintf(stderr, "citanje baze od %d proteina iz faste: %lf s\n", dbSeqCnt, toSeconds(clock() - t));
    t = clock();
    for (int i = 0; i < dbSeqCnt; i++) {
        minimizer::addMinimizers(db[i], dbLengths[i], i, W, K, indexTable, BANDS);
    }
    fprintf(stderr, "indeksacija baze: %lf s\n", toSeconds(clock() - t));
    t = clock();

    int cnt = 0;
    while (ReadFASTA(ffq, &seq, &name, &len)) {
        cnt++;
        queries.push_back(seq);

        string n = "";
        for (char *c = name; *c && *c != ' ' && *c != '\n' && *c != '\r'; c++) {
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
    for (int band = 0; band < BANDS; band++) {
        for (int i = 0; i < indexTable[band].size(); i++) {
            sort(indexTable[band][i].begin(), indexTable[band][i].end());
        }
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
    double filtering = 0;
    double opalfunc = 0;
    fprintf(stderr, "cijela indeksacija i sav posao napravljen u: %lf s\n", toSeconds(qs - start));
    int rt = (int)dbNames.size() * REDUCE_TO_PERC / 100;
    for (int q = 0; q < queries.size(); q++) {
        int reduceTo = std::min(rt, REDUCE_TO_PER_QUERY_LEN / queryLength[q]);

        t = clock();
        vector<vector<minimizer::Minimizer>> mins = minimizer::computeForSequence(queries[q], queryLength[q], W, K, BANDS);
        vector<vector<minimizer::Minimizer>> seqMin(BANDS, vector<minimizer::Minimizer>());
        for (int band = 0; band < BANDS; band++) {
            for (minimizer::Minimizer mini: mins[band]) {
                if (freqHashes[band].find(mini.h) == freqHashes[band].end())
                    seqMin[band].push_back(mini);
            }
        }
        vector<int> similar = filter::getSimilar(seqMin, indexTable, BANDS, reduceTo);
        //sortiranje similar sekvence tako da prvo dolaze manji proteini
        std::sort(similar.begin(), similar.end(), [&](int a, int b) { return dbLengths[a] < dbLengths[b]; });

        double filteringDone = toSeconds(clock() - t);
        filtering += filteringDone;
        fprintf(stderr, "prostor reduciran na %d proteina u: %lf s\n", similar.size(), filteringDone);


        t = clock();
        int opStart = clock();
        fprintf(stderr, "minimizera ima: %d\n", seqMin.size());




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

        double opalfuncDone = toSeconds(clock() - t);
        opalfunc += opalfuncDone;
        fprintf(stderr, "kraj opala u: %lf s\n", opalfuncDone);

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
                        if (j == 0 || results[id]->alignment[j - 1] != OPAL_ALIGN_DEL)
                            gapCnt++;
                        break;
                    case OPAL_ALIGN_INS:
                        if (j == 0 || results[id]->alignment[j - 1] != OPAL_ALIGN_INS)
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

    double cijeloOdgovaranje = toSeconds(clock() - qs);
    fprintf(stderr, "cijelo odgovaranje: %lf , filtriranje: %lf, opalfunc %lf\n", cijeloOdgovaranje, filtering,
            opalfunc);

    fprintf(stderr, "queryije odgovorio u: %lf s \n", toSeconds(clock() - qs));
    fprintf(stderr, "total kraj: %lf s \n", toSeconds(clock() - start));
    fprintf(stderr, "FINAL: %lld\n", finalResult);
    return 0;
}
