#include <iostream>
#include <unordered_map>
#include <map>
#include "minimizer.h"
#include "lis.h"
#include <algorithm>
#include <vector>
#include <cassert>
#include "opal.h"
#include "ScoreMatrix.hpp"
#include <ctime>
#include <unordered_set>

#define W 10
#define K 3
#define OUTPUT_RESULT 5

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

int main() {
    clock_t start = clock();
    FASTAFILE *ffp;
    ffp = OpenFASTA("../uniprot_sprot3MB.fasta");
    FASTAFILE *ffq;
    ffq = OpenFASTA("../uniprot_sprotS.fasta");


    char *seq, *name;
    int len;
    int dbSeqCnt = 0;
    vector<char*> db;
    vector<char*> dbNames;
    vector<int> dbLengths;

    int queryCnt = 0;
    vector<char*> queries;
    vector<char*> queryName;
    vector<int> queryLength;

    minimizer::IndexTable indexTable;

    clock_t t = clock();
    while (ReadFASTA(ffp, &seq, &name, &len)) {
        db.push_back(seq);
        for (char* c = name; ; c++) {
            if (*c == ' ') {
                *c = 0;
                break;
            }
        }
        dbNames.push_back(name);
        dbLengths.push_back(len);
        dbSeqCnt++;
    }

    CloseFASTA(ffp);

    printf("citanje baze od %d proteina iz faste: %lf s\n", dbSeqCnt, toSeconds(clock() - t));
    t = clock();
    for (int i = 0; i < dbSeqCnt; i++) {
        minimizer::addMinimizers(db[i], dbLengths[i], i, W, K, indexTable);
    }
    printf("indeksacija baze: %lf s\n", toSeconds(clock() - t));
    t = clock();

    while (ReadFASTA(ffq, &seq, &name, &len)) {
        queries.push_back(seq);

        for (char* c = name; ; c++) {
            if (*c == ' ') {
                *c = 0;
                break;
            }
        }
        queryName.push_back(name);
        queryLength.push_back(len);
    }
    CloseFASTA(ffq);

    printf("citanje %d queryija iz faste: %lf s\n", queries.size(), toSeconds(clock() - t));
    t = clock();

    reduceIndexTable(indexTable, (1));
    printf("reduciranje cestih minimizera: %lf s\n", toSeconds(clock() - t));
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

    for (int q = 0; q < queries.size(); q++) {
        t = clock();
        vector<minimizer::Minimizer> mins = minimizer::computeForSequence(queries[q], queryLength[q], W, K);
        vector<minimizer::Minimizer> seqMin;
        for (minimizer::Minimizer mini: mins) {
            if (freqHashes.find(mini.h) == freqHashes.end())
                seqMin.push_back(mini);
        }
        vector<int> similar = lis::getSimilar(seqMin, indexTable);
        printf("prostor reduciran na %d proteina u: %lf s\n", similar.size(), toSeconds(clock() - t));
        t = clock();
        int opStart = clock();
        printf("minimizera ima: %d\n", seqMin.size());





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
                                            gapOpen, gapExt, ScoreMatrix::getBlosum62().getMatrix(), alphabetLength,
                                            results,
                                            OPAL_SEARCH_SCORE, OPAL_MODE_SW, OPAL_OVERFLOW_BUCKETS);

        printf("kraj opala u: %lf s\n", toSeconds(clock() - t));

        t = clock();
        delete opalQuery;
        for (int i = 0; i < dbOpalLength; i++)
            delete opalDb[i];
        delete opalDbSeqsLen;
        delete opalDb;

// Print scores
        int i = 0;
        vector<pair<int, int>> result;
        for (int cand: similar) {
            result.push_back({results[i]->score, cand});
            i++;
        }

        sort(result.begin(), result.end());
        int cnt = 0;
        printf("\n\n");
        for (int i = (int)result.size() - 1; i >= 0 && cnt < OUTPUT_RESULT; i--) {
            cnt++;
            printf("%s\t%s\t%d\n", queryName[q], dbNames[result[i].second], result[i].first);
        }
        printf("\n\n");

        for (int i = 0; i < dbOpalLength; i++)
            delete results[i];
        delete results;
    }

    printf("total kraj: %lf s \n", toSeconds(clock() - start));
    return 0;
}