#include <iostream>
#include <chrono>
#include <unordered_map>
#include <map>
#include <map>
#include "minimizer.h"
#include "lis.h"
#include <algorithm>
#include <vector>
#include <cassert>
#include "opal.h"
#include "ScoreMatrix.hpp"

using namespace std;

extern "C" {
    #include "fasta.h"
}

void reduceIndexTable(minimizer::IndexTable& indexTable, double p) {
    assert(p <= 100 && p >= 0);
    int indexCnt = 0;
    for (auto it = indexTable.begin(); it != indexTable.end(); it++) {
        indexCnt += it->size();
    }
    int limit = (p/100)*indexCnt;

    for (auto it = indexTable.begin(); it != indexTable.end(); it++) {
        if (it->size() > limit) {
            it->clear();
        }
    }

}

int main() {
    FASTAFILE *ffp;
    ffp = OpenFASTA("../uniprot_sprot.fasta");

    chrono::steady_clock::time_point start = chrono::steady_clock::now();

    char *seq, *name;
    int len;
    int cnt = 0;
    vector<char*> db;
    vector<char*> dbNames;
    vector<int> dbLengths;
    char *query;
    char *queryName;
    int queryLength;

    minimizer::IndexTable indexTable;

    vector<minimizer::Minimizer> seqMin;
    while (ReadFASTA(ffp, &seq, &name, &len)) {
        db.push_back(seq);
        dbNames.push_back(name);
        dbLengths.push_back(len);
        minimizer::addMinimizers(seq, len, cnt, 10, 3, indexTable);
        if (cnt == 7) {
            seqMin = minimizer::computeForSequence(seq, len, 10, 3);
            query = seq;
            queryName = name;
            queryLength = len;
        }
        cnt++;
    }

    CloseFASTA(ffp);

    chrono::steady_clock::time_point indexingEnd = chrono::steady_clock::now();
    cout << "indexing finished in: " << chrono::duration_cast<chrono::seconds>(indexingEnd - start).count() << endl;

    reduceIndexTable(indexTable, (1));
    chrono::steady_clock::time_point reduceingEnd = chrono::steady_clock::now();
    cout << "reduceing end in: " << chrono::duration_cast<chrono::seconds>(reduceingEnd - indexingEnd).count() << endl;
    vector<int> similar = lis::getSimilar(seqMin, indexTable);

    chrono::steady_clock::time_point lisEnd = chrono::steady_clock::now();
    cout << "lis end in: " << chrono::duration_cast<chrono::seconds>(lisEnd - reduceingEnd).count() << endl;
    cout << "total end in: " << chrono::duration_cast<chrono::seconds>(lisEnd - start).count() << endl;
    cout << similar.size() << endl;


    chrono::steady_clock::time_point opalStart = chrono::steady_clock::now();


    // OPAL -------------------------------------------------------------------------------

    int alphabetLength = 24;
    int gapOpen = 3;
    int gapExt = 1;
    int scoreMatrix[24 * 24] = {
         5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -2, -1, -1, -3, -1,  1,  0, -3, -2,  0, -2, -1, -1, -5,
        -2,  7, -1, -2, -4,  1,  0, -3,  0, -4, -3,  3, -2, -3, -3, -1, -1, -3, -1, -3, -1,  0, -1, -5,
        -1, -1,  7,  2, -2,  0,  0,  0,  1, -3, -4,  0, -2, -4, -2,  1,  0, -4, -2, -3,  5,  0, -1, -5,
        -2, -2,  2,  8, -4,  0,  2, -1, -1, -4, -4, -1, -4, -5, -1,  0, -1, -5, -3, -4,  6,  1, -1, -5,
        -1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -3, -3, -1, -5,
        -1,  1,  0,  0, -3,  7,  2, -2,  1, -3, -2,  2,  0, -4, -1,  0, -1, -1, -1, -3,  0,  4, -1, -5,
        -1,  0,  0,  2, -3,  2,  6, -3,  0, -4, -3,  1, -2, -3, -1, -1, -1, -3, -2, -3,  1,  5, -1, -5,
         0, -3,  0, -1, -3, -2, -3,  8, -2, -4, -4, -2, -3, -4, -2,  0, -2, -3, -3, -4, -1, -2, -1, -5,
        -2,  0,  1, -1, -3,  1,  0, -2, 10, -4, -3,  0, -1, -1, -2, -1, -2, -3,  2, -4,  0,  0, -1, -5,
        -1, -4, -3, -4, -2, -3, -4, -4, -4,  5,  2, -3,  2,  0, -3, -3, -1, -3, -1,  4, -4, -3, -1, -5,
        -2, -3, -4, -4, -2, -2, -3, -4, -3,  2,  5, -3,  3,  1, -4, -3, -1, -2, -1,  1, -4, -3, -1, -5,
        -1,  3,  0, -1, -3,  2,  1, -2,  0, -3, -3,  6, -2, -4, -1,  0, -1, -3, -2, -3,  0,  1, -1, -5,
        -1, -2, -2, -4, -2,  0, -2, -3, -1,  2,  3, -2,  7,  0, -3, -2, -1, -1,  0,  1, -3, -1, -1, -5,
        -3, -3, -4, -5, -2, -4, -3, -4, -1,  0,  1, -4,  0,  8, -4, -3, -2,  1,  4, -1, -4, -4, -1, -5,
        -1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3, -4, 10, -1, -1, -4, -3, -3, -2, -1, -1, -5,
         1, -1,  1,  0, -1,  0, -1,  0, -1, -3, -3,  0, -2, -3, -1,  5,  2, -4, -2, -2,  0,  0, -1, -5,
         0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  2,  5, -3, -2,  0,  0, -1, -1, -5,
        -3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1,  1, -4, -4, -3, 15,  2, -3, -5, -2, -1, -5,
        -2, -1, -2, -3, -3, -1, -2, -3,  2, -1, -1, -2,  0,  4, -3, -2, -2,  2,  8, -1, -3, -2, -1, -5,
         0, -3, -3, -4, -1, -3, -3, -4, -4,  4,  1, -3,  1, -1, -3, -2,  0, -3, -1,  5, -3, -3, -1, -5,
        -2, -1,  5,  6, -3,  0,  1, -1,  0, -4, -4,  0, -3, -4, -2,  0,  0, -5, -3, -3,  6,  1, -1, -5,
        -1,  0,  0,  1, -3,  4,  5, -2,  0, -3, -3,  1, -1, -4, -1,  0, -1, -2, -2, -3,  1,  5, -1, -5,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5,
        -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  1
    };

    unsigned char blosumAlphabet[24] = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
                                  'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*'};
    unsigned char toBlosumIndex[26];
    for (int i = 0; i < 26; i++)
        toBlosumIndex[i] = 23;
    for (int i = 0; i < 23; i++) {
        toBlosumIndex[blosumAlphabet[i] - 'A'] = i;

    }
// Query

    unsigned char* opalQuery =  new unsigned char[sizeof(unsigned char) * queryLength];
    for (int i = 0; i < queryLength; i++) {
        opalQuery[i] = toBlosumIndex[query[i] - 'A'];
    }

// Database
    int dbOpalLength = similar.size();

    unsigned char** opalDb = new unsigned char*[sizeof(unsigned char*) * dbOpalLength];
    int* opalDbSeqsLen = new int[sizeof(int) * dbOpalLength];
    int opdb = 0;
    for (int candidate: similar) {
        unsigned char* opalCandidate = new unsigned char[sizeof(unsigned char) * dbLengths[candidate]];
        for (int i = 0; i < dbLengths[candidate]; i++) {
            opalCandidate[i] = toBlosumIndex[db[candidate][i] - 'A'];
        }
        opalDb[opdb] = opalCandidate;
        opalDbSeqsLen[opdb] = dbLengths[candidate];
        opdb++;
    }

// Results for each sequence in database
    OpalSearchResult** results  = new OpalSearchResult*[sizeof(OpalSearchResult*) * dbOpalLength];
    for (int i = 0; i < dbOpalLength; i++) {
        results[i] = new OpalSearchResult();
        opalInitSearchResult(results[i]);
    }


// Do calculation!
    int resultCode = opalSearchDatabase(opalQuery, queryLength, opalDb, dbOpalLength, opalDbSeqsLen,
                                        gapOpen, gapExt, scoreMatrix, alphabetLength, results, OPAL_SEARCH_SCORE,
                                        OPAL_MODE_SW, OPAL_OVERFLOW_BUCKETS);


    chrono::steady_clock::time_point opalEnd = chrono::steady_clock::now();
    cout << "opal end in: " << chrono::duration_cast<chrono::seconds>(opalEnd - opalStart).count() << endl;
// Print scores
    int i = 0;
    for (int cand: similar) {
        printf("cand: %d, score: %d\n", cand, results[i]->score);
        i++;
    }


    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "total end in: " << chrono::duration_cast<chrono::seconds>(end - start).count() << endl;
    return 0;
}