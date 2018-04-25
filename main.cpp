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
    //cout << "Hello, World!" << endl;
    //vector<minimizer::Minimizer> vec = minimizer::computeMinimizers("ACATACAGAGFGRGTHFGSA", 20, 10, 3);
    //for (minimizer::Minimizer mini: vec) {
     //   cout << mini.h << " " << mini.position << " " << endl;
    //}

    FASTAFILE *ffp;
    ffp = OpenFASTA("../uniprot_sprot3MB.fasta");

    chrono::steady_clock::time_point start = chrono::steady_clock::now();

    char *seq, *name;
    int len;
    int cnt = 0;

    minimizer::IndexTable indexTable;

    vector<minimizer::Minimizer> seqMin;
    while (ReadFASTA(ffp, &seq, &name, &len)) {
        minimizer::addMinimizers(seq, len, cnt, 10, 3, indexTable);
        if (cnt == 7) {
            seqMin = minimizer::computeForSequence(seq, len, 10, 3);
        }
        cnt++;
        free(seq);
        free(name);
    }

    CloseFASTA(ffp);

    chrono::steady_clock::time_point indexingEnd = chrono::steady_clock::now();
    cout << "indexing finished in: " << chrono::duration_cast<chrono::seconds>(indexingEnd - start).count() << endl;
//    vector<pair<int, int>> hashToCnt;
 //   for (auto it = indexTable.begin(); it != indexTable.end(); it++) {

  //      hashToCnt.push_back({it->second.size(), it->first});
   // }
    //sort(hashToCnt.begin(), hashToCnt.end());
//    cout << hashToCnt.size() << endl;
    //for (pair<int, int> p: hashToCnt) {
 //       cout << p.second << " " << p.first << endl;
   // }

    //cout << "----------------" << endl << endl;


    reduceIndexTable(indexTable, (1));
    chrono::steady_clock::time_point reduceingEnd = chrono::steady_clock::now();
    cout << "reduceing end in: " << chrono::duration_cast<chrono::seconds>(reduceingEnd - indexingEnd).count() << endl;
    //hashToCnt.clear();
    //for (auto it = indexTable.begin(); it != indexTable.end(); it++) {

     //   hashToCnt.push_back({it->second.size(), it->first});
   // }
    //sort(hashToCnt.begin(), hashToCnt.end());
  //  cout << hashToCnt.size() << endl;
    //for (pair<int, int> p: hashToCnt) {
   //     cout << p.second << " " << p.first << endl;
   // }
    cout << "ukupno u bazi proteina ima: " << cnt << endl;
    vector<int> ret = lis::getSimilar(seqMin, indexTable);
    chrono::steady_clock::time_point lisEnd = chrono::steady_clock::now();
    cout << "lis end in: " << chrono::duration_cast<chrono::seconds>(lisEnd - reduceingEnd).count() << endl;
    cout << "total end in: " << chrono::duration_cast<chrono::seconds>(lisEnd - start).count() << endl;
    cout << ret.size() << endl;
       // cout << i << endl;



    int alphabetLength = 4;
    int gapOpen = 3;
    int gapExt = 1;
    int scoreMatrix[16] = {
            2, -1, -3, 0,
            -1, 4, -5, -1,
            -3, -5, 1, -10,
            0, -1, -10, 4
    };

// Query
    int queryLength = 10;
    unsigned char query[10] = {0,1,3,2,1,0,3,0,1,1};

// Database
    int dbLength = 4;
    unsigned char dbSeq1[14] = {1,3,2,3,0,0,1,0,2,2,1,2,3,2};
    unsigned char dbSeq2[12] = {2,1,1,3,2,0,0,2,2,0,2,1};
    unsigned char dbSeq3[13] = {0,0,2,1,0,3,1,1,2,3,2,1,0};
    unsigned char dbSeq4[9] = {2,3,3,3,1,1,2,2,0};
    unsigned char* db[4] = { dbSeq1, dbSeq2, dbSeq3, dbSeq4 };
    int dbSeqsLengths[4] = {14, 12, 13, 9};

    cout << "here maybe" << endl;
// Results for each sequence in database
    OpalSearchResult* results[4];
    for (int i = 0; i < 4; i++) {
        results[i] = new OpalSearchResult();
        opalInitSearchResult(results[i]);
    }

    cout << "okay here" << endl;
// Do calculation!
    int resultCode = opalSearchDatabase(query, queryLength, db, dbLength, dbSeqsLengths,
                                        gapOpen, gapExt, scoreMatrix, alphabetLength, results, OPAL_SEARCH_SCORE,
                                        OPAL_MODE_SW, OPAL_OVERFLOW_BUCKETS);

    cout << "here to" << endl;
// Print scores
    printf("%d %d %d %d\n", results[0]->score, results[1]->score, results[2]->score, results[3]->score);
    return 0;
}