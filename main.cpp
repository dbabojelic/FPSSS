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

using namespace std;
typedef unordered_map<minimizer::hashType, vector<minimizer::Index>> IndexTable;

extern "C" {
    #include "fasta.h"
}

void reduceIndexTable(IndexTable& indexTable, double p) {
    assert(p <= 100 && p >= 0);
    int indexCnt = 0;
    for (auto it = indexTable.begin(); it != indexTable.end(); it++) {
        indexCnt += it->second.size();
    }
    int limit = (p/100)*indexCnt;

    for (auto it = indexTable.begin(); it != indexTable.end(); ) {
        if (it->second.size() > limit) {
            it = indexTable.erase(it);
        }
        else {
            it++;
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
    ffp = OpenFASTA("../uniprot_sprot.fasta");

    chrono::steady_clock::time_point start = chrono::steady_clock::now();

    char *seq, *name;
    int len;
    int cnt = 0;
    unordered_map<minimizer::hashType,
            vector<minimizer::Index>> indexTable;

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
    return 0;
}