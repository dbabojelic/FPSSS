#include <iostream>
#include <unordered_map>
#include <map>
#include <map>
#include "minimizer.h"
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

    char *seq, *name;
    int len;
    int cnt = 0;
    unordered_map<minimizer::hashType,
            vector<minimizer::Index>> indexTable;

    while (ReadFASTA(ffp, &seq, &name, &len)) {
        minimizer::addMinimizers(seq, len, cnt, 10, 3, indexTable);
        free(seq);
        free(name);
    }
    CloseFASTA(ffp);
    vector<pair<int, int>> hashToCnt;
    for (auto it = indexTable.begin(); it != indexTable.end(); it++) {

        hashToCnt.push_back({it->second.size(), it->first});
    }
    sort(hashToCnt.begin(), hashToCnt.end());
//    cout << hashToCnt.size() << endl;
    for (pair<int, int> p: hashToCnt) {
 //       cout << p.second << " " << p.first << endl;
    }

    cout << "----------------" << endl << endl;


    reduceIndexTable(indexTable, (1));
    hashToCnt.clear();
    for (auto it = indexTable.begin(); it != indexTable.end(); it++) {

        hashToCnt.push_back({it->second.size(), it->first});
    }
    sort(hashToCnt.begin(), hashToCnt.end());
  //  cout << hashToCnt.size() << endl;
    for (pair<int, int> p: hashToCnt) {
   //     cout << p.second << " " << p.first << endl;
    }
    return 0;
}