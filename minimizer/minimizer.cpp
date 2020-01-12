//
// Created by dario on 02.11.17..
//

#include <queue>
#include <assert.h>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include <memory>
#include <vector>
#include "minimizer.h"

using namespace std;


namespace {
    const int BASE = 11;
    const int A = 0;
    const int KR = 1;
    const int EDNQ = 2;
    const int C = 3;
    const int G = 4;
    const int H = 5;
    const int ILVM = 6;
    const int FYW = 7;
    const int P = 8;
    const int ST = 9;
    const int OTHER = 10000;

    const int murphy[] = {A, OTHER, C, EDNQ, EDNQ, FYW, G, H, ILVM, OTHER, KR, ILVM, ILVM, EDNQ, OTHER, P, EDNQ, KR, ST,
                          ST,
                          OTHER, ILVM, FYW, OTHER, FYW, OTHER};

    inline int value(char c) {
        return murphy[c - 'A'];
    }

    void push(minimizer::Minimizer triple, std::deque<minimizer::Minimizer> &dq) {
        while (!dq.empty() && triple < dq.back()) {
            dq.pop_back();
        }
        dq.push_back(triple);
    }

    void insertInIndexTables(minimizer::IndexTable *indexTables, const int BANDS, int k, int targetIndex, int position,
                             const char *target) {
        assert(k % BANDS == 0);
        const int BAND_SIZE = k / BANDS;
        for (int band = 0; band < BANDS; band++) {
            minimizer::hashType tmp_hash = 0;

            int start = position + band * BAND_SIZE;
            int end = position + band * BAND_SIZE + BAND_SIZE;
            for (int pos = start; pos < end; pos++) {
                tmp_hash *= BASE;
                tmp_hash += value(target[pos]);
            }

            if (tmp_hash < indexTables[band].size())
                indexTables[band][tmp_hash].push_back(minimizer::Index(targetIndex, position));
        }

    }

    void processState(std::deque<minimizer::Minimizer> &dq, minimizer::IndexTable *indexTable, const int BANDS,
                      int targetIndex, int &lastPositionTaken, int k, const char *target) {
        minimizer::Minimizer front = dq.front();
        dq.pop_front();

        if (lastPositionTaken < front.position) {
//            indexTable[front.h].push_back(minimizer::Index(targetIndex, front.position));
            insertInIndexTables(indexTable, BANDS, k, targetIndex, front.position, target);
            lastPositionTaken = front.position;
        }
        while (!dq.empty() && dq.front().h == front.h) {
            front = dq.front();
            dq.pop_front();
            if (lastPositionTaken < front.position) {
//                indexTable[front.h].push_back(minimizer::Index(targetIndex, front.position));
                insertInIndexTables(indexTable, BANDS, k, targetIndex, front.position, target);
                lastPositionTaken = front.position;
            }
        }
        dq.push_front(front);
    }

    void pop(int position, std::deque<minimizer::Minimizer> &dq) {
        while (!dq.empty() && dq.front().position == position)
            dq.pop_front();
    }

    bool sortByPosition(const minimizer::Minimizer &m1, const minimizer::Minimizer &m2) {
        return m1.position < m2.position;
    }

    minimizer::hashType modMul(minimizer::hashType a, minimizer::hashType b) {
        return (a * 1LL * b) % minimizer::MOD;
    }

    minimizer::hashType modAdd(minimizer::hashType a, minimizer::hashType b) {
        minimizer::hashType tmp = a + b;
        if (tmp > minimizer::MOD)
            tmp -= minimizer::MOD;
        return tmp;
    }

    minimizer::hashType modSub(minimizer::hashType a, minimizer::hashType b) {
        if (a < b)
            a += minimizer::MOD;
        return a - b;
    }


} // namespace

namespace minimizer {


    void addMinimizers(const char *target, int targetLen, int targetIndex, int w, int k,
                       IndexTable *indexTable, const int BANDS) {
        int n = targetLen;

        assert(k % BANDS == 0);
        const int BAND_SIZE = k / BANDS;


        if (n < k) {
            return; //ne postoji niti jedan kmer
        }

        if (n < k + w - 1) {//ne postoji ni jedan window od w kmera
            w = n - k + 1; // smanji velicinu trazenog prozora na najvise sta mozes, da se nadje barem jedan minimizer
        }
        hashType maxHashes = 1;
        for (int i = 0; i < BAND_SIZE; i++)
            maxHashes *= BASE;

        for (int i = 0; i < BANDS; i++) {
            while (indexTable[i].size() < maxHashes) {
                indexTable[i].push_back(vector<minimizer::Index>());
            }
        }

//        while (indexTable.size() < maxHashes) {
//            indexTable.push_back(vector<minimizer::Index>());
//        }


        std::deque<Minimizer> dqMin;

        hashType lastPower = 1;
        hashType tmpHash = 0;
        hashType tmpPot = 1;
        int lastPositionTaken = -1;

        for (int i = 0; i < k - 1; i++) {
            lastPower = modMul(lastPower, BASE);
        }

        for (int i = 0; i < k; i++) {
            tmpHash = modMul(tmpHash, BASE);
            tmpHash = modAdd(tmpHash, value(target[i]));

            tmpPot = modMul(tmpPot, BASE);
        }

        // queue s maksimumom algoritam
        for (int i = 0; i < w; i++) {
            Minimizer mp1(tmpHash, i);
            push(mp1, dqMin);
            tmpHash = modSub(tmpHash, modMul(lastPower, value(target[i])));
            tmpHash = modMul(tmpHash, BASE);
            tmpHash = modAdd(tmpHash, value(target[i + k]));
        }

        processState(dqMin, indexTable, BANDS, targetIndex, lastPositionTaken, k, target);

        for (int i = w; i < n - k + 1; i++) {
            pop(i - w, dqMin);
            Minimizer mp1(tmpHash, i);
            push(mp1, dqMin);
            processState(dqMin, indexTable, BANDS, targetIndex, lastPositionTaken, k, target);

            tmpHash = modSub(tmpHash, modMul(lastPower, value(target[i])));
            tmpHash = modMul(tmpHash, BASE);
            tmpHash = modAdd(tmpHash, value(target[i + k]));
        }
    }

    std::vector<std::vector<Minimizer>>
    computeForSequence(const char *target, int targetLen, int w, int k, const int BANDS) {
        IndexTable* table = new IndexTable[BANDS];
        addMinimizers(target, targetLen, -1, w, k, table, BANDS);
        vector<vector<Minimizer>> ret (BANDS, vector<Minimizer>());


        for (int band = 0; band < BANDS; band++) {
            for (int i = 0; i < table[band].size(); i++) {
                for (Index el: table[band][i]) {
                    ret[band].push_back(Minimizer(i, el.position));
                }
            }
            sort(ret[band].begin(), ret[band].end(), sortByPosition);
        }

        return ret;
    }


} // namespace minimizer
