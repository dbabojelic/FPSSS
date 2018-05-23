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

    const int murphy[] = {A, OTHER, C, EDNQ, EDNQ, FYW, G, H, ILVM, OTHER, KR, ILVM, ILVM, EDNQ, OTHER, P, EDNQ, KR, ST, ST,
                            OTHER, ILVM, FYW, OTHER, FYW, OTHER};
    inline int value(char c) {
        return murphy[c - 'A'];
    }

    void push(minimizer::Minimizer triple, std::deque<minimizer::Minimizer>& dq) {
        while (!dq.empty() && triple < dq.back()) {
            dq.pop_back();
        }
        dq.push_back(triple);
    }

    void processState(std::deque<minimizer::Minimizer>& dq, minimizer::IndexTable& indexTable,
                      int targetIndex, int& lastPositionTaken) {
        minimizer::Minimizer front = dq.front();
        dq.pop_front();

        if (lastPositionTaken < front.position && front.h < indexTable.size()) {
            indexTable[front.h].push_back(minimizer::Index(targetIndex,front.position));
            lastPositionTaken = front.position;
        }
        while (!dq.empty() && dq.front().h == front.h) {
            front = dq.front();
            dq.pop_front();
            if (lastPositionTaken < front.position && front.h < indexTable.size()) {
                indexTable[front.h].push_back(minimizer::Index(targetIndex,front.position));
                lastPositionTaken = front.position;
            }
        }
        dq.push_front(front);
    }

    void pop(int position, std::deque<minimizer::Minimizer>& dq) {
        while (!dq.empty() && dq.front().position == position)
            dq.pop_front();
    }

    bool sortByPosition(const minimizer::Minimizer& m1, const minimizer::Minimizer& m2) {
        return m1.position < m2.position;
    }


} // namespace

namespace minimizer {
    

    void addMinimizers(const char* target, int targetLen, int targetIndex, int w, int k, 
                       IndexTable& indexTable) {
        int n = targetLen;


        if (n < k) {
            return; //ne postoji niti jedan kmer
        }

        if (n < k + w - 1) {//ne postoji ni jedan window od w kmera
            w = n - k + 1; // smanji velicinu trazenog prozora na najvise sta mozes, da se nadje barem jedan minimizer
        }
        int maxHashes = 1;
        for (int i = 0; i < k; i++)
            maxHashes *= BASE;

        while(indexTable.size() < maxHashes) {
            indexTable.push_back(vector<minimizer::Index>());
        }


        std::deque<Minimizer> dqMin;

        hashType lastPower = 1;
        hashType tmpHash = 0;
        hashType tmpPot = 1;
        int lastPositionTaken = -1;

        for (int i = 0; i < k - 1; i++)
            lastPower *= BASE;

        for (int i = 0; i < k; i++) {
            tmpHash *= BASE;
            tmpHash += value(target[i]);
            tmpPot *= BASE;
        }

        // queue s maksimumom algoritam
        for (int i = 0; i < w; i++) {
            Minimizer mp1(tmpHash, i);
            push(mp1, dqMin);
            tmpHash -= lastPower * value(target[i]);
            tmpHash *= BASE;
            tmpHash += value(target[i + k]);
        }

        processState(dqMin, indexTable, targetIndex, lastPositionTaken);

        for (int i = w; i < n - k + 1; i++) {
            pop(i - w, dqMin);
            Minimizer mp1(tmpHash, i);
            push(mp1, dqMin);
            processState(dqMin, indexTable, targetIndex, lastPositionTaken);

            tmpHash -= lastPower * value(target[i]);
            tmpHash *= BASE;
            tmpHash += value(target[i + k]);
        }
    }

    std::vector<Minimizer> computeForSequence(const char *target, int targetLen, int w, int k) {
        IndexTable table;
        addMinimizers(target, targetLen, -1, w, k, table);
        vector<Minimizer> ret;

        for (int i = 0; i < table.size(); i++) {
            for (Index el: table[i]) {
                ret.push_back(Minimizer(i, el.position));
            }
        }

        sort(ret.begin(), ret.end(), sortByPosition);
        return ret;
    }


} // namespace minimizer
