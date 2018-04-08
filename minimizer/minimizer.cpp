//
// Created by dario on 02.11.17..
//

#include <queue>
#include <assert.h>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include "minimizer.h"

namespace {
    std::unordered_map<minimizer::hashType, int> hashCnt;
    const int BASE = 31;

    int value(char c) {
        return c - 'A';
    }

    void push(const minimizer::Minimizer& triple, std::deque<minimizer::Minimizer>& dq) {
        while (!dq.empty() && triple < dq.back()) {
          dq.pop_back();
        }
        dq.push_back(triple);
    }

    void processState(std::deque<minimizer::Minimizer>& dq,
                      std::unordered_map<minimizer::hashType, std::vector<minimizer::Index>>& indexTable,
                      int targetIndex, int& lastPositionTaken) {
        assert(!dq.empty());
        minimizer::Minimizer& front = dq.front();
        dq.pop_front();

        if (lastPositionTaken < front.position) {
            indexTable[front.h].push_back(minimizer::Index(targetIndex,front));
            hashCnt[front.h]++;
            lastPositionTaken = front.position;
        }
        while (!dq.empty() && dq.front().h == front.h) {
            front = dq.front();
            dq.pop_front();
            if (lastPositionTaken < front.position) {
                indexTable[front.h].push_back(minimizer::Index(targetIndex,front));
                hashCnt[front.h]++;
                lastPositionTaken = front.position;
            }
        }
        dq.push_front(front);
    }

    void pop(int position, std::deque<minimizer::Minimizer>& dq) {
        while (!dq.empty() && dq.front().position == position)
            dq.pop_front();
    }


} // namespace

namespace minimizer {

    void addMinimizers(const char* target, int targetLen, int targetIndex, int w, int k, 
                       std::unordered_map<hashType, std::vector<Index>>& indexTable) {
        int n = targetLen;


        if (n < k) {
            return; //ne postoji niti jedan kmer
        }

        if (n < k + w - 1) {//ne postoji ni jedan window od w kmera
            w = n - k + 1; // smanji velicinu trazenog prozora na najvise sta mozes, da se nadje barem jedan minimizer
        }


        std::deque<Minimizer> dqMin,dqMax;

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
            Minimizer mp1 = Minimizer(tmpHash, i);
            push(mp1, dqMin);
            tmpHash -= lastPower * value(target[i]);
            tmpHash *= BASE;
            tmpHash += value(target[i + k]);
        }

        processState(dqMin, indexTable, targetIndex, lastPositionTaken);

        for (int i = w; i < n - k + 1; i++) {
            pop(i - w, dqMin);
            Minimizer mp1 = Minimizer(tmpHash, i);
            push(mp1, dqMin);
            processState(dqMin, indexTable, targetIndex, lastPositionTaken);

            tmpHash -= lastPower * value(target[i]);
            tmpHash *= BASE;
            tmpHash += value(target[i + k]);
        }
    }

    std::vector<minimizer::Minimizer> reduceMinimizers(std::vector<Minimizer>& minimizers) {
        std::vector<Minimizer> ret;
        for (Minimizer mini: minimizers) {
            if (hashCnt[mini.h] < 34)
                ret.push_back(mini);
        }
        return ret;
    }




} // namespace minimizer
