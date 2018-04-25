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

typedef unordered_map<minimizer::hashType, std::vector<minimizer::Index>> IndexTable;

namespace {
    const int BASE = 31;
    std::unordered_map<char, int> murphy = {{'A', 0},
                                            {'K', 1}, {'R', 1},
                                            {'E', 2}, {'D', 2}, {'N', 2}, {'Q', 2},
                                            {'C', 3},
                                            {'G', 4},
                                            {'H', 5},
                                            {'I', 6}, {'L', 6}, {'V', 6}, {'M', 6},
                                            {'F', 7}, {'Y', 7}, {'W', 7},
                                            {'P', 8},
                                            {'S', 9}, {'T', 9}};

    int value(char c) {
        return murphy[c];
    }

    void push(shared_ptr<minimizer::Minimizer> triple, std::deque<shared_ptr<minimizer::Minimizer>>& dq) {
        while (!dq.empty() && *triple < *dq.back()) {
            dq.pop_back();
        }
        dq.push_back(triple);
    }

    void processState(std::deque<shared_ptr<minimizer::Minimizer>>& dq, IndexTable& indexTable,
                      int targetIndex, int& lastPositionTaken) {
        assert(!dq.empty());
        shared_ptr<minimizer::Minimizer> front = dq.front();
        dq.pop_front();

        if (lastPositionTaken < front->position) {
            indexTable[front->h].push_back(minimizer::Index(targetIndex,front->position));
            lastPositionTaken = front->position;
        }
        while (!dq.empty() && dq.front()->h == front->h) {
            front = dq.front();
            dq.pop_front();
            if (lastPositionTaken < front->position) {
                indexTable[front->h].push_back(minimizer::Index(targetIndex,front->position));
                lastPositionTaken = front->position;
            }
        }
        dq.push_front(front);
    }

    void pop(int position, std::deque<shared_ptr<minimizer::Minimizer>>& dq) {
        while (!dq.empty() && dq.front()->position == position)
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


        std::deque<shared_ptr<Minimizer>> dqMin,dqMax;

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
            shared_ptr<Minimizer> mp1 = make_shared<Minimizer>(tmpHash, i);
            push(mp1, dqMin);
            tmpHash -= lastPower * value(target[i]);
            tmpHash *= BASE;
            tmpHash += value(target[i + k]);
        }

        processState(dqMin, indexTable, targetIndex, lastPositionTaken);

        for (int i = w; i < n - k + 1; i++) {
            pop(i - w, dqMin);
            shared_ptr<Minimizer> mp1 = make_shared<Minimizer>(tmpHash, i);
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

        for (auto it = table.begin(); it != table.end(); it++) {
            for (Index el: it->second) {
                ret.push_back(Minimizer(it->first, el.position));
            }
        }

        sort(ret.begin(), ret.end(), sortByPosition);
        return ret;
    }


} // namespace minimizer
