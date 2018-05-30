//
// Created by dario on 12.12.17..
//

#include "lis.h"
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <vector>
using namespace std;

namespace {
    int ceilIndex(std::vector<int> &v, int l, int r, int key) {
        while (r-l > 1) {
            int m = l + (r-l)/2;
            if (v[m] >= key)
                r = m;
            else
                l = m;
        }

        return r;
    }

    int lisF(std::vector<int>& v) {
        int n = (int)v.size();
        if (n == 0)
            return 0;

        std::vector<int> tail(n, 0);
        int length = 1; // always points empty slot in tail

        tail[0] = v[0];
        for (int i = 1; i < n; i++) {
            if (v[i] < tail[0])
                // new smallest value
                tail[0] = v[i];
            else if (v[i] > tail[length-1])
                // v[i] extends largest subsequence
                tail[length++] = v[i];
            else
                // v[i] will become end candidate of an existing subsequence or
                // Throw away larger elements in all LIS, to make room for upcoming grater elements than v[i]
                // (and also, v[i] would have already appeared in one of LIS, identify the location and replace it)
                tail[ceilIndex(tail, -1, length-1, v[i])] = v[i];
        }

        return length;
    }
}

namespace lis {

    std::vector<int> getSimilar(std::vector<minimizer::Minimizer> v1, minimizer::IndexTable &indexTable,
                                std::vector<int>& lens, int qLen, int cnt) {
        clock_t start = clock();
        unordered_map<int, vector<int>> seqsForLis;

        clock_t lisStart = clock();
        for (minimizer::Minimizer mini: v1) {
            int prosli = -1;
            for (auto& element: indexTable[mini.h]) {
//                if (element.sequenceIndex == prosli)
//                    continue;
//                if (std::abs(element.position - mini.position) > 0.2 * v1.size())
//                    continue;
                prosli = element.sequenceIndex;
                seqsForLis[element.sequenceIndex].push_back(element.position - mini.position);
            }
        }

        clock_t lisEnd = clock();
        vector<int> ret;
        vector<pair<int, int>> candidates;
        int sz = v1.size();
        for (auto& seq: seqsForLis) {
//            int lisValue = lisF(seq.second);
//            if (lisValue < 1)
//                continue;
//            double lisValueDenom = std::max(qLen, lens[seq.first]) * 1. / std::min(qLen, lens[seq.first]);
            sort(seq.second.begin(), seq.second.end());
            int diff = 50;
            int p1 = 0;
            int p2 = 0;
            int sz = seq.second.size();
            int lisValue = 0;
            while (p2 < sz) {
                while (seq.second[p2] - seq.second[p1] > diff)
                    p1++;
                lisValue = max(lisValue, p2 - p1 + 1);
                p2++;
            }
            candidates.push_back({lisValue, seq.first});
        }
        sort(candidates.begin(), candidates.end());
        for (int i = candidates.size() - 1; i >= 0 && cnt > 0; i--, cnt--) {
            ret.push_back(candidates[i].second);
        }
        clock_t end = clock();
        return ret;
    }
}
