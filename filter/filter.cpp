//
// Created by dario on 12.12.17..
//

#include "filter.h"
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <vector>
using namespace std;

namespace filter {

    std::vector<int> getSimilar(std::vector<std::vector<minimizer::Minimizer>> v1, minimizer::IndexTable *indexTable, int BANDS, int reduceTo) {
        unordered_map<int, vector<int>> NPRs;

        for (int band = 0; band < BANDS; band++) {
            for (minimizer::Minimizer mini: v1[band]) {
                int prosli = -1;
                for (auto &element: indexTable[band][mini.h]) {
                    prosli = element.sequenceIndex;
                    NPRs[element.sequenceIndex].push_back(element.position - mini.position);
                }
            }
        }

        vector<int> ret;
        vector<pair<int, int>> candidates;
        int sz = v1.size();
        for (auto& npr: NPRs) {
            sort(npr.second.begin(), npr.second.end());
            int maxDiff = 50;
            int p1 = 0;
            int p2 = 0;
            int sz = npr.second.size();
            int maxSegment = 0;
            while (p2 < sz) {
                while (npr.second[p2] - npr.second[p1] > maxDiff)
                    p1++;
                maxSegment = max(maxSegment, p2 - p1 + 1);
                p2++;
            }
            candidates.push_back({maxSegment, npr.first});
        }

        sort(candidates.begin(), candidates.end());
        for (int i = candidates.size() - 1; i >= 0 && reduceTo > 0; i--, reduceTo--) {
            ret.push_back(candidates[i].second);
        }
        return ret;
    }
}
