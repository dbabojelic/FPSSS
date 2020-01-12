//
// Created by dario on 12.12.17..
//

#ifndef PROJEKT_LIS_H
#define PROJEKT_LIS_H

#include <minimizer.h>
#include <unordered_map>


namespace filter {

    std::vector<int> getSimilar(std::vector<std::vector<minimizer::Minimizer>> v1,
                                                 minimizer::IndexTable* indexTable, int BANDS, int reduceTo);
} // namespace filter

#endif //PROJEKT_LIS_H
