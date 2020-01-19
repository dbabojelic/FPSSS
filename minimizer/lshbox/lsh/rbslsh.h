//////////////////////////////////////////////////////////////////////////////
/// Copyright (C) 2014 Gefu Tang <tanggefu@gmail.com>. All Rights Reserved.
///
/// This file is part of LSHBOX.
///
/// LSHBOX is free software: you can redistribute it and/or modify it under
/// the terms of the GNU General Public License as published by the Free
/// Software Foundation, either version 3 of the License, or(at your option)
/// any later version.
///
/// LSHBOX is distributed in the hope that it will be useful, but WITHOUT
/// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
/// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
/// more details.
///
/// You should have received a copy of the GNU General Public License along
/// with LSHBOX. If not, see <http://www.gnu.org/licenses/>.
///
/// @version 0.1
/// @author Gefu Tang & Zhifeng Xiao
/// @date 2014.6.30
//////////////////////////////////////////////////////////////////////////////

/**
 * @file rbslsh.h
 *
 * @brief Locality-Sensitive Hashing Scheme Based on Random Bits Sampling.
 */
#ifndef RBSLSH_H
#define RBSLSH_H
#include <map>
#include <vector>
#include <random>
#include <iostream>
#include <functional>
#include <ctime>
#include <algorithm>
namespace lshbox
{
/**
 * Locality-Sensitive Hashing Scheme Based on Random Bits Sampling.
 *
 *
 * For more information on random bits sampling based LSH, see the following reference.
 *
 *     P. Indyk and R. Motwani. Approximate Nearest Neighbor - Towards Removing
 *     the Curse of Dimensionality. In Proceedings of the 30th Symposium on Theory
 *     of Computing, 1998, pp. 604-613.
 *
 *     A. Gionis, P. Indyk, and R. Motwani. Similarity search in high dimensions
 *     via hashing. Proceedings of the 25th International Conference on Very Large
 *     Data Bases (VLDB), 1999.
 */
class rbsLsh
{
public:
    struct Parameter
    {
        /// Hash table size
        unsigned M;
        /// Number of hash tables
        unsigned L;
        /// Dimension of the vector, it can be obtained from the instance of Matrix
        unsigned D;
        /// Binary code bytes
        unsigned N;
        /// The Difference between upper and lower bound of each dimension
        unsigned C;
    };
    rbsLsh() {}
    rbsLsh(const Parameter &param_)
    {
        reset(param_);
    }
    ~rbsLsh() {}
    /**
     * Reset the parameter setting
     *
     * @param param_ A instance of rbsLsh::Parametor, which contains the necessary
     * parameters
     */
    inline void reset(const Parameter &param_);
    /**
     * get the hash value of a vector.
     *
     * @param k     The idx of the table
     * @param domin The pointer to the vector
     * @return      The hash value
     */
    inline unsigned getHashVal(unsigned k, const unsigned *domin);

    Parameter getParameters() {
        return param;
    }
private:
    Parameter param;
    std::vector<std::vector<unsigned> > rndBits;
    std::vector<std::vector<unsigned> > rndArray;
};
}

// ------------------------- implementation -------------------------
inline void lshbox::rbsLsh::reset(const Parameter &param_)
{
    param = param_;
    rndBits.resize(param.L);
    rndArray.resize(param.L);
    std::mt19937 rng(unsigned(time(NULL)));
    std::uniform_int_distribution<unsigned> usBits(0, param.D * param.C - 1);
    for (std::vector<std::vector<unsigned> >::iterator iter = rndBits.begin(); iter != rndBits.end(); ++iter)
    {
        while (iter->size() != param.N)
        {
            unsigned target = usBits(rng);
            if (std::find(iter->begin(), iter->end(), target) == iter->end())
            {
                iter->push_back(target);
            }
        }
        std::sort(iter->begin(), iter->end());
    }
    std::uniform_int_distribution<unsigned> usArray(0, param.M - 1);
    for (std::vector<std::vector<unsigned> >::iterator iter = rndArray.begin(); iter != rndArray.end(); ++iter)
    {
        for (unsigned i = 0; i != param.N; ++i)
        {
            iter->push_back(usArray(rng));
        }
    }
}

inline unsigned lshbox::rbsLsh::getHashVal(unsigned k, const unsigned *domin)
{
    unsigned sum(0), seq(0);
    for (std::vector<unsigned>::iterator it = rndBits[k].begin(); it != rndBits[k].end(); ++it)
    {
        if ((*it % param.C) <= unsigned(domin[*it / param.C]))
        {
            sum += rndArray[k][seq];
        }
        ++seq;
    }
    unsigned hashVal = sum % param.M;
    return hashVal;
}

#endif