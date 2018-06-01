//
// Created by dario on 02.11.17..
//

#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <vector>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <vector>



namespace minimizer {
    typedef int hashType;

    struct Minimizer {
        hashType h;
        int position;

        Minimizer(hashType _h, int _position)
                : h(_h), position(_position) {}

        inline bool operator<(const Minimizer& other) const {
            return h < other.h;
        }

        inline bool operator>(const Minimizer& other) const {
            return h > other.h;
        }

    };

    struct Index {
        int sequenceIndex;
        int position;

        Index(int _sequenceIndex, int _position) : sequenceIndex(_sequenceIndex), position(_position) {}

        inline bool operator<(const Index& other) const {
            return sequenceIndex < other.sequenceIndex;
        }



    };

    typedef std::vector<std::vector<minimizer::Index>> IndexTable;

    // stavlja minimizere iz target-a u odgovarajuci vektoru u mapi hasheva
    void addMinimizers(const char* target, int targetLen, int targetIndex, int w, int k, 
                        IndexTable& indexTable);

    std::vector<Minimizer> computeForSequence(const char* target, int targetLen, int w, int k);


} // namespace minimizer

#endif //MINIMIZER_H
