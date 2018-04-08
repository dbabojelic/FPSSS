#include <iostream>
#include <unordered_map>
#include <map>
#include "minimizer.h"
#include <algorithm>


extern "C" {
    #include "fasta.h"
}

int main() {
    //std::cout << "Hello, World!" << std::endl;
    //std::vector<minimizer::Minimizer> vec = minimizer::computeMinimizers("ACATACAGAGFGRGTHFGSA", 20, 10, 3);
    //for (minimizer::Minimizer mini: vec) {
     //   std::cout << mini.h << " " << mini.position << " " << std::endl;
    //}

    FASTAFILE *ffp;
    ffp = OpenFASTA("../uniprot_sprot.fasta");

    char *seq, *name;
    int len;
    int cnt = 0;
    std::unordered_map<minimizer::hashType,
            std::vector<minimizer::Index>> indexTable;

    while (ReadFASTA(ffp, &seq, &name, &len)) {
        minimizer::addMinimizers(seq, len, cnt, 10, 3, indexTable);
        free(seq);
        free(name);
        cnt++;
    }
    CloseFASTA(ffp);
    std::cout << cnt << std::endl;
    return 0;
}