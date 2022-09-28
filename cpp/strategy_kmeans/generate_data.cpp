#include "strategies/lloyd_kmeans_strategy.cpp"
#include <string>
//#include <iostream>
//#include "../utils/csv.hpp"

#ifdef DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 0
#endif
#define DEBUGPRINT(fmt, ...) \
            do { if (DEBUG_TEST) fprintf(stderr, fmt, __VA_ARGS__); } while (0)


int main(int argc, char **argv) { 
    //Handle args
    std::string data_file_name = argv[1];

    DEBUGPRINT("DATA PATH: %s", data_file_name.c_str());
    //TODO: Ground truth label comparison not supported for now
    

    //RUN DCT FILTERS
    //TODO: I have no expertise with with this so calling the python impl for now. Should be fine as it is mainly torch.
    std::string run_line = "../run_python_dct.sh ";
    run_line += data_file_name;
    system(run_line.c_str());


    //LOAD CSV DATA
    std::string csv_file = data_file_name;
    csv_file += "_dct.csv";







    
   
    //RUN KMEANS


}


