#include <string>
#include "dataset.cpp"


int main(int argc, char **argv) { 
    std::string data_file_name = argv[1];
    int n = 168;
    int d = 64;

    Dataset data(n, d);
    data.load_datafile(data_file_name);

    data.print_datasample();

    return 0;
};