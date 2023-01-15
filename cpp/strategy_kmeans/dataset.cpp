#ifndef DATASET_
#define DATASET_

#include <array>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
using std::array;
using std::string;


class Dataset {
    private:
    // 
    int n, d;
    std::ifstream data_file;
    //Data as one memoryblock: data[i][j] written as data[i*d+j]
    double* data;
    
    public:
    explicit Dataset(int _n, int _d){
        n = _n;
        d = _d;
    
        data = new double[n*d];
    };
    explicit Dataset(){
        n = 0;
        d = 0;
    
        data = new double[0];
    };
    ~Dataset(){    
        delete data; 
    };

    
    
    public:
    double get(int i, int j) {
        return data[i*d+j];
    };

    void set(int i, int j, double val) {
        data[i*d+j] = val;
    };

    double* get_data_pointer() {
        return data;
    }
    

    void load_datafile(string data_file_name){
        std::ifstream data_file(data_file_name.c_str());
        int data_count = std::count(std::istreambuf_iterator<char>(data_file),
                                std::istreambuf_iterator<char>(), '\n');
        if (data_count < n) {
            std::cout << "NOT ENOUGH DATA!! data_count: " << data_count << " n: " << n << "\n";
            exit(3);
        }

        data_file.clear();
        data_file.seekg(0, std::ios::beg);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < d; j++) {
                data_file >> data[i*d+j];
            }
        }
    };

    void print_datasample() {
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                std::cout << data[i*d+j] << " ";
            }
            std::cout << "\n";
        }
        std::cout << std::endl;
    };
};
#endif