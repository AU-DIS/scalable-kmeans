#include <string>
#include "dataset.cpp"
#include "kmeans_runner.cpp"
#include "strategies/lloyd_kmeans_strategy.cpp"
#include "benchmark/benchmark.h"


static void BM_Kmeans_lloyd(benchmark::State& state) { 
    std::string data_file_name = "/mnt/c/Users/kaspe/OneDrive/Skrivebord/Reps/scalable-kmeans/Data/steinn_14_jun/processed/misfit_Bi5d_full_h5_512_dct.txt";
    int n = 168;
    int d = 512*512;
    int k = 5;

    //; lloyd_strategy;
    KmeansRunner runner(std::make_unique<LloydKmeansStrategy>());
    
    runner.init_run(n,d,k,data_file_name);

    for (auto _ : state) { 
        runner.run_kmeans(n,d,k);
    }
}

int main(int argc, char **argv) { 
    //std::string data_file_name = argv[1];
    //int n = 168;
    //int d = 32*32;
    //int k = 5;

    //; lloyd_strategy;
    /*KmeansRunner runner(std::make_unique<LloydKmeansStrategy>());

    runner.init_run(n,d,k,data_file_name);
    runner.run_kmeans();
    runner.save_result(n,"/mnt/c/Users/kaspe/OneDrive/Skrivebord/Reps/scalable-kmeans/cpp/strategy_kmeans/result.txt");
    */
    BENCHMARK(BM_Kmeans_lloyd)->Iterations(20)->Unit(benchmark::kMillisecond);
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    //Dataset data(n, d);
    //data.load_datafile(data_file_name);

    //data.print_datasample();

    return 0;
};


