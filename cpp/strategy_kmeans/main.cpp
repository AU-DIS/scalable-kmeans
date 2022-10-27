#include <string>
#include "dataset.cpp"
#include "kmeans_runner.cpp"
#include "strategies/lloyd_kmeans_strategy.cpp"
#include "strategies/marigold_kmeans_strategy.cpp"
#include "strategies/stepwise_kmeans_strategy.cpp"
#include "strategies/elkham_kmeans_strategy.cpp"
#include "benchmark/benchmark.h"

std::string data_file_name;
int n;
int d;
int k;

static void BM_Kmeans_lloyd(benchmark::State& state) { 
    //lloyd_strategy;
    KmeansRunner runner(std::make_unique<LloydKmeansStrategy>());

    runner.init_run(n,d,k,data_file_name);

    for (auto _ : state) { 
        runner.run_kmeans(n,d,k);
        
    }
}

static void BM_Kmeans_MARIGOLD(benchmark::State& state) { 

    //MARIGOLD_strategy;
    KmeansRunner runner(std::make_unique<MARIGOLDKmeansStrategy>());
    
    runner.init_run(n,d,k,data_file_name);

    for (auto _ : state) { 
        runner.run_kmeans(n,d,k);
    }
}

static void BM_Kmeans_StepWise(benchmark::State& state) { 


    //StepWise_strategy;
    KmeansRunner runner(std::make_unique<StepWiseKmeansStrategy>());
    
    runner.init_run(n,d,k,data_file_name);

    for (auto _ : state) { 
        runner.run_kmeans(n,d,k);
    }
}

static void BM_Kmeans_ElkHam(benchmark::State& state) { 
    
    //ElkHam_strategy;
    KmeansRunner runner(std::make_unique<ElkHamKmeansStrategy>());
    
    runner.init_run(n,d,k,data_file_name);

    for (auto _ : state) { 
        runner.run_kmeans(n,d,k);
    }
}
 




int main(int argc, char **argv) { 
    data_file_name = "/mnt/c/Users/kaspe/OneDrive/Skrivebord/Reps/scalable-kmeans/Data/steinn_14_jun/processed/gr_flake_256_h5_dct.txt";
    n = 1911;
    d = 256*256;
    k = 40;

    //; lloyd_strategy;
    /*KmeansRunner runner(std::make_unique<LloydKmeansStrategy>());

    runner.init_run(n,d,k,data_file_name);
    runner.run_kmeans();
    runner.save_result(n,"/mnt/c/Users/kaspe/OneDrive/Skrivebord/Reps/scalable-kmeans/cpp/strategy_kmeans/result.txt");
    */
    //BENCHMARK(BM_Kmeans_lloyd)->Iterations(1)->Unit(benchmark::kMillisecond);
    BENCHMARK(BM_Kmeans_MARIGOLD)->Iterations(1)->Unit(benchmark::kMillisecond);
    BENCHMARK(BM_Kmeans_StepWise)->Iterations(1)->Unit(benchmark::kMillisecond);
    BENCHMARK(BM_Kmeans_ElkHam)->Iterations(1)->Unit(benchmark::kMillisecond);
    
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    //Dataset data(n, d);
    //data.load_datafile(data_file_name);

    //data.print_datasample();

    return 0;
};


