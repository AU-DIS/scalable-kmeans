#include <string>
#include "dataset.cpp"
#include "kmeans_runner.cpp"
#include "strategies/lloyd_kmeans_strategy.cpp"
#include "strategies/marigold_kmeans_strategy.cpp"
#include "strategies/stepwise_kmeans_strategy.cpp"
#include "strategies/elkham_kmeans_strategy.cpp"
#include "strategies/elkan_kmeans_strategy.cpp"
#include "strategies/hamerly_kmeans_strategy.cpp"
#include "strategies/gstar_kmeans_strategy.cpp"
#include "benchmark/benchmark.h"

std::string data_file_name;
int n;
int d;
int k;
KmeansRunner runner(std::make_unique<LloydKmeansStrategy>());

static void BM_Kmeans_lloyd(benchmark::State& state) { 
    //lloyd_strategy;
    runner.set_strategy(std::make_unique<LloydKmeansStrategy>());


    for (auto _ : state) { 
        runner.run_kmeans(state.range(0),state.range(1),state.range(2));
    
        
    }
}

static void BM_Kmeans_gstar(benchmark::State& state) { 
    //lloyd_strategy;
    runner.set_strategy(std::make_unique<GstarKmeansStrategy>());

    for (auto _ : state) { 
        runner.run_kmeans(state.range(0),state.range(1),state.range(2));
    
        
    }
}

static void BM_Kmeans_MARIGOLD(benchmark::State& state) { 

    //MARIGOLD_strategy;
    //KmeansRunner runner(std::make_unique<MARIGOLDKmeansStrategy>());
    runner.set_strategy(std::make_unique<MARIGOLDKmeansStrategy>());
    for (auto _ : state) { 
        runner.run_kmeans(state.range(0),state.range(1),state.range(2));
    }
}

static void BM_Kmeans_StepWise(benchmark::State& state) { 


    //StepWise_strategy;
    runner.set_strategy(std::make_unique<StepWiseKmeansStrategy>());
    

    for (auto _ : state) { 
        runner.run_kmeans(state.range(0),state.range(1),state.range(2));
    }
}

static void BM_Kmeans_ElkHam(benchmark::State& state) { 
    

    runner.set_strategy(std::make_unique<ElkHamKmeansStrategy>());
    

    for (auto _ : state) { 
        runner.run_kmeans(state.range(0),state.range(1),state.range(2));
    }
}
 
static void BM_Kmeans_Elkan(benchmark::State& state) { 

    runner.set_strategy(std::make_unique<ElkanKmeansStrategy>());
    

    for (auto _ : state) { 
        runner.run_kmeans(state.range(0),state.range(1),state.range(2));
    }
}

static void BM_Kmeans_Hamerly(benchmark::State& state) { 
    

    runner.set_strategy(std::make_unique<HamerlyKmeansStrategy>());

    for (auto _ : state) { 
        runner.run_kmeans(state.range(0),state.range(1),state.range(2));
    }
}

static void BM_load_data_flake(benchmark::State& state) { 
    data_file_name = "/mnt/c/Users/kaspe/OneDrive/Skrivebord/Reps/scalable-kmeans/Data/steinn_14_jun/processed/gr_flake_"+std::to_string(int(sqrt(state.range(1))))+"_h5_dct.txt";
    //std::cout << data_file_name;
    for (auto _ : state) { 
        runner.init_run(state.range(0),state.range(1),state.range(2),data_file_name);
    }
}

static void BM_load_data_Bi5d(benchmark::State& state) { 
    data_file_name = "/mnt/c/Users/kaspe/OneDrive/Skrivebord/Reps/scalable-kmeans/Data/steinn_14_jun/processed/misfit_Bi5d_"+std::to_string(int(sqrt(state.range(1))))+"_h5_dct.txt";
    //std::cout << data_file_name;
    for (auto _ : state) { 
        runner.init_run(state.range(0),state.range(1),state.range(2),data_file_name);
    }
}

static void BM_load_data_Se3d(benchmark::State& state) { 
    data_file_name = "/mnt/c/Users/kaspe/OneDrive/Skrivebord/Reps/scalable-kmeans/Data/steinn_14_jun/processed/misfit_Se3d_"+std::to_string(int(sqrt(state.range(1))))+"_h5_dct.txt";
    //std::cout << data_file_name;
    for (auto _ : state) { 
        runner.init_run(state.range(0),state.range(1),state.range(2),data_file_name);
    }
}

static void BM_load_data_VB(benchmark::State& state) { 
    data_file_name = "/mnt/c/Users/kaspe/OneDrive/Skrivebord/Reps/scalable-kmeans/Data/steinn_14_jun/processed/misfit_VB_"+std::to_string(int(sqrt(state.range(1))))+"_h5_dct.txt";
    //std::cout << data_file_name;
    for (auto _ : state) { 
        runner.init_run(state.range(0),state.range(1),state.range(2),data_file_name);
    }
}

static void BM_load_data_Flickr(benchmark::State& state) { 
    data_file_name = "/mnt/c/Users/kaspe/OneDrive/Skrivebord/Reps/scalable-kmeans/Data/steinn_14_jun/processed/my_file_"+std::to_string(int(sqrt(state.range(1))))+"_dct.txt";
    //std::cout << data_file_name;
    for (auto _ : state) { 
        runner.init_run(state.range(0),state.range(1),state.range(2),data_file_name);
    }
}

int main(int argc, char **argv) { 
    n = 168;
    d = 128*128;
    k = 5;

    //; lloyd_strategy;
    /*KmeansRunner runner(std::make_unique<LloydKmeansStrategy>());

    runner.init_run(n,d,k,data_file_name);
    runner.run_kmeans();
    runner.save_result(n,"/mnt/c/Users/kaspe/OneDrive/Skrivebord/Reps/scalable-kmeans/cpp/strategy_kmeans/result.txt");
    */
    
    //BENCHMARK(BM_load_data_flake)->Args({n,d,k})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);
    BENCHMARK(BM_load_data_VB)->Args({n,d,k})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);
    //BENCHMARK(BM_load_data_Flickr)->Args({n,d,k})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);
    
    BENCHMARK(BM_Kmeans_lloyd)->Args({n,d,5})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_gstar)->Args({n,d,5})->Unit(benchmark::kMillisecond)->Iterations(1);
    
    //BENCHMARK(BM_Kmeans_Hamerly)->Args({n,d,5})->Args({n,d,40})->Unit(benchmark::kMillisecond)->Iterations(1);
    
    
    //BENCHMARK(BM_Kmeans_Elkan)->Args({n,d,5})->Args({n,d,40})->Unit(benchmark::kMillisecond)->Iterations(1);
    //BENCHMARK(BM_Kmeans_StepWise)->Args({n,d,5})->Args({n,d,40})->Unit(benchmark::kMillisecond)->Iterations(1);
    //BENCHMARK(BM_Kmeans_ElkHam)->Args({n,d,5})->Args({n,d,40})->Unit(benchmark::kMillisecond)->Iterations(1);
    
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({n,d,5})->Args({n,d,40})->Unit(benchmark::kMillisecond)->Iterations(1);
    


    //BENCHMARK(BM_Kmeans_lloyd)->Args({n,d,5})->Args({n,d,10})->Args({n,d,15})->Args({n,d,20})->Args({n,d,25})->Args({n,d,30})->Args({n,d,35})->Args({n,d,40})->Unit(benchmark::kMillisecond)->Iterations(1);
    //BENCHMARK(BM_Kmeans_Elkan)->Args({n,d,5})->Args({n,d,10})->Args({n,d,15})->Args({n,d,20})->Args({n,d,25})->Args({n,d,30})->Args({n,d,35})->Args({n,d,40})->Unit(benchmark::kMillisecond)->Iterations(1);
    
    //BENCHMARK(BM_Kmeans_StepWise)->Args({n,d,5})->Args({n,d,10})->Args({n,d,15})->Args({n,d,20})->Args({n,d,25})->Args({n,d,30})->Args({n,d,35})->Args({n,d,40})->Unit(benchmark::kMillisecond)->Iterations(1);
    //BENCHMARK(BM_Kmeans_ElkHam)->Args({n,d,5})->Args({n,d,10})->Args({n,d,15})->Args({n,d,20})->Args({n,d,25})->Args({n,d,30})->Args({n,d,35})->Args({n,d,40})->Unit(benchmark::kMillisecond)->Iterations(1);

    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({n,d,5})->Args({n,d,10})->Args({n,d,15})->Args({n,d,20})->Args({n,d,25})->Args({n,d,30})->Args({n,d,35})->Args({n,d,40})->Unit(benchmark::kMillisecond)->Iterations(1);
    

    //Scale DCT
    /*d=8*8;
    BENCHMARK(BM_load_data_Se3d)->Args({n,d,k})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);
    //BENCHMARK(BM_Kmeans_lloyd)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_Elkan)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_StepWise)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_ElkHam)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_MARIGOLD)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    
    d=16*16;
    BENCHMARK(BM_load_data_Se3d)->Args({n,d,k})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);
    //BENCHMARK(BM_Kmeans_lloyd)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_Elkan)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_StepWise)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_ElkHam)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_MARIGOLD)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    
    d=32*32;
    BENCHMARK(BM_load_data_Se3d)->Args({n,d,k})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);
    //BENCHMARK(BM_Kmeans_lloyd)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_Elkan)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_StepWise)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_ElkHam)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_MARIGOLD)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    
    d=64*64;
    BENCHMARK(BM_load_data_Se3d)->Args({n,d,k})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);
    //BENCHMARK(BM_Kmeans_lloyd)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_Elkan)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_StepWise)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_ElkHam)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_MARIGOLD)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    
    d=128*128;
    BENCHMARK(BM_load_data_Se3d)->Args({n,d,k})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);
    //BENCHMARK(BM_Kmeans_lloyd)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_Elkan)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_StepWise)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_ElkHam)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_MARIGOLD)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    
    d=256*256;
    BENCHMARK(BM_load_data_Se3d)->Args({n,d,k})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);
    //BENCHMARK(BM_Kmeans_lloyd)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_Elkan)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_StepWise)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_ElkHam)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_MARIGOLD)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    
    d=512*512;
    BENCHMARK(BM_load_data_Se3d)->Args({n,d,k})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);
    //BENCHMARK(BM_Kmeans_lloyd)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_Elkan)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_StepWise)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_ElkHam)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_MARIGOLD)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    
    d=1024*1024;
    BENCHMARK(BM_load_data_Se3d)->Args({n,d,k})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);
    //BENCHMARK(BM_Kmeans_lloyd)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_Elkan)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_StepWise)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_ElkHam)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_MARIGOLD)->Args({n,d,k})->Unit(benchmark::kMillisecond)->Iterations(1);
    */
    
    //BENCHMARK(BM_Kmeans_StepWise)->Args({n,d,5})->Args({n,d,10})->Args({n,d,15})->Args({n,d,20})->Args({n,d,25})->Args({n,d,30})->Args({n,d,35})->Args({n,d,40})->Unit(benchmark::kMillisecond)->Iterations(1);   
    //BENCHMARK(BM_Kmeans_ElkHam)->Args({n,d,5})->Args({n,d,10})->Args({n,d,15})->Args({n,d,20})->Args({n,d,25})->Args({n,d,30})->Args({n,d,35})->Args({n,d,40})->Unit(benchmark::kMillisecond)->Iterations(1);
    //BENCHMARK(BM_Kmeans_Elkan)->Args({n,d,5})->Args({n,d,10})->Args({n,d,15})->Args({n,d,20})->Args({n,d,25})->Args({n,d,30})->Args({n,d,35})->Args({n,d,40})->Unit(benchmark::kMillisecond)->Iterations(1);
    
    //BENCHMARK(BM_Kmeans_ElkHam)->Args({n,d,25})->Args({n,d,30})->Args({n,d,35})->Args({n,d,40})->Unit(benchmark::kMillisecond)->Iterations(1);
    //BENCHMARK(BM_Kmeans_Hamerly)->Iterations(1)->Unit(benchmark::kMillisecond);
    //BENCHMARK(BM_Kmeans_Elkan)->Args({n,d,5})->Args({n,d,10})->Args({n,d,15})->Args({n,d,20})->Args({n,d,25})->Args({n,d,30})->Args({n,d,35})->Args({n,d,40})->Unit(benchmark::kMillisecond)->Iterations(1);
    //BENCHMARK(BM_Kmeans_Elkan)->Args({n,d,40})->Unit(benchmark::kMillisecond)->Iterations(1);
    
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    //Dataset data(n, d);
    //data.load_datafile(data_file_name);

    //data.print_datasample();

    return 0;
};


