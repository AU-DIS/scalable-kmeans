#ifdef __WIN32__
#ifdef BUILD_LIB
#define LIB_FUNC __declspec(dllexport)
#else
#define LIB_FUNC __declspec(dllimport)
#endif
#else
#define LIB_FUNC       // Blank
#endif

#include "kmeans_runner.cpp"
#include "strategies/marigold_kmeans_strategy.cpp"

KmeansRunner runner(std::make_unique<MARIGOLDKmeansStrategy>());

LIB_FUNC int* run(const char* file_name, int n, int d, int k) {   
    runner.init_run(n,d,k,file_name);
    return runner.run_wo_clear_kmeans(n,d,k);
}

LIB_FUNC int* run(double* data_pnt, int n, int d, int k) {   
    runner.init_run(n,d,k,data_pnt);
    return runner.run_wo_clear_kmeans(n,d,k);
}

LIB_FUNC void clear() {
    runner.clear();
}

/*int main(int argc, char **argv) { 
    const char* path = "/mnt/c/Users/kaspe/OneDrive/Skrivebord/Reps/scalable-kmeans/data/steinn_14_jun/processed/misfit_VB_128_h5_dct.txt";
    run(path,168,128*128,5); 
    return 0;
}*/