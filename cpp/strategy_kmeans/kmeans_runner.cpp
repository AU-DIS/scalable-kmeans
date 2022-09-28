#include "interfaces/interface_kmeans.hpp"
#include "dataset.cpp"
#include <string>
#include <memory>

class KmeansRunner {

    private:
        std::unique_ptr<KmeansStrategy> kmeans_strategy_;
        std::unique_ptr<Dataset> data;// = std::make_unique<Dataset>();

    public:
        explicit KmeansRunner(std::unique_ptr<KmeansStrategy> &&kmeans_strategy = {}) : kmeans_strategy_(std::move(kmeans_strategy)) {

        };

        //Allow replacement of strategy
        void set_strategy(std::unique_ptr<KmeansStrategy> &&kmeans_strategy) {
             
        };

        void init_run(int n, int d, int k, std::string data_file_name){
            data = std::move(std::make_unique<Dataset>(n,d));
            data->load_datafile(data_file_name);
        };


        void run_kmeans() {
            // std::unique_ptr<Dataset>(new Dataset(n,d));
            kmeans_strategy_->run(data.get());
        };

};