#include "interfaces/interface_kmeans.hpp"
#include "dataset.cpp"
#include <string>
#include <memory>

class KmeansRunner {

    private:
        std::unique_ptr<KmeansStrategy> kmeans_strategy_;
        std::unique_ptr<Dataset> data;// = std::make_unique<Dataset>();
        int* result_labels;
        

    public:
        explicit KmeansRunner(std::unique_ptr<KmeansStrategy> &&kmeans_strategy = {}) : kmeans_strategy_(std::move(kmeans_strategy)) {

        };

        //Allow replacement of strategy
        void set_strategy(std::unique_ptr<KmeansStrategy> &&kmeans_strategy) {
            kmeans_strategy_ = std::move(kmeans_strategy);
        };

        void init_run(int n, int d, int k, std::string data_file_name){
            data = std::move(std::make_unique<Dataset>(n,d));
            data->load_datafile(data_file_name);
            
        };

        void init_run(int n, int d, int k, double* data_pnt){
            data = std::move(std::make_unique<Dataset>(n,d,data_pnt));
            
            //data->load_datafile(data_pnt);
            std::cout << "Print sample from c++: " << std::endl;
            data->print_datasample(n, d);
        };


        void run_kmeans(int n, int d, int k) {
            // std::unique_ptr<Dataset>(new Dataset(n,d));
            kmeans_strategy_->init(1000,n, d, k, data.get());
            result_labels = kmeans_strategy_->run(data.get());
            kmeans_strategy_->clear();
        };

        int* run_wo_clear_kmeans(int n, int d, int k) {
            // std::unique_ptr<Dataset>(new Dataset(n,d));
            kmeans_strategy_->init(1000,n, d, k, data.get());
            result_labels = kmeans_strategy_->run(data.get());
            return result_labels;
        };

        void clear() {
            kmeans_strategy_->clear();
        }

        void save_result(int n, std::string label_file_name) {
            std::ofstream label_file;
            label_file.open(label_file_name.c_str());

            for (int i = 0; i < n; i++) {
                label_file << result_labels[i] << "\n";
            }
            label_file.close();
        }

};