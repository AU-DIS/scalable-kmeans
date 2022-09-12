#include "interfaces/interface_kmeans.hpp"
#include <memory>

class KmeansRunner {

    private:
        std::unique_ptr<KmeansStrategy> kmeans_strategy_;

    public:
        explicit KmeansRunner(std::unique_ptr<KmeansStrategy> &&kmeans_strategy = {}) : kmeans_strategy_(std::move(kmeans_strategy)) {

        };

        //Allow replacement of strategy
        void set_strategy(std::unique_ptr<KmeansStrategy> &&kmeans_strategy) {
            kmeans_strategy_ = std::move(kmeans_strategy);
        }


        void run_kmeans() {
            kmeans_strategy_->run();
        };

};