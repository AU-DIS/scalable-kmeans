#include "../interfaces/interface_kmeans.hpp"


class LloydKmeansStrategy : public KmeansStrategy {
    public:
        void run(Dataset* data) {
            data->print_datasample();
            //Write lloyd
        };
    private:
        void init() {
            //Init centroids

        }
};