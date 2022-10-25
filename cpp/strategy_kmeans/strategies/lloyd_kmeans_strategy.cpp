#include "../interfaces/interface_kmeans.hpp"
#include "../kmeans_utils/utils.cpp"
#include <cstring>
#include <algorithm>
//#include <limits>


class LloydKmeansStrategy : public KmeansStrategy {
    public:
        int* run(Dataset* data) {
            //data->print_datasample();
            //Write lloyd
            int iter = 0;
            bool converged = false;
            
            while ((iter < max_inter) && (!converged)) {
                /*if (iter < 2) {
                     for (int i = 0; i < k; i++) {
                    for (int j = 0; j < d; j++) {
                        std::cout << centroids[i*d+j] << " ";
                    }       
                    std::cout << "\n";
                }
                std::cout << std::endl;

                }*/
               

                
                //std::cout << iter << std::endl;
                //Calculate all distances
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < k; j++) {
                        distances[i*k+j] = Euclidian_distance(i, j, d, k, data_ptr, centroids);
                        if (distances[i*k+j] < distances[i*k+labels[i]]) {
                            labels[i] = j;
                        }
                    }
                }
                //Update labels
                /*for (int i = 0; i < n; i++) {
                    for (int j = 0; j < k; j++) {
                        //std::cout << distances[i*k+j] << " " << distances[i*k+labels[i]] << "\n"; 
                        if (distances[i*k+j] < distances[i*k+labels[i]]) {
                            labels[i] = j;
                        }
                    }
                }*/
                
                converged = Recalculate(data_ptr, centroids, old_centroids, cluster_count, labels, div, n, k, d);

                /*
                //Save centroids for convergence test    
                memcpy(old_centroids, centroids, sizeof(double)*k*d);
                
                //Wipe memory for new centroid calculations
                memset(centroids, 0.0, sizeof(double)*k*d);
                //for(int i = 0; i < k; i++) {
                //    for (int j = 0; j < 10; j++) {
                //        std::cout << centroids[i*k+j] << " ";
                //    }
                //    std::cout << "\n";
                //}
                //std::cout << std::endl;

                //std::fill(centroids, centroids+k*d, 0);
                
                //Count size of clusters and add pos to centroid
                memset(cluster_count, 0, sizeof(double)*k);
                //std::fill(cluster_count, cluster_count+k, 0);
                for (int i = 0; i < n; i++) {
                    cluster_count[labels[i]]++;
                    for (int j = 0; j < d; j++) { 
                        centroids[labels[i]*d+j] += data_ptr[i*d+j];
                    }
                }

                //for(int i = 0; i < k; i++) {
                //    for (int j = 0; j < 10; j++) {
                //        std::cout << centroids[i*d+j] << " ";
                //    }
                //    std::cout << "\n";
                //}
                //std::cout << std::endl;

                //Calculate new centroid positions
                for (int i = 0; i < k; i++) {  
                    if (cluster_count[i] > 0) {
                        for (int j = 0; j < d; j++) {
                            centroids[i*d+j] /= cluster_count[i]; 
                        }
                    } else {
                        for (int j = 0; j < d; j++) {
                            centroids[i*d+j] = old_centroids[i*d+j];
                        }
                    }
                }

                for (int i = 0; i < k; i++) {      
                std::cout << cluster_count[i] << " "; 
                }
                std::cout << std::endl;
                

                //Check convergence
                converged = true;
                for (int i = 0; i < k; i++){
                    for (int j = 0; j < d; j++) {
                        double check_val = abs(old_centroids[i*d+j]-centroids[i*d+j]);
                        if (!(check_val == 0)) {
                            converged = false;
                            break;
                        }
                    }
                }*/
                iter++;
            }

            for (int i = 0; i < k; i++) {      
                std::cout << cluster_count[i] << " "; 
            }
            std::cout << std::endl;
            std::cout << iter << std::endl;

            return labels;
        };

        void init(int _max_iter, int _n, int _d, int _k, Dataset* _data) {
            
            max_inter = _max_iter;
            n = _n;
            d = _d;
            k = _k;
            data_ptr = _data->get_data_pointer();

            //Init labels
            labels = new int[n];
            std::fill(labels, labels+n, 0); 

            //Init cluster_counts
            cluster_count = new double[k];
            

            div = new double[k];

            //Init distances
            distances = new double[n*k];
            std::fill(distances, distances+n*k, std::numeric_limits<double>::max());
            //memset(distances, std::numeric_limits<double>::max(), sizeof(double)*n*k);

            //Init centroids  
            centroids = new double[k*d];
            old_centroids = new double[k*d];
            //double* data_ptr = data_ptr->get_data_pointer();
            /*for (int i = 0; i < k; i++) {
                for (int j = 0; j < d; j++) {
                    centroids[i*d+j] = data_ptr[i*d+j];
                }
            }*/
            memcpy(centroids, data_ptr , sizeof(double)*k*d);
            /*std::cout << "Printing centroids first 10 d\n";
            for (int i = 0; i < k; i++) {
                for (int j = 0; j < 10; j++) {
                    std::cout << centroids[i*d+j] << " ";
                }
                std::cout << "\n";
            }
            std::cout << std::endl;*/
        }
    private:

        int max_inter;
        int n;
        int d;
        int k;
        double* div;
        double* centroids;
        double* old_centroids;
        double* cluster_count;

        //x to c [x*k+c]
        double* distances;
        int* labels;

        double* data_ptr;// = data->get_data_pointer();

};