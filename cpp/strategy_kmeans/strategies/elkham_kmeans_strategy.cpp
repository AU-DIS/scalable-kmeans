#include "../interfaces/interface_kmeans.hpp"
#include "../kmeans_utils/utils.cpp"
#include <cstring>
#include <algorithm>

class ElkHamKmeansStrategy : public KmeansStrategy {
    public:
        int* run(Dataset* data) {
            int iter = 0;
            bool converged = false;

            for (int i = 0; i < k; i++) {
                for (int j = i; j < k; j++) {
                    c_to_c[i][j] = 0;
                    // THEY'RE THE SAME
                    c_to_c[j][i] = 0; 
                }
            }
            
            while ((iter < max_inter) && (!converged)) {
                //assign to centroids
                for (int i = 0; i < n; i++) {
                    double val = near[labels[i]] < l_hamerly[i] ? l_hamerly[i] : near[labels[i]]; 
                    if (!(u_elkan[i] > val)) continue;
                    double distance = Euclidian_distance(i, labels[i], d, k, data_ptr, centroids, feature_cnt);  
                    l_elkan[i][labels[i]] = distance;
                    u_elkan[i] = distance;
                    for (int j = 0; j < k; j++) {
                        if (j == labels[i]) continue;
                        double val2 = std::max(l_elkan[i][j], 0.5 * c_to_c[labels[i]][j]);
                        if (u_elkan[i] > val2) {
                            double distance2 = Euclidian_distance(i, j, d, k, data_ptr, centroids, feature_cnt);  
                            l_elkan[i][j] = distance2;
                            if (distance2 < u_elkan[i]) {
                                labels[i] = j;
                                u_elkan[i] = distance2;
                            }
                        }
                    }
                }
                converged = Recalculate(data_ptr, centroids, old_centroids, cluster_count, labels, div, n, k, d, feature_cnt);
                if (!converged) {
                    Update_bounds(data_ptr, centroids, c_to_c, centroid_ss, l_elkan, u_elkan, l_hamerly, labels, div, near, n, k, d, feature_cnt);
                    iter++;
                }
            }   

            for (int j = 0; j < k; j++) {
                std::cout << cluster_count[j] << " ";
            }
            std::cout << std::endl;
            std::cout << "Iter:" << iter << " Feature_cnt: " << feature_cnt << std::endl;              

            return labels;
        };

        void clear() {
            for (int i = 0; i < n; i++) {
                delete[] l_elkan[i];
            }
            delete[] l_elkan;

            delete[] l_hamerly;
            
            delete[] u_elkan;

            delete[] near;

            delete[] div;
            
            for (int i = 0; i < k; i++) {
                delete[] c_to_c[i];
            }
            delete[] c_to_c; 
            
            for (int i = 0; i < n; i++) {
                delete[] data_ss[i];
            }
            delete[] data_ss;
            
            for (int i = 0; i < k; i++) {
                delete[] centroid_ss[i];
            }
            delete[] centroid_ss;

            delete[] labels;
      
            delete[] cluster_count;
                
            //Init centroids  
            delete[] centroids;
            delete[] old_centroids;
        }


        void init(int _max_iter, int _n, int _d, int _k, Dataset* _data) {
            
            max_inter = _max_iter;
            n = _n;
            d = _d;
            k = _k;
            data_ptr = _data->get_data_pointer();

            feature_cnt = 0;

            //stepwise levels
            L = log10(d)/log10(4);

            //bounds
            l_elkan = new double*[n];
            for (int i = 0; i < n; i++) {
                l_elkan[i] = new double[k];
                std::fill(l_elkan[i], l_elkan[i]+k, 0.0);
            }

            l_hamerly = new double[n];
            std::fill(l_hamerly, l_hamerly+n, 0.0);


            u_elkan = new double[n];
            std::fill(u_elkan, u_elkan+n, std::numeric_limits<double>::max());
            
            near = new double[k];
            std::fill(near, near+k, 0);

            div = new double[k];

            //c_to_c
            c_to_c = new double*[k];
            for (int i = 0; i < k; i++) {
                c_to_c[i] = new double[k];
            }


            //squared
            data_ss = new double*[n];
            for (int i = 0; i < n; i++) {
                data_ss[i] = new double[L+2];
            }

            centroid_ss = new double*[k];
            for (int i = 0; i < k; i++) {
                centroid_ss[i] = new double[L+2];
            }

            //Init labels
            labels = new int[n];
            std::fill(labels, labels+n, 0); 

            //Init cluster_counts
            cluster_count = new double[k];

            //Init centroids  
            centroids = new double[k*d];
            old_centroids = new double[k*d];

            memcpy(centroids, data_ptr , sizeof(double)*k*d); //Initial dentroids
        }
    private:

        int max_inter;
        int n;
        int d;
        int k;
        int L;
        double* centroids;
        double* old_centroids;
        double* cluster_count;

        long long feature_cnt;

        double** l_elkan;
        double* l_hamerly;
        double* u_elkan;

        double* near;

        double* div;

        double** c_to_c;
        double** data_ss;
        double** centroid_ss;

        int* labels;

        double* data_ptr;// = data->get_data_pointer();

};