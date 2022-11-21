#include "../interfaces/interface_kmeans.hpp"
#include "../kmeans_utils/utils.cpp"
#include <cstring>
#include <algorithm>

class HamerlyKmeansStrategy : public KmeansStrategy {
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
                    u_elkan[i] = distance;

                    if (u_elkan[i] > val) {              
                        double smallest = std::numeric_limits<double>::max();
                        for (int j = 0; j < k; j++) {
                            double dist_ = Euclidian_distance(i, j, d, k, data_ptr, centroids, feature_cnt);
                            if (dist_ <= smallest) {
                                l_hamerly[i] = smallest;
                                labels[i] = j;
                                smallest = dist_;
                                u_elkan[i] = dist_;                                    
                            }
                        }
                    }
                }
                //MOVE CENTERS
                converged = Recalculate(data_ptr, centroids, old_centroids, cluster_count, labels, div, n, k, d, feature_cnt);
                if (!converged) {
                    Update_bounds_hamerly();
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

        void Update_bounds_hamerly() {//double data[], double centroids[], double* c_to_c[], double* centroids_ss[], double* l_elkan[], double u_elkan[], double l_hamerly[], int labels[], double div[], double near[], int n, int k, int d, long long &feature_cnt) {
            //For all x in X
            for (int i = 0; i < n; i++) { 
                int sec_smallest_id = labels[i] == 0 ? 1 : 0; 
                //For all c in C
                double sec_smallest = std::numeric_limits<double>::max();
                for (int j = 0; j < k; j++) {
                    //l_elkan(x, c) <-- max{0, l_elkan(x, c) - div[c]}  
                    if ((labels[i] != j) && (div[j] = sec_smallest)) {
                        sec_smallest = div[j];
                        sec_smallest_id = j; 
                    }  
                }
                //u_elkan <- u_elkan + div[a[x]]
                u_elkan[i] += div[labels[i]];
                //l_hamerly <- min{l_elkan(x, c != a[x])}
                l_hamerly[i] -= div[sec_smallest_id];
            }
            //Calculate centroid to centroid distances
            for (int i = 0; i < k; i++) {
                c_to_c[i][i] = 0;
                
                for (int j = i+1; j < k; j++) {
                    double tmp = 0; 
                    for (int f = 0; f < d; f++) {
                        tmp += ((centroids[i*d+f] - centroids[j*d+f]) *
                            (centroids[i*d+f] - centroids[j*d+f]));
                    }
                    //feature_cnt += d;
                    if(tmp < 0.0) tmp = 0.0;
                    tmp = sqrt(tmp);

                    
                    //We can save distances for later use
                    c_to_c[i][j] = (tmp);
                    // THEY'RE THE SAME
                    c_to_c[j][i] = c_to_c[i][j];
                }
                double smallest = DBL_MAX;
                double near_id = -1 ;
                for (int j = 0; j < k; j++) {
                if (i == j) continue;
                if (c_to_c[i][j] < smallest) {  
                        smallest = c_to_c[i][j];  
                        near[i] = 0.5 * smallest;
                        near_id = j;
                        //near[j] = 0.5 * smallest;
                    } 
                }
                //std::cout << i << " is near " <<  near_id << "\n";       
            }

            //END: Updated l_elkan, u_elkan, l_hamerly, near, c_to_c
        };

        void clear() {
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

            l_hamerly = new double[n];
            std::fill(l_hamerly, l_hamerly+n, 0.0);


            u_elkan = new double[n];
            std::fill(u_elkan, u_elkan+n, std::numeric_limits<double>::max());
            
            near = new double[k];
            std::fill(near, near+k, 0);

            div = new double[k];

            //c_to_c
            c_to_c = new double*[k];//[new double[k]];
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
            //double* data_ptr = data_ptr->get_data_pointer();

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

        double** dots;

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