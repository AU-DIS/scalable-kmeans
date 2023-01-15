#include "../interfaces/interface_kmeans.hpp"
#include "../kmeans_utils/utils.cpp"
#include <cstring>
#include <algorithm>

class StepWiseKmeansStrategy : public KmeansStrategy {
    public:
        int* run(Dataset* data) {

            int iter = 0;
            bool converged = false;


            // calculate square data 
            //calculate_data_squares(data_ptr, data_ss, n, d)
            Calculate_squared_botup(d, n, data_ptr, data_ss, l_pow);    
            
            while ((iter < max_inter) && (!converged)) {
                //calculate square centroids
                Calculate_squared_botup(d, k, centroids, centroid_ss, l_pow);


                //assign to centroids
                for (int i = 0; i < n; i++) {                  
                    labels[i] = SetLabel(i);//, d, k, data_ptr, centroids, data_ss, centroid_ss, dots, L, feature_cnt);                   
                }
                converged = Recalculate(data_ptr, centroids, old_centroids, cluster_count, labels, div, n, k, d, feature_cnt);
                iter++;
            }   

            for (int j = 0; j < k; j++) {
                std::cout << cluster_count[j] << " ";
            }
            std::cout << std::endl;
            std::cout << "Iter:" << iter << " Feature_cnt: " << feature_cnt << std::endl;          

            return labels;
        };

        int SetLabel(const int x) {
            int l = 0;
            int a = -1;
            
            std::fill(LB, LB+k, 0.0);
            UB_min = std::numeric_limits<double>::max();
            
            std::fill(mask, mask+k, 1);

            
            for (int j = 0; j < k; j++) {
                dist[j] = data_ss[x][0]+centroid_ss[j][0];
            }

            int mask_sum = k;

            while (l <= L && mask_sum > 1) {
                for (int j = 0; j < k; j++) {
                    if (mask[j] != 1) continue;

                    if (UB_min < LB[j]) {
                        mask[j] = 0;
                    } else {
                        
                        //DistToLevel(int x, int c, int d, double data[], double centroids[], double* data_ss[], double* centroid_ss[], double* dots[], int l, int L, double &UB, double &LB)
                        DistToLevel_bot(x, j, d, data_ptr, centroids, data_ss, centroid_ss, l, L, dist[j], UB, LB[j], feature_cnt, l_pow);
                        
                        if (UB < UB_min) {
                            a = j;
                            UB_min = UB;
                        }
                    }
                }
                mask_sum = 0;
                for (int j = 0; j < k; j++) {
                    mask_sum += mask[j];
                }
                l++;
            }
        
            return a;
        }

        //void clear() {}
        void clear() {
            
            delete[] LB;
            delete[] div;
            delete[] dist;
            delete[] mask;

            
            for (int i = 0; i < k; i++) {
                delete[] c_to_c[i];
            }
            delete[] c_to_c; 

            delete[] l_pow;

            
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
            L = ceil(log10(d)/log10(4));
            LB = new double[k];

            div = new double[k];
            dist = new double[k];
            mask = new int[k];

            //c_to_c
            c_to_c = new double*[k];//[new double[k]];
            for (int i = 0; i < k; i++) {
                c_to_c[i] = new double[k];
            }

            l_pow = new int[L+1];
            for (int i = 0; i <= L; i++) {
                if (i == L && log10(d)/log10(4) < L) {
                    l_pow[i] = sqrt(d);
                } else {
                    l_pow[i] = int(pow(2,i));
                } 
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

        double* LB;
        double UB;
        double UB_min;
        double *dist;
        int *mask;

        long long feature_cnt;

        //double** dots;

        int* l_pow;

        double* div;

        double** c_to_c;
        double** data_ss;
        double** centroid_ss;

        //x to c [x*k+c]
        //double* distances;
        int* labels;

        double* data_ptr;// = data->get_data_pointer();

};