#include "../interfaces/interface_kmeans.hpp"
#include "../kmeans_utils/utils.cpp"
#include <cstring>
#include <algorithm>

//#include <limits>

//TODO: dissasemble in further strategies to encapsule params where they belong, using smart pointers to avoid ref loss.

class StepWiseKmeansStrategy : public KmeansStrategy {
    public:
        int* run(Dataset* data) {
            //data->print_datasample();
            //Write lloyd
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
            double* LB = new double[k];
            std::fill(LB, LB+k, 0.0);
            double UB_min = std::numeric_limits<double>::max();
            double UB;
            int *mask = new int[k];
            std::fill(mask, mask+k, 1);

            double *dist = new double[k];
            for (int j = 0; j < k; j++) {
                dist[j] = data_ss[x][0]+centroid_ss[j][0];
            }

            int mask_sum = k;

            while (l <= L && mask_sum > 1) {
                for (int j = 0; j < k; j++) {
                    if (mask[j] != 1) continue;
                    //if (a == j) continue; 

                    if (UB_min < LB[j]) {
                        mask[j] = 0;
                    } else {
                        
                        //DistToLevel(int x, int c, int d, double data[], double centroids[], double* data_ss[], double* centroid_ss[], double* dots[], int l, int L, double &UB, double &LB)
                        DistToLevel_bot(x, j, d, data_ptr, centroids, data_ss, centroid_ss, l, L, dist[j], UB, LB[j], feature_cnt, l_pow);
                        //auto val_ = Euclidian_distance(x,j,d,k,data,centroids);
                        //UB = val_;
                        //LB[j] = val_;
                        
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
            delete[] mask;
            delete[] dist;

            return a;
        }

        void clear() {}
        /*void clear() {
            

            delete div;

            
            for (int i = 0; i < k; i++) {
                delete c_to_c[i];
            }
            delete c_to_c; 

            delete l_pow;

            
            
            //delete data_ss;

            
            
            //delete centroid_ss;

            
            for (int i = 0; i < n; i++) {
                delete dots[i];
            }
            delete dots;

            delete labels;

           
            delete cluster_count;
            
            //Init distances
            delete distances; 
           
            //Init centroids  
            delete centroids;
            delete old_centroids;
            

        }*/


        void init(int _max_iter, int _n, int _d, int _k, Dataset* _data) {
            
            max_inter = _max_iter;
            n = _n;
            d = _d;
            k = _k;
            data_ptr = _data->get_data_pointer();
            feature_cnt = 0;

            //stepwise levels
            L = log10(d)/log10(4);

            div = new double[k];

            //c_to_c
            c_to_c = new double*[k];//[new double[k]];
            for (int i = 0; i < k; i++) {
                c_to_c[i] = new double[k];
            }

            l_pow = new int[L+1];
            for (int i = 0; i <= L; i++) {
                l_pow[i] = int(pow(2,i));
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

            //dots
            dots = new double*[n];
            for (int i = 0; i < n; i++) {
                dots[i] = new double[k];
                std::fill(dots[i], dots[i]+k, 0);
            }

            //Init labels
            labels = new int[n];
            std::fill(labels, labels+n, 0); 

            //Init cluster_counts
            cluster_count = new double[k];
            

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
            memcpy(centroids, data_ptr , sizeof(double)*k*d); //Initial dentroids
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
        int L;
        double* centroids;
        double* old_centroids;
        double* cluster_count;

        long long feature_cnt;

        double** dots;

        int* l_pow;

        double* div;

        double** c_to_c;
        double** data_ss;
        double** centroid_ss;

        //x to c [x*k+c]
        double* distances;
        int* labels;

        double* data_ptr;// = data->get_data_pointer();

};