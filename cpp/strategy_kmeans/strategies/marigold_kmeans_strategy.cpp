#include "../interfaces/interface_kmeans.hpp"
#include "../kmeans_utils/utils.cpp"
#include <cstring>
#include <algorithm>
//#include <limits>

//TODO: dissasemble in further strategies to encapsule params where they belong, using smart pointers to avoid ref loss.

class MARIGOLDKmeansStrategy : public KmeansStrategy {
    public:
        int* run(Dataset* data) {
            //data->print_datasample();
            //Write lloyd
            int iter = 0;
            bool converged = false;
           


            // calculate square data 
            //calculate_data_squares(data_ptr, data_ss, n, d)
            Calculate_squared_botup(d, n, data_ptr, data_ss, l_pow);

            
            for (int i = 0; i < k; i++) {
                for (int j = i; j < k; j++) {
                    /*double tmp = 0; //centroids_ss[i][0] + centroids_ss[j][0];
                    for (int f = 0; f < d; f++) {
                        //TODO: this does not use squares when it could
                        tmp += ((centroids[i*d+f] - centroids[j*d+f]) *
                            (centroids[i*d+f] - centroids[j*d+f]));
                    }
                    //feature_cnt += d;
                    if(tmp < 0.0) tmp = 0.0;
                    tmp = sqrt(tmp);*/
                    //We can save distances for later use
                    c_to_c[i][j] = 0; //(tmp);
                    // THEY'RE THE SAME
                    c_to_c[j][i] = 0; //c_to_c[i][j];
                }
            }
            //Recalculate(data_ptr, centroids, old_centroids, cluster_count, labels, div, n, k, d);
            //Update_bounds(data_ptr, centroids, c_to_c, centroid_ss, l_elkan, u_elkan, l_hamerly, labels, div, near, n, k, d);
            
            while ((iter < max_inter) && (!converged)) {
                //calculate square centroids
                Calculate_squared_botup(d, k, centroids, centroid_ss, l_pow);    

                //assign to centroids
                for (int i = 0; i < n; i++) {
                    double val = near[labels[i]] < l_hamerly[i] ? l_hamerly[i] : near[labels[i]]; 
                    if (u_elkan[i] > val) {
                         MG_SetLabel(i); 
                    }
                }
                converged = Recalculate(data_ptr, centroids, old_centroids, cluster_count, labels, div, n, k, d, feature_cnt);
                if (!converged) {
                    //TODO: refactor location of .. you know the drill 
                    Update_bounds(data_ptr, centroids, c_to_c, centroid_ss, l_elkan, u_elkan, l_hamerly, labels, div, near, n, k, d, feature_cnt);                   
                    
                }
                iter++;
            }   

            for (int j = 0; j < k; j++) {
                std::cout << cluster_count[j] << " ";
            }
            std::cout << std::endl;
            std::cout << "Iter:" << iter << " Feature_cnt: " << feature_cnt << std::endl;
                

            return labels;
        };

        void MG_SetLabel(const int x) {
            int l = 0;
            int *mask = new int[k];
            std::fill_n(mask, k, 1);

            double *dist = new double[k];
            for (int j = 0; j < k; j++) {
                dist[j] = data_ss[x][0]+centroid_ss[j][0];
            }
            
            
            double val;
            double UB, LB;
            
            int mask_sum = k;

            while (l <= L && mask_sum > 1) {
                for (int j = 0; j < k; j++) {
                    if (mask[j] != 1) continue;  

                    //Elkan prune
                    val = std::max(l_elkan[x][j], 0.5 * c_to_c[labels[x]][j]);// l_elkan[x][j] < 0.5 * c_to_c[labels[x]][j] ? 0.5 * c_to_c[labels[x]][j] : l_elkan[x][j];  
                    if (u_elkan[x] < val) {     //Elkan check
                        mask[j] = 0;            //Mark as pruned centroid
                    } else {
                        //DistToLevel params (int x, int c, int d, double data[], double centroids[], double* data_ss[], double* centroid_ss[], double* dots[], int l, int L)
                        DistToLevel_bot(x, j, d, data_ptr, centroids, data_ss, centroid_ss, l, L, dist[j], UB, LB, feature_cnt, l_pow);
                        LB = sqrt(std::max(0.0, LB));                
                        if (LB > l_elkan[x][j]) {
                            
                            if (LB > l_elkan[x][j]) {
                                l_elkan[x][j] = LB; //Keep maximum LB per c
                            }   
                        }
                        
                        UB = sqrt(std::max(0.0, UB));
                        if (UB < u_elkan[x]) {
                            labels[x] = j;
                            u_elkan[x] = UB; //Keep minimum UB across c
                        }       
                    } 
                }
                mask_sum = 0;
                for (int j = 0; j < k; j++) {
                    mask_sum += mask[j];
                }
                l++;
            }

            //delete[] mask;
            //delete[] dist;
            //END: Updated labels, l_elkan[x][.], u_elkan[x]
        }

        void clear() {}
        /*void clear() {
            for (int i = 0; i < n; i++) {
                delete l_elkan[i];
            }
            delete l_elkan;

            delete l_hamerly;
            
            delete u_elkan;

            delete near;

            delete div;

            
            for (int i = 0; i < k; i++) {
                delete c_to_c[i];
            }
            delete c_to_c; 

            

            
            for (int i = 0; i < n; i++) {
                delete data_ss[i];
            }
            delete data_ss;

            for (int i = 0; i < k; i++) {
                delete centroid_ss[i];
            }
            delete centroid_ss;

            
            //for (int i = 0; i < n; i++) {
            //    delete[] dots[i];
            //}
            //delete[] dots;

            delete labels;

           
            delete cluster_count;
            
           
            //Init centroids  
            delete centroids;
            delete old_centroids;
            
            delete l_pow;
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

            //bounds
            l_elkan = new double*[n];
            for (int i = 0; i < n; i++) {
                l_elkan[i] = new double[k];
                std::fill(l_elkan[i], l_elkan[i]+k, 0.0);
            }

            l_hamerly = new double[n];
            std::fill(l_hamerly, l_hamerly+n, 0);


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
            //dots = new double*[n];
            //for (int i = 0; i < n; i++) {
            //    dots[i] = new double[k];
            //    std::fill(dots[i], dots[i]+k, 0);
            //}

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

        //double** dots;

        double** l_elkan;
        double* l_hamerly;
        double* u_elkan;


        int* l_pow;

        double* near;

        double* div;

        double** c_to_c;
        double** data_ss;
        double** centroid_ss;

        //x to c [x*k+c]
        //double* distances;
        int* labels;

        double* data_ptr;// = data->get_data_pointer();

};