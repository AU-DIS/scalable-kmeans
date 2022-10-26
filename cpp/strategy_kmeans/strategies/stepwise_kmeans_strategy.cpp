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

            double* LB = new double[k];
            double UB_min;// = std::numeric_limits<double>::max();
            double UB;
            int *mask = new int[k];
            int l;
            int a;
            int mask_sum;

            int pow_ll = pow(2,l-1);
            int pow_l = pow(2,l);
            int d_sqrt = sqrt(d);

            double dist;
            double margin;
            // calculate square data 
            //calculate_data_squares(data_ptr, data_ss, n, d)
            Calculate_squared(d, n, data_ptr, data_ss);

            /*for (int i = 0; i < k; i++) {
                for (int j = i; j < k; j++) {
                    double tmp = 0; //centroids_ss[i][0] + centroids_ss[j][0];
                    for (int f = 0; f < d; f++) {
                        //TODO: this does not use squares when it could
                        tmp += ((centroids[i*d+f] - centroids[j*d+f]) *
                            (centroids[i*d+f] - centroids[j*d+f]));
                    }
                    if(tmp < 0.0) tmp = 0.0;
                    tmp = sqrt(tmp);
                    //We can save distances for later use
                    c_to_c[i][j] = sqrt(tmp);
                    // THEY'RE THE SAME
                    c_to_c[j][i] = c_to_c[i][j];
                }
            }*/

            
            
            while ((iter < max_inter) && (!converged)) {
                //calculate square centroids
                Calculate_squared(d, k, centroids, centroid_ss);

                /*if (iter < 2) {
                     for (int i = 0; i < k; i++) {
                    for (int j = 0; j < d; j++) {
                        std::cout << centroids[i*d+j] << " ";
                    }       
                    std::cout << "\n";
                }
                std::cout << std::endl;

                }*/

                //assign to centroids
                for (int i = 0; i < n; i++) {
                    l = 0;
                    a = -1;
                    UB_min = std::numeric_limits<double>::max();
                    std::fill(LB, LB+k, 0.0);
                    std::fill(mask, mask+k, 1);

                     
                    mask_sum = k;

                    while (l <= L && mask_sum > 1) {
                        for (int j = 0; j < k; j++) {
                            if (mask[j] != 1) continue;
                            //if (a == j) continue; 

                            if (UB_min < LB[j]) {
                                mask[j] = 0;
                            } else {
                                
                                //DistToLevel(int x, int c, int d, double data[], double centroids[], double* data_ss[], double* centroid_ss[], double* dots[], int l, int L, double &UB, double &LB)
                                    //Calculate dots  
                                
                                //L is an int of log4(d), hence rounded down. 
                                //Hence when log4(d) is not in natural, d_sqrt will be bigger and the correct end for the final level instead of 2^l.

                                /*dots[x][c] = 0;
                                for (int l_ = 0; l_ < std::min(d_sqrt,(int) pow(2,l)); l_++) {
                                    for (int l_2 = 0; l_2 < std::min(d_sqrt,(int) pow(2,l)); l_2++) {
                                        dots[x][c] += data[x*d+l_*d_sqrt+l_2]*centroids[c*d+l_*d_sqrt+l_2];
                                    }
                                }*/
                                pow_ll = pow(2,l-1);
                                pow_l = pow(2,l);
                                
                                if (l==0) {
                                    dots[i][j] = data_ptr[i*d+0]*centroids[j*d+0];
                                } else {
                                    //dots saved from previous level, hence only add dots from this level.
                                    //adding new cols from from known rows
                                    for (int l_ = 0; l_ < pow_ll; l_++) {
                                        for (int l_2 = pow_ll; l_2 < pow_l ; l_2++) {
                                            dots[i][j] += data_ptr[i*d+l_*d_sqrt+l_2]*centroids[j*d+l_*d_sqrt+l_2]; 
                                        }
                                    }
                                    //TODO: sqrt stuff for d != 2^x
                                    //add full new rows
                                    for (int l_ = pow_ll; l_ < pow_l; l_++) {
                                        for (int l_2 = 0; l_2 < pow_l; l_2++) {
                                            dots[i][j] += data_ptr[i*d+l_*d_sqrt+l_2]*centroids[j*d+l_*d_sqrt+l_2]; 
                                        }
                                    }
                                }
                                
                                
                                
                                dist = data_ss[i][L] + centroid_ss[j][L] - 2*dots[i][j]; 

                                margin = 2 * sqrt((data_ss[i][L]-data_ss[i][l])*(centroid_ss[j][L]-centroid_ss[j][l]));


                                
                                LB[j] = sqrt(std::max(0.0,dist - margin));
                                UB = sqrt(std::max(0.0,dist + margin));
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
                    
                    labels[i] = a;//SetLabel(i, d, k, data_ptr, centroids, data_ss, centroid_ss, dots, L);                   
                }
                converged = Recalculate(data_ptr, centroids, old_centroids, cluster_count, labels, div, n, k, d);
                iter++;
            }   

            for (int j = 0; j < k; j++) {
                std::cout << cluster_count[j] << " ";
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

            //stepwise levels
            L = log10(d)/log10(4);

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

        double** dots;

        double* div;

        double** c_to_c;
        double** data_ss;
        double** centroid_ss;

        //x to c [x*k+c]
        double* distances;
        int* labels;

        double* data_ptr;// = data->get_data_pointer();

};