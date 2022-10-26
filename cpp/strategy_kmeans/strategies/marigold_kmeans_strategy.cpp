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
            int l;
            int *mask = new int[k];

            int pow_ll = pow(2,l-1);
            int pow_l = pow(2,l);

            double UB, LB;
            int mask_sum;
            int d_sqrt = sqrt(d);

            double dist;
            double margin;
            double val_2;


            // calculate square data 
            //calculate_data_squares(data_ptr, data_ss, n, d)
            Calculate_squared(d, n, data_ptr, data_ss);

            for (int i = 0; i < k; i++) {
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
            }
            //Recalculate(data_ptr, centroids, old_centroids, cluster_count, labels, div, n, k, d);
            //Update_bounds(data_ptr, centroids, c_to_c, centroid_ss, l_elkan, u_elkan, l_hamerly, labels, div, near, n, k, d);

            while ((iter < max_inter) && (!converged)) {
                //calculate square centroids
                Calculate_squared(d, k, centroids, centroid_ss);

                //assign to centroids
                for (int i = 0; i < n; i++) {
                    double val = near[labels[i]] < l_hamerly[i] ? l_hamerly[i] : near[labels[i]]; 
                    if (u_elkan[i] > val) {
                        //TODO: refactor placement of implementation to avoid bassillion arguments
                        //params = (int x, int d, int k, double data[],  double centroids[], double* data_ss[], double* centroid_ss[], double* dots[], int L, int labels[], double* l_elkan[], double u_elkan[], double* c_to_c[])
                            l = 0;
                            std::fill_n(mask, k, 1);
                            /*for (int i = 0; i < k; i++) {
                                mask[i] = 1;
                            }*/
                            
                            mask_sum = k;

                            while (l <= L && mask_sum > 1) {
                                for (int j = 0; j < k; j++) {
                                    if (mask[j] != 1) continue;  
                                    
                                    //Elkan prune
                                    val_2 = std::max(l_elkan[i][j], 0.5 * c_to_c[labels[i]][j]);// l_elkan[x][j] < 0.5 * c_to_c[labels[x]][j] ? 0.5 * c_to_c[labels[x]][j] : l_elkan[x][j];  
                                    if (u_elkan[i] < val_2) {     //Elkan check
                                        mask[j] = 0;            //Mark as pruned centroid
                                    } else {
                                        //DistToLevel params (int x, int c, int d, double data[], double centroids[], double* data_ss[], double* centroid_ss[], double* dots[], int l, int L)
                                        
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


                                        
                                        LB = sqrt(std::max(0.0,dist - margin));
                                        UB = sqrt(std::max(0.0,dist + margin));
                                        //auto val_ = Euclidian_distance(x,j,d,k,data,centroids);
                                        //UB = val_;
                                        //LB = val_;
                                        if (LB > l_elkan[i][j]) {
                                            l_elkan[i][j] = LB; //Keep maximum LB per c
                                        }
                                        if (UB < u_elkan[i]) {
                                            labels[i] = j;
                                            u_elkan[i] = UB; //Keep minimum UB across c
                                        } 

                                    } 
                                }
                                mask_sum = 0;
                                for (int j = 0; j < k; j++) {
                                    mask_sum += mask[j];
                                }
                                l++;
                            }    
                    }
                }
                converged = Recalculate(data_ptr, centroids, old_centroids, cluster_count, labels, div, n, k, d);
                if (!converged) {
                    //TODO: refactor location of .. you know the drill 
                    //(double data[], double centroids[], double* c_to_c[], double* centroids_ss[], double* l_elkan[], double u_elkan[], double l_hamerly[], int labels[], double div[], double near[], int n, int k, int d) 
                    Update_bounds(data_ptr, centroids, c_to_c, centroid_ss, l_elkan, u_elkan, l_hamerly, labels, div, near, n, k, d);                   
                }
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

            /*for (int i = 0; i < n; i++) {
                l_hamerly[i] = 0;
                u_elkan[i] = std::numeric_limits<double>::max();
            }*/
            
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

        double** l_elkan;
        double* l_hamerly;
        double* u_elkan;

        double* near;

        double* div;

        double** c_to_c;
        double** data_ss;
        double** centroid_ss;

        //x to c [x*k+c]
        double* distances;
        int* labels;

        double* data_ptr;// = data->get_data_pointer();

};