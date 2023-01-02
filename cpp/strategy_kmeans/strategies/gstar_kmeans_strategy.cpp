#include "../interfaces/interface_kmeans.hpp"
#include "../kmeans_utils/utils.cpp"
#include <cstring>
#include <algorithm>
#include <set>
#include <tuple>
#include <vector>
#include <algorithm>
#include <utility>
#include <math.h>


class GstarKmeansStrategy : public KmeansStrategy {
    public:
        int* run(Dataset* data) {
            int iter = 0;
            bool converged = false;
            
            //Init old centers done in recalculate

            //Label each point using Alg 6
            first_center_assign();

            //Update centers with recalculate. Ignore convergence output 
            Recalculate(data_ptr, centroids, old_centroids, cluster_count, labels, div, n, k, d, feature_cnt);

            while ((iter < max_inter) && (!converged)) {              
                //Label each point using alg 3 
                center_assign();

                //Update centers with recalculate
                converged = Recalculate(data_ptr, centroids, old_centroids, cluster_count, labels, div, n, k, d, feature_cnt);
                iter++;
            }

            for (int i = 0; i < k; i++) {      
                std::cout << cluster_count[i] << " "; 
            }
            std::cout << std::endl;
            std::cout << iter << std::endl;

            return labels;
        };

        void center_assign() {
            //Alg 3
            //For all centroids C, sort other O by distance to C
            for (int c = 0; c < k; c++) {
                std::vector<std::pair<double, int>> list_of_c;
                for (int o = 0; o < k; o++) {
                    if (o == c) continue;
                    list_of_c.push_back(std::make_pair(c_to_c[c][o], o));
                }
                //Sort pairs in vector form, dump keys() to memory array. 
                //Distances can still be found in c_to_c.
                sort(list_of_c.begin(), list_of_c.end());
                for (int j = 0; j < k-1; j++) {
                    sorted_to_c[c][j] = list_of_c[j].second;
                }    
            }

            //Main assign loop
            //for each cluster center C
            for (int c = 0; c < k; c++) {
                int D = sorted_to_c[c][k/2];


                //2d coordinations
                for (int i = 0; i < k; i++) {
                    //TODO: should i == c be excluded?
                    std::tuple<double, double> cords = get_2D_coordinations(c_to_c[c][D], c_to_c[i][c], c_to_c[i][D]);
                    centers_2d_cd[c][i*2] = std::get<0>(cords);
                    centers_2d_cd[c][i*2+1] = std::get<1>(cords);
                    
                    //calculate dist to div
                    double tmp = 0;
                    for (int f = 0; f < d; f++) {
                        tmp += ((old_centroids[c*d+f] - centroids[i*d+f]) *
                            (old_centroids[c*d+f] - centroids[i*d+f]));
                    }
                    //feature_cnt += d;
                    if(tmp < 0.0) tmp = 0.0;
                    tmp = sqrt(tmp);

                    std::tuple<double, double> cords2 = get_2D_coordinations(div[c], tmp ,c_to_c[i][c]);
                    centers_2d_deltac[c][i*2] = std::get<0>(cords2);
                    centers_2d_deltac[c][i*2+1] = std::get<1>(cords2);
                }
                  
            }

            //NOTE: I moved this loop out of the main loop as we have to loop over all points anyways
            //We can hence skip grouping points by label at the cost of memory.

            //for each point p assigned C
            for (int p = 0; p < n; p++) {
                
                int D = sorted_to_c[labels[p]][k/2];
                
                //calculate p to D
                double pD = 0;
                for (int f = 0; f < d; f++) {
                    pD += ((data_ptr[p*d+f] - centroids[D*d+f]) *
                        (data_ptr[p*d+f] - centroids[D*d+f]));
                }
                //feature_cnt += d;
                if(pD < 0.0) pD = 0.0;
                pD = sqrt(pD);

                std::tuple<double, double> cords_p_2D_cd = get_2D_coordinations(c_to_c[labels[p]][D], _min[p], pD);
                
                //calculate pOldc
                double pOldc = 0;
                for (int f = 0; f < d; f++) {
                    pOldc += ((data_ptr[p*d+f] - old_centroids[labels[p]*d+f]) *
                        (data_ptr[p*d+f] - old_centroids[labels[p]*d+f]));
                }
                //feature_cnt += d;
                if(pOldc < 0.0) pOldc  = 0.0;
                pOldc  = sqrt(pOldc );

                std::tuple<double, double> cords_p_2D_deltac = get_2D_coordinations(div[labels[p]], pOldc, pD);

                //TODO Replace Near Neighbors search for candidates with KD tree
                //NNS in linear 2D scan
                //TODO: Move Set min into loop and see if we save time
                double chosen = _min[p];
                for (int j = 0; j < k; j++) {
                    double j_x_cd = centers_2d_cd[labels[p]][j*2]; 
                    double j_y_cd = centers_2d_cd[labels[p]][j*2+1];
                    double j_x_deltac = centers_2d_deltac[labels[p]][j*2]; 
                    double j_y_deltac = centers_2d_deltac[labels[p]][j*2+1];

                    if (sqrt((std::get<0>(cords_p_2D_cd)-j_x_cd)*(std::get<0>(cords_p_2D_cd)-j_x_cd)+(std::get<1>(cords_p_2D_cd)-j_y_cd)*(std::get<1>(cords_p_2D_cd)-j_y_cd)) < _min[p]) {
                        if (sqrt((std::get<0>(cords_p_2D_deltac)-j_x_deltac)*(std::get<0>(cords_p_2D_deltac)-j_x_deltac)+(std::get<1>(cords_p_2D_deltac)-j_y_deltac)*(std::get<1>(cords_p_2D_deltac)-j_y_deltac)) < _min[p]) {
                            //full distance comparison
                            //calculate dist to div
                    
                            double tmp = 0;
                            for (int f = 0; f < d; f++) {
                                tmp += ((data_ptr[p*d+f] - centroids[j*d+f]) *
                                    (data_ptr[p*d+f] - centroids[j*d+f]));
                            }
                            //feature_cnt += d;
                            if(tmp < 0.0) tmp = 0.0;
                                tmp = sqrt(tmp);

                            if (chosen > tmp) {
                                chosen = tmp;
                                labels[p] = j;
                            } 

                        } 
                    }
                }
                _min[p] = chosen;
  
            }

        }

        void first_center_assign() {
            //calculate c_to_c
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
            }

            //Choose random C_r from k
            int C_r = 0;
            //NOTE: Set to 0 to keep initialisation the same for testing

            //Find D, the k/2-th nearest center to C_r
            std::vector<std::pair<double, int>> list_of_c;
            for (int o = 0; o < k; o++) {
                if (o == C_r) continue;
                list_of_c.push_back(std::make_pair(c_to_c[C_r][o], o));
            }
            //Sort pairs in vector form, dump keys() to memory array. 
            //Distances can still be found in c_to_c.
            sort(list_of_c.begin(), list_of_c.end());

            int D = list_of_c[k/2].second;


            //For all centers i, find 2D coordinate on plane iDC_r  
            for (int j = 0; j < k; j++) {
                std::tuple<double, double> cords = get_2D_coordinations(c_to_c[C_r][D], c_to_c[j][C_r], c_to_c[j][D]); 
                centroids_2D[j*2+0] = std::get<0>(cords);
                centroids_2D[j*2+1] = std::get<1>(cords);
            }

            //Loop over all points
            for (int p = 0; p < n; p++) {
                //Compute |p, C_r|
                double tmp = 0;
                for (int f = 0; f < d; f++) {
                    tmp += ((data_ptr[p*d+f] - centroids[C_r*d+f]) *
                        (data_ptr[p*d+f] - centroids[C_r*d+f]));
                }
                //feature_cnt += d;
                if(tmp < 0.0) tmp = 0.0;
                tmp = sqrt(tmp);

                double pC_r = tmp;
                

                //Compute |p, D|
                double tmp2 = 0;
                for (int f = 0; f < d; f++) {
                    tmp2 += ((data_ptr[p*d+f] - centroids[D*d+f]) *
                        (data_ptr[p*d+f] - centroids[D*d+f]));
                }
                //feature_cnt += d;
                if(tmp2 < 0.0) tmp2 = 0.0;
                tmp2 = sqrt(tmp2);

                double pD = tmp2;

                //Get 2d cords
                std::tuple<double, double> cords = get_2D_coordinations(c_to_c[C_r][D], pC_r, pD); 
                double p2Dx = std::get<0>(cords);
                double p2Dy = std::get<1>(cords);

                //Find candidates with 2D nearest neighbor search to centroids
                //std::set<int> candidates;   
                //NOTE: In their code they loop all except C_r and D

                //Init label to closest of C_r and D
                labels[p] = C_r ? pC_r < pD : D; 

                //Init min to dist to label
                _min[p] = pC_r ? pC_r < pD : pD;

                //Init f to unassigned and declare the distance
                int F = -1;
                double pF;
 

                //For all candidates
                //NOTE: In their code they loop all excep C_r and D
                for (int O = 0; O < k; O++) {
                    if (O == C_r || O == D) continue;

                    if (F != -1) {
                        //Determine angle using alg 5
                        double theta = determine_thetha(pC_r, pD, pF, c_to_c[F][C_r], c_to_c[F][D], c_to_c[F][O], c_to_c[C_r][D], c_to_c[C_r][O], c_to_c[D][O]); 
                        
                        //Get 3D coord with alg 4
                        std::tuple<double, double, double> cords3d = get_3D_coordinations(c_to_c[C_r][D], pC_r, pD, theta);
                        double p3Dx = std::get<0>(cords3d);
                        double p3Dy = std::get<1>(cords3d);
                        double p3Dz = std::get<2>(cords3d);

                        //Calculate 3d distance as an lower bound
                        double LB = sqrt(((p3Dx-centroids_2D[O*2])*(p3Dx-centroids_2D[O*2]))+((p3Dy-centroids_2D[O*2+1])*(p3Dy-centroids_2D[O*2+1]))+((p3Dz-0)*(p3Dz-0))); //TODO

                        if (LB > _min[p]) {
                            continue;
                        }
                    }

                    //Calculate pO
                    double tmp = 0;
                    for (int f = 0; f < d; f++) {
                        tmp += ((data_ptr[p*d+f] - centroids[O*d+f]) *
                            (data_ptr[p*d+f] - centroids[O*d+f]));
                    }
                    //feature_cnt += d;
                    if(tmp < 0.0) tmp = 0.0;
                    tmp = sqrt(tmp);

                    double pO = tmp; 

                    //Set F if undefined
                    if (F == -1) {
                        F = O;
                        pF = pO;
                    } 

                    if (pO < _min[p]) {
                        _min[p] = pO;
                        labels[p] = O;
                    }
                }
            }
        }

        std::tuple<double, double> get_2D_coordinations(double M_N, double q_M, double q_N) {
            //Alg 1
            
            double cosine = ((q_N*q_N)+(q_M*q_M)-(M_N*M_N)) / 2*q_N*q_M;
            double x = cosine*q_M;
            double y = sqrt( (M_N*M_N) - (x*x));

            return {x, y};
        }

        std::tuple<double, double, double> get_3D_coordinations(double M_N, double q_M, double q_N, double theta) {
            //Alg 4
            
            double cosine = ((q_N*q_N)+(q_M*q_M)-(M_N*M_N)) / 2*q_N*q_M;
            double x = cosine*q_M;
            double y = sqrt( (M_N*M_N) - (x*x)) * sin(theta);
            double z = sqrt( (M_N*M_N) - (x*x)) - (sqrt( (M_N*M_N) - (x*x)) * sin(theta)) ;

            return {x, y, z};
        }

        double determine_thetha(double pM, double pN, double pF, double FM, double FN, double FO, double MN, double MO, double NO) {
            //Algorithm 5
            //Find FP_map
            std::tuple<double, double> p_2d = get_2D_coordinations(MN, pM, pN);
            std::tuple<double, double> F_2d = get_2D_coordinations(MN, FM, FN);
            double HpHF = abs(std::get<0>(p_2d)-std::get<0>(F_2d));
            double FP_map = sqrt((pF*pF)-(HpHF*HpHF));

            //Find angle1
            double angle1 = acos(((std::get<1>(F_2d)*std::get<1>(F_2d))+(std::get<1>(p_2d)*std::get<1>(p_2d))-(FP_map*FP_map)) / 2*std::get<1>(p_2d)*std::get<1>(F_2d));            

            //Find angle2
            std::tuple<double, double> O_2d = get_2D_coordinations(MN, MO, NO);
            double HOHF = abs(std::get<0>(O_2d)-std::get<0>(F_2d));
            double OF_map = sqrt((FO*FO)-(HOHF*HOHF));
            double angle2 = acos((std::get<1>(F_2d)*std::get<1>(F_2d))+(std::get<1>(O_2d)*std::get<1>(O_2d))-(OF_map*OF_map)/ 2*std::get<1>(F_2d)*std::get<1>(O_2d));

            //Find theta
            double theta = abs(angle1-angle2);
            
            return theta;
        }

        void clear() {
            delete[] div;        

            delete[] labels;
        
            delete[] cluster_count;

            delete[] distances; 
           
            delete[] centroids;
            delete[] old_centroids;

            delete[] centroids_2D;

            for (int i = 0; i < k; i++) {
                delete[] c_to_c[i];
            }
            delete[] c_to_c;

            for (int i = 0; i < k; i++) {
                delete[] sorted_to_c[i];
            }
            delete[] sorted_to_c;

            for (int i = 0; i < k; i++) {
                delete[] centers_2d_cd[i];
            }
            delete[] centers_2d_cd;

            for (int i = 0; i < k; i++) {
                delete[] centers_2d_deltac[i];
            }
            delete[] centers_2d_deltac;

            delete[] _min;



        }

        void init(int _max_iter, int _n, int _d, int _k, Dataset* _data) {
            
            max_inter = _max_iter;
            n = _n;
            d = _d;
            k = _k;
            data_ptr = _data->get_data_pointer();
            feature_cnt = 0;

            //Init labels
            labels = new int[n];
            std::fill(labels, labels+n, 0); 

            //Init cluster_counts
            cluster_count = new double[k];
            

            div = new double[k];

            //Init distances
            distances = new double[n*k];
            std::fill(distances, distances+n*k, std::numeric_limits<double>::max());

            //Init centroids  
            centroids = new double[k*d];
            old_centroids = new double[k*d];

            memcpy(centroids, data_ptr , sizeof(double)*k*d);

            centroids_2D = new double[k*2];

            //c_to_c
            c_to_c = new double*[k];//[new double[k]];
            for (int i = 0; i < k; i++) {
                c_to_c[i] = new double[k];
            } 

            sorted_to_c = new double*[k];//[new double[k]];
            for (int i = 0; i < k; i++) {
                sorted_to_c[i] = new double[k-1];
            }

            centers_2d_cd = new double*[k];
            for (int i = 0; i < k; i++) {
                centers_2d_cd[i] = new double[k*2];
            }

            centers_2d_deltac = new double*[k];
            for (int i = 0; i < k; i++) {
                centers_2d_deltac[i] = new double[k*2];
            }

            _min = new double[n];

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
        long long feature_cnt;

        //x to c [x*k+c]
        double* distances;
        int* labels;

        double* data_ptr;// = data->get_data_pointer();
        double* centroids_2D;

        double** c_to_c;
        double** sorted_to_c;

        double** centers_2d_cd;
        double** centers_2d_deltac;

        double* _min;

};