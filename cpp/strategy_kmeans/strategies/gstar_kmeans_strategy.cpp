#include "../interfaces/interface_kmeans.hpp"
#include "../kmeans_utils/utils.cpp"
#include <cstring>
#include <algorithm>
#include <set>
#include <tuple>


class GstarKmeansStrategy : public KmeansStrategy {
    public:
        int* run(Dataset* data) {
            int iter = 0;
            bool converged = false;
            
            //TODO: Init old centers (done in init function?)

            //Label each point using Alg 6
            first_center_assign();

            //Update centers with recalculate
            //TODO

            while ((iter < max_iter) && (!converged)) {              
                //Label each point using alg 3 
                center_assign()

                //Update centers with recalculate
                //TODO
                
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
            for (int j = 0; j < k; j++) {
 

            }
        }

        void first_center_assign() {
            //calculate c_to_c
            //TODO

            //Choose random C_r from k
            int C_r = 0;
            //TODO

            //Find D, the k/2-th nearest center to C_r
            int D; 
            //TODO

            //For all centers i, find 2D coordinate on plane iDC_r  
            for (int j = 0; j < k; j++) {
                std::tuple<double, double> cords = get_2D_coordinations(c_to_c[C_r][D], c_to_c[j][C_r], c_to_c[j][D]); 
                centroids_2D[j*2+0] = std::get<0>(cords);
                centroids_2D[j*2+1] = std::get<1>(cords);
            }

            //Loop over all points
            for (int p = 0; p < n; p++) {
                //Compute |p, C_r|
                double pC_r;
                //TODO

                //Compute |p, D|
                double pD;
                //TODO

                std::tuple<double, double> cords = get_2D_coordinations(c_to_c[C_r][D], pC_r, pD); 
                double p2Dx = std::get<0>(cords);
                double p2Dy = std::get<1>(cords);

                //Find candidates with 2D nearest neighbor search to centroids
                std::set<int> candidates;
                //TODO: Check their code 

                //Init label to closest of C_r and D
                labels[p] = C_r ? pC_r < pD : D; 

                //Init min to dist to label
                double _min = pC_r ? pC_r < pD : pD;

                //Init f to unassigned and declare the distance
                int F = -1;
                double pF;
 

                //For all candidates
                for (int O : candidates) {
                    

                    if (F != -1) {
                        //Determine angle using alg 5
                        double theta = determine_thetha(pC_r, pD, pF, c_to_c[F][C_r], c_to_c[F][D], c_to_c[F][O], c_to_c[C_r][D], c_to_c[C_r][O], c_to_c[D][O]); 
                        
                        //Get 3D coord with alg 4
                        std::tuple<double, double, double> cords3d = get_3D_coordination(c_to_c[C_r][D], pC_r, pD, theta);
                        double p3Dx = std::get<0>(cords3d);
                        double p3Dy = std::get<1>(cords3d);
                        double p3Dz = std::get<2>(cords3d);

                        //Calculate 3d distance as an lower bound
                        double LB; //TODO

                        if (LB > _min) {
                            continue;
                        }
                    }
                    //Calculate pO
                    double pO; 
                    //TODO

                    //Set F if undefined
                    if (F = -1) {
                        F = O;
                        pF = pO;
                    } 

                    if (pO < _min) {
                        _min = pO;
                        labels[p] = O;
                    }
                }
            }
        }

        std::tuple<double, double> get_2D_coordinations(double C_rD, double pC_r, double pD) {
            double x;
            double y;

                //TODO:

            return {x, y};
        }

        std::tuple<double, double> get_3D_coordinations(double C_rD, double pC_r, double pD, double theta) {
            double x;
            double y;
            double z;

                //TODO:

            return {x, y, z};
        }

        double determine_thetha(double pC_r, double pD, double pF, double FC_r, double FD, double FO, double C_rD, double C_rO, double DO) {
            //Algorithm 5
            double theta;

            //TODO
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

};