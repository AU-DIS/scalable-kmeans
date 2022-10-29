#ifndef UTILS_CPP
#define UTILS_CPP
#include <cfloat>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <iostream>
//#include <memory>
#include <cstring>

//TODO: In a pefect world this file does not exist and theese functions nicely exist within relevant strategies.

double Euclidian_distance(int x, int c, int d, int k, double data[], double centroids[], long long &feature_cnt) {
    //x datapoint index
    //c centroid index

    //d features/dimension
    double dist = 0;
    for (int j = 0; j < d; j++) {
        dist += (data[x*d+j]-centroids[c*d+j])*(data[x*d+j]-centroids[c*d+j]);
    };
    //std::cout << dist << std::endl;
    feature_cnt += d;
    if(dist < 0.0) dist = 0.0;
    return sqrt(dist);
}

double Squared_euclidian_distance(int x, int c, int d, int k, double data[], double centroids[], long long &feature_cnt) {
    //x datapoint index
    //c centroid index

    //d features/dimension
    double dist = 0;
    for (int j = 0; j < d; j++) {
        dist += (data[x*d+j]-centroids[c*d+j])*(data[x*d+j]-centroids[c*d+j]);
    };
    //std::cout << dist << std::endl;
    feature_cnt += d;
    if(dist < 0.0) dist = 0.0;
    return dist;
}


                                                                        //TODO: This parsing is insane and can be cleaned clean up, by defining this and dist to level within a strategy scope.   
void Update_bounds(double data[], double centroids[], double* c_to_c[], double* centroids_ss[], double* l_elkan[], double u_elkan[], double l_hamerly[], int labels[], double div[], double near[], int n, int k, int d, long long &feature_cnt) {
    //For all x in X
    for (int i = 0; i < n; i++) { 
        int smallest_id = labels[i] == 0 ? 1 : 0; 
        //For all c in C
        for (int j = 0; j < k; j++) {
            //l_elkan(x, c) <-- max{0, l_elkan(x, c) - div[c]} 
            double val = l_elkan[i][j]-div[j];  
            l_elkan[i][j] = 0 < val ? val : 0;
            if ((labels[i] != j) && (l_elkan[i][j] <= l_elkan[i][smallest_id])) {
                smallest_id = j; 
            }  
        }
        //u_elkan <- u_elkan + div[a[x]]
        u_elkan[i] += div[labels[i]];
        //l_hamerly <- min{l_elkan(x, c != a[x])}
        l_hamerly[i] = l_elkan[i][smallest_id];
    }
    //Calculate centroid to centroid distances
    for (int i = 0; i < k; i++) {
        c_to_c[i][i] = 0;
        
        for (int j = i+1; j < k; j++) {
            double tmp = 0; //centroids_ss[i][0] + centroids_ss[j][0];
            for (int f = 0; f < d; f++) {
                //TODO: this does not use squares when it could, ot could it?
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
    //std::cout << std::endl;
    //Nearest other centroid
    //THIS CANNOT BE DONE IN C_TO_C LOOP! IF SO IT IS ORDER DEPENDENT!
    /*for (int i = 0; i < k; i++) {      
        double smallest = DBL_MAX;     
        for (int j = 0; j < k; j++) {
           if (i == j) continue;
           if (c_to_c[i][j] < smallest) {  
                smallest = c_to_c[i][j];  
                near[i] = 0.5 * smallest;
            } 
        }
    }*/

    //END: Updated l_elkan, u_elkan, l_hamerly, near, c_to_c
};


bool Recalculate(double data[], double centroids[], double old_centroids[], double cluster_count[], int labels[], double div[], int n, int k, int d, long long &feature_cnt) {
    bool converged = true;
    //With the power of clever indexing we can use memcpy, to save old centroids. (If we are even smarter we swap the pointers and save a pass)   
    memcpy(old_centroids, centroids, sizeof(double)*k*d);
    //Wipe memory for new centroid calculations
    memset(centroids, 0.0, sizeof(double)*k*d);

    memset(cluster_count, 0, sizeof(double)*k);
    
    //Count size of clusters and add pos to centroid
    for (int i = 0; i < n; i++) {
        cluster_count[labels[i]]++;
        for (int j = 0; j < d; j++) { 
            centroids[labels[i]*d+j] += data[i*d+j];
        }
    }

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

    /*for (int i = 0; i < k; i++) {      
        std::cout << cluster_count[i] << " "; 
        }
    std::cout << std::endl;*/

    //calculate div
    for (int j = 0; j < k; j++) {
        double tmp = 0.0;
        for (int f = 0; f < d; f++) {
            tmp += ((centroids[j*d+f] - old_centroids[j*d+f]) *
                    (centroids[j*d+f] - old_centroids[j*d+f]));
        }
        feature_cnt += d;
        if(tmp < 0.0) tmp = 0.0;
        div[j] = sqrt(tmp);
        if (div[j] > 0) { //This convergence check is slow and can possibly be done better using the old check
            converged = false;
        }
    }



    //END: div[.] updated
    return converged;
}


std::tuple<double, double> DistToLevel(const int x, const int c, const int d, const double data[], const double centroids[], const double *const data_ss[], const double *const centroid_ss[], const int l, const int L, double* dots[], double &UB, double &LB, long long &feature_cnt) {
    //Calculate dots  
    int d_sqrt = sqrt(d);
    //L is an int of log4(d), hence rounded down. 
    //Hence when log4(d) is not in natural, d_sqrt will be bigger and the correct end for the final level instead of 2^l.

    /*dots[x][c] = 0;
    for (int l_ = 0; l_ < std::min(d_sqrt,(int) pow(2,l)); l_++) {
        for (int l_2 = 0; l_2 < std::min(d_sqrt,(int) pow(2,l)); l_2++) {
            dots[x][c] += data[x*d+l_*d_sqrt+l_2]*centroids[c*d+l_*d_sqrt+l_2];
        }
    }*/
    int pow_ll = pow(2,l-1);
    int pow_l = pow(2,l);
    
    if (l==0) {
        dots[x][c] = data[x*d+0]*centroids[c*d+0];
    } else {
        //dots saved from previous level, hence only add dots from this level.
        //adding new cols from from known rows
        for (int l_ = 0; l_ < pow_ll; l_++) {
            for (int l_2 = pow_ll; l_2 < pow_l ; l_2++) {
                dots[x][c] += data[x*d+l_*d_sqrt+l_2]*centroids[c*d+l_*d_sqrt+l_2]; 
            }
        }
        //TODO: sqrt stuff for d != 2^x
        //add full new rows
        for (int l_ = pow_ll; l_ < pow_l; l_++) {
            for (int l_2 = 0; l_2 < pow_l; l_2++) {
                dots[x][c] += data[x*d+l_*d_sqrt+l_2]*centroids[c*d+l_*d_sqrt+l_2]; 
            }
        }   
    }
    feature_cnt += (2*pow_l)-(2*pow_ll);
    
    
    double dist = data_ss[x][L] + centroid_ss[c][L] - 2*dots[x][c]; 

    double margin = 2 * sqrt((data_ss[x][L]-data_ss[x][l])*(centroid_ss[c][L]-centroid_ss[c][l]));


    
    LB = dist - margin;//sqrt(std::max(0.0,dist - margin));
    UB = dist + margin;//sqrt(std::max(0.0,dist + margin));

    return {LB, UB};
};

//TODO: This parsing is insane can can be cleaned clean up, by defining this and dist to level within a strategy scope. 
/*void old_MG_SetLabel(int x, int d, int k, double data[],  double centroids[], double* data_ss[], double* centroid_ss[], double* dots[], int L, int labels[], double* l_elkan[], double u_elkan[], double* c_to_c[], long long &feature_cnt) {
    int l = 0;
    int *mask = new int[k];
    std::fill_n(mask, k, 1);
    //for (int i = 0; i < k; i++) {
    //    mask[i] = 1;
    //}
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
                DistToLevel(x, j, d, data, centroids, data_ss, centroid_ss, l, L, dots, UB, LB,  feature_cnt);
                                    
                if (LB > l_elkan[x][j]) {
                    LB = sqrt(std::max(0.0, LB));
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
    //TODO: free mask?
    //END: Updated labels, l_elkan[x][.], u_elkan[x]
}*/

/*void old_MG_SetLabel_test(int x, int d, int k, double data[],  double centroids[], double* data_ss[], double* centroid_ss[], double* dots[], int L, int labels[], double* l_elkan[], double u_elkan[], double l_hamerly[], double* c_to_c[], long long &feature_cnt) {
    int l = 0;
    int *mask = new int[k];
    std::fill_n(mask, k, 1);
    //for (int i = 0; i < k; i++) {
    //    mask[i] = 1;
    //}
    double val;
    double UB, LB;
    //double* LB = new double[k];
    //std::fill_n(LB, k, 0);
    //double LB_min = ; 
    double UB_min = std::numeric_limits<double>::max(); // ? std::numeric_limits<double>::max() : u_elkan[x]*u_elkan[x];
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
                DistToLevel(x, j, d, data, centroids, data_ss, centroid_ss, l, L, dots, UB, LB, feature_cnt);
                                    
                if (LB > l_elkan[x][j]) {
                    LB = sqrt(std::max(0.0, LB));
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

    //TODO: free mask?
    //END: Updated labels, l_elkan[x][.], u_elkan[x]
}*/




/*int old_SetLabel(int x, int d, int k, double data[], double centroids[], double* data_ss[], double* centroid_ss[], double* dots[], int L, long long &feature_cnt) {
    int l = 0;
    int a = -1;
    double* LB = new double[k];
    std::fill(LB, LB+k, 0.0);
    double UB_min = std::numeric_limits<double>::max();
    double UB;
    int *mask = new int[k];
    std::fill(mask, mask+k, 1);

    int mask_sum = k;

    while (l <= L && mask_sum > 1) {
        for (int j = 0; j < k; j++) {
            if (mask[j] != 1) continue;
            //if (a == j) continue; 

            if (UB_min < LB[j]) {
                mask[j] = 0;
            } else {
                
                //DistToLevel(int x, int c, int d, double data[], double centroids[], double* data_ss[], double* centroid_ss[], double* dots[], int l, int L, double &UB, double &LB)
                DistToLevel(x, j, d, data, centroids, data_ss, centroid_ss, l, L, dots, UB, LB[j], feature_cnt);
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
    

    return a;
}*/


void Calculate_squared(int d, int elements, double raw[], double* squared[]) {
    int L = log10(d)/log10(4);
    int d_sqrt = sqrt(d);
    int two_p_level_m1, two_p_level;

    for (int e = 0; e < elements; e++) {
        squared[e][0] = raw[e*d+0] * raw[e*d+0];

        for (int l = 1; l <= L; l++){
            squared[e][l] = squared[e][l-1];
            //dots saved from previous level, hence only add dots from this level.
            //adding new cols from from known rows
            for (int l_ = 0; l_ < pow(2,l-1); l_++) {
                for (int l_2 = pow(2,l-1); l_2 < pow(2, l); l_2++) {
                    squared[e][l] += raw[e*d+l_*d_sqrt+l_2]*raw[e*d+l_*d_sqrt+l_2]; 
                }
            }
            //TODO: sqrt stuff for d != 2^x
            //add full new rows
            for (int l_ = pow(2,l-1); l_ < pow(2,l); l_++) {
                for (int l_2 = 0; l_2 < pow(2, l); l_2++) {
                    squared[e][l] += raw[e*d+l_*d_sqrt+l_2]*raw[e*d+l_*d_sqrt+l_2]; 
                }
            }
      
        } 
    }
}

#endif