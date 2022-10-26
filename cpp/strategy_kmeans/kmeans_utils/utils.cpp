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

double Euclidian_distance(int x, int c, int d, int k, double data[], double centroids[]) {
    //x datapoint index
    //c centroid index

    //d features/dimension
    double dist = 0;
    for (int j = 0; j < d; j++) {
        dist += (data[x*d+j]-centroids[c*d+j])*(data[x*d+j]-centroids[c*d+j]);
    };
    //std::cout << dist << std::endl;
    if(dist < 0.0) dist = 0.0;
    return sqrt(dist);
}

void Update_bounds(double data[], double centroids[], double* c_to_c[], double* centroids_ss[], double* l_elkan[], double u_elkan[], double l_hamerly[], int labels[], double div[], double near[], int n, int k, int d) {
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
    //Nearest other centroid
    for (int i = 0; i < k; i++) {
        double smallest = DBL_MAX;
        for (int j = 0; j < k; j++) {
            if (i == j) continue;
            if (c_to_c[i][j] < smallest) {  
                smallest = c_to_c[i][j];  
                near[i] = 0.5 * smallest; 
            }
        }
    }

    //END: Updated l_elkan, u_elkan, l_hamerly, near, c_to_c
};


bool Recalculate(double data[], double centroids[], double old_centroids[], double cluster_count[], int labels[], double div[], int n, int k, int d) {
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

    for (int i = 0; i < k; i++) {      
        std::cout << cluster_count[i] << " "; 
        }
    std::cout << std::endl;

    //calculate div
    for (int j = 0; j < k; j++) {
        double tmp = 0.0;
        for (int f = 0; f < d; f++) {
            tmp += ((centroids[j*d+f] - old_centroids[j*d+f]) *
                    (centroids[j*d+f] - old_centroids[j*d+f]));
        }
        if(tmp < 0.0) tmp = 0.0;
        div[j] = sqrt(tmp);
        if (div[j] > 0) { //This convergence check is slow and can possibly be done better using the old check
            converged = false;
        }
    }

    //END: div[.] updated
    return converged;
}


std::tuple<double, double> DistToLevel(int x, int c, int d, double data[], double centroids[], double* data_ss[], double* centroid_ss[], double* dots[], int l, int L, double &UB, double &LB) {
    //Calculate dots  
    //TODO: is this dot correct or should it loop to l instead?
    //dots[x][c] = 0;
    dots[x][c] = 0;
    int d_sqrt = sqrt(d);
    for (int l_ = 0; l_ < std::min(d_sqrt,(int) pow(2,l)); l_++) {
        //TODO: min with d_sqrt in both loops
        for (int l_2 = 0; l_2 < std::min(d_sqrt,(int) pow(2,l)); l_2++) {
            dots[x][c] += data[x*d+l_*d_sqrt+l_2]*centroids[c*d+l_*d_sqrt+l_2];
        }
    }
    //dots[x][c] += data[x*d+l] * centroids[c*d+l];
    
    
    
    double dist = data_ss[x][L] + centroid_ss[c][L] - 2*dots[x][c]; 

    double margin = 2 * sqrt(data_ss[x][L]-data_ss[x][l]) * sqrt(centroid_ss[c][L]-centroid_ss[c][l]);


    
    LB = sqrt(std::max(0.0,dist - margin));
    UB = sqrt(std::max(0.0,dist + margin));

    return {LB, UB};
};

//TODO: This parsing is insane can can be cleaned clean up, by defining this and dist to level within a strategy scope. 
void MG_SetLabel(int x, int d, int k, double data[],  double centroids[], double* data_ss[], double* centroid_ss[], double* dots[], int L, int labels[], double* l_elkan[], double u_elkan[], double* c_to_c[]) {
    int l = 0;
    int *mask = new int[k];
    std::fill_n(mask, k, 1);
    /*for (int i = 0; i < k; i++) {
        mask[i] = 1;
    }*/

    int mask_sum = k;

    while (l <= L && mask_sum > 1) {
        for (int j = 0; j < k; j++) {
            if (mask[j] != 1) continue;  
            //if (j == labels[x]) continue; //TODO: how to treat the assigned point? is this correct?
            
            //Elkan prune
            double test = 0.5 * c_to_c[labels[x]][j];
            double test2 = l_elkan[x][j];
            double val = std::max(test2, test);// l_elkan[x][j] < 0.5 * c_to_c[labels[x]][j] ? 0.5 * c_to_c[labels[x]][j] : l_elkan[x][j];  
            if (u_elkan[x] < val) {     //Elkan check
                mask[j] = 0;            //Mark as pruned centroid
            } else {
                //DistToLevel params (int x, int c, int d, double data[], double centroids[], double* data_ss[], double* centroid_ss[], double* dots[], int l, int L)
                double UB, LB;
                DistToLevel(x, j, d, data, centroids, data_ss, centroid_ss, dots, l, L, UB, LB);
                //auto val_ = Euclidian_distance(x,j,d,k,data,centroids);
                //UB = val_;
                //LB = val_;
                if (LB > l_elkan[x][j]) {
                    l_elkan[x][j] = LB; //Keep maximum LB per c
                }
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
}

int SetLabel(int x, int d, int k, double data[], double centroids[], double* data_ss[], double* centroid_ss[], double* dots[], int L) {
    int l = 0;
    int a = -1;
    double* LB = new double[k];
    std::fill(LB, LB+k, 0.0);
    double UB_min = std::numeric_limits<double>::max();
    
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
                double UB;
                //DistToLevel(int x, int c, int d, double data[], double centroids[], double* data_ss[], double* centroid_ss[], double* dots[], int l, int L, double &UB, double &LB)
                DistToLevel(x, j, d, data, centroids, data_ss, centroid_ss, dots, l, L, UB, LB[j]);
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
    if (a == -1) {
        
        std::cout << "SetLabel did not set a for item " << x << " " << mask_sum <<std::endl;  
    }
    return a;
}

/*void Calculate_squared(int d, int elements, double raw[], double* squared[]) {

}*/

void Calculate_squared(int d, int elements, double raw[], double* squared[]) {
    int L = log10(d)/log10(4);
    int d_sqrt = sqrt(d);
    int two_p_level_m1, two_p_level;

    for (int e = 0; e < elements; e++) {
        squared[e][0] = raw[e*d+0] * raw[e*d+0];
        for (int l = 1; l <= L; l++){
            squared[e][l] = 0;
            for (int l_ = 0; l_ < pow(2,l); l_++) {
                for (int l_2 = 0; l_2 < pow(2,l); l_2++) {
                    squared[e][l] += raw[e*d+l_*d_sqrt+l_2]*raw[e*d+l_*d_sqrt+l_2];
                    //squared[0][0] = [e*d+0*d_sqrt+0]^2 (0,0)^2
                    //squared[0][1] = [e*d+0*d_sqrt+1]^2 (0,1)^2
                    //squared[1][0] = [e*d+1*d_sqrt+0]^2 (1,8)^2
                    //squared[1][1] = [e*d+1*d_sqrt+1]^2 (1,9)^2
                }
            }// x_sq[e][(e*d,0)^2,(e*d,0)^2+(e*d,1)^2+(e*d+d_sqrt,0)^2+(d_sqrt,1)^2,....]
        } 
    }
}

/*void Calculate_squared_v1(int d, int elements, double raw[], double* squared[]) {
    int L = log10(d)/log10(4);
    int d_sqrt = sqrt(d);
    int two_p_level_m1, two_p_level;

    
    for (int e = 0; e < elements; e++){
        squared[e][L+1] = 0;
        
        for (int l = L; l > 0; l--) {
            squared[e][l] = squared[e][l+1];

            two_p_level_m1 = int(pow(2, l - 1));
            two_p_level = std::min(int(pow(2, l)), d_sqrt);

            for (int i = 0; i < two_p_level_m1; i++) {
                for (int j = i * d_sqrt + two_p_level_m1; j < i * d_sqrt + two_p_level; j++)
                squared[e][l] += (raw[e*d+j] * raw[e*d+j]);
            }
            for (int i = 0; i < std::min(two_p_level_m1, d_sqrt - two_p_level_m1); i++) {
                for (int j = d_sqrt * (two_p_level_m1+i); j < d_sqrt * (two_p_level_m1 + i) + two_p_level; j++)
                squared[e][l] += (raw[e*d+j] * raw[e*d+j]);
            }
        }
        squared[e][0] = squared[e][1] + (raw[e*d+0] * raw[e*d+0]);
    }
}*/
#endif