#include <cfloat>
#include <cmath>
#include <tuple>
#include <algorithm>

double Euclidian_distance(int x, int c, int d, int k, double data[], double centroids[]) {
    //x datapoint index
    //c centroid index

    //d features/dimension
    double dist = 0;
    for (int j = 0; j < d; j++) {
        dist += (data[x*d+j]-centroids[c*d+j])*(data[x*d+j]-centroids[c*d+j]);
    };
    //std::cout << dist << std::endl;
    return dist;
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
            if ((labels[i] != j) && (l_elkan[i][j] < l_elkan[i][smallest_id])) {
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
            double tmp = centroids_ss[i][0] + centroids_ss[j][0];
            for (int f = 0; f < d; f++) {
                tmp -= (2 * centroids[i*d+f] * centroids[j*d+f]);
            }
            if(tmp < 0.0) tmp = 0.0;
            tmp = sqrt(tmp);
            //We can save distances for later use
            c_to_c[i][j] = sqrt(tmp);
            // THEY'RE THE SAME
            c_to_c[j][i] = c_to_c[j][i];
        }
    }
    //Nearest other centroid
    for (int i = 0; i < k; i++) {
        double smallest = DBL_MAX;
        for (int j = 0; j < k; i++) {
            if (i == j) continue;
            if (c_to_c[i][j] < smallest) {   
                near[i] = smallest; 
            }
        }
    }

    //END: Updated l_elkan, u_elkan, l_hamerly, near, c_to_c
};


//void calculate_squared_centroids()

//void calculate_squared_data()

std::tuple<double, double> DistToLevel(int x, int c, int d, double data[], double centroids[], double* data_ss[], double* centroid_ss[], double* dots[], int l, int L) {
    //Calculate dots  
    //TODO: init dots to 0 when defined
    dots[x][c] += data[x*d+l] * centroids[c*d+l];
    double dist = data_ss[x][l] + centroid_ss[c][l] - 2*dots[x][c]; 
    double margin = 2 * sqrt(data_ss[x][L]-data_ss[x][l]) * sqrt(centroid_ss[x][L]-centroid_ss[x][l]);

    double LB = dist - margin;
    double UB = dist + margin;

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
            if (j == labels[x]) continue;
            
            //Elkan prune
            double val = l_elkan[x][j] < 0.5 * c_to_c[labels[x]][j] ? 0.5 * c_to_c[labels[x]][j] : l_elkan[x][j];  
            if (u_elkan[x] < val) {     //Elkan check
                mask[j] = 0;            //Mark as pruned centroid
            } else {
                //DistToLevel params (int x, int c, int d, double data[], double centroids[], double* data_ss[], double* centroid_ss[], double* dots[], int l, int L)
                auto [UB, LB] = DistToLevel(x, j, d, data, centroids, data_ss, centroid_ss, dots, l, L);
                
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