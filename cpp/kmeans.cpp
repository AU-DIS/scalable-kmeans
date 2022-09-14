#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <vector>
#include <iostream>
#include <utility>
#include <queue>
#include <string>
#include <chrono>
#include <string.h>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <limits.h>
#include <float.h>
// #include <cblas.h>



#ifdef DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 0
#endif
#define DEBUGPRINT(fmt, ...) \
            do { if (DEBUG_TEST) fprintf(stderr, fmt, __VA_ARGS__); } while (0)


// #define D 784 // data dimensionality
// #define K 10 // k variable in kmeans
// #define N 69997 // count of data points
#define D 144 // data dimensionality
#define K 10 // k variable in kmeans
#define N 167 // count of data points
#define MAX_ITERATIONS 50 // maximum number of iterations to do before giving up convergence

// double data_arr[N][D]; // the data itself 
// int labels[N]; // where the labels of the algorithm will be eventually
// double centroids[K][D]; // centroid of clusters
// double distances[N][K]; // distance of data[i] to centroid[j]

double **data_arr;
int *labels;
double **centroids;
double **distances;

// STEPWISE EXTRAS
double **data_arr_ss;
double **centroids_ss;
bool **is_candidate;
double **upper_bounds;
double **lower_bounds;
double **incremental_dots;
// STEPWISE EXTRAS end


// HAMERLY EXTRAS
double *hamerly_upper_bounds;
double *hamerly_lower_bounds;
double *closest_centroid_distance;
double *centroid_movement;
// HAMERLY EXTRAS end

// ELKAN EXTRAS
double **elkan_lower_bounds;
double **centroid_to_centroid_distances;
// ELKAN EXTRAS end


// helper variables
double **old_centroids;
int *assigned;

// returns euc_dist of data[i] and centroid[j]
double euclidean_distance(int i, int j) {
    double dist = 0;
    for (int index = 0; index < D; index++) {
        dist += ((data_arr[i][index] - centroids[j][index]) * (data_arr[i][index] - centroids[j][index]));
    }
    return dist;
}

// lloyd, simple
void kmeans() {
    std::cout << "in kmeans..." << std::endl;
    int folan = 0;
    int filan = 0; // for loop usage
    // set initial centroids
    // TODO
    // FIX: for now let's set them to the first k points --> these are a really bad init case, gonna try the random selection 
    for(folan = 0; folan < K; folan++){
        for(filan = 0; filan < D; filan++){
            centroids[folan][filan] = data_arr[folan][filan];
        }
    }
    // TODO : check if I can do this with memcpy, DONE but not tested
    // memcpy(centroids, data_arr, sizeof(double) * K * D);
    std::cout << "copied init centroids" << std::endl;
    // srand(time(0)); 
    // int rand_index = 0;
    // for(folan = 0; folan < K; folan++){
    //     // generate random index
    //     rand_index = rand() % N;
    //     for(filan = 0; filan < D; filan++){
    //         centroids[folan][filan] = data_arr[rand_index][filan];
    //     }
    // }
    // std::cout << "filled centroids randomly for init" << std::endl;



    bool has_converged = false;
    int cluster_counts[K];

    // loop over max_iter
    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
        std::cout << "iteration " << iter << "..." << std::endl;
        // assign points to closest centroid
        // TODO : DONE
        for (folan = 0; folan < N; folan++) {
            for (filan = 0; filan < K; filan++) {
                distances[folan][filan] = euclidean_distance(folan, filan);
            }
        }
        std::cout << "calculated distances" << std::endl;
        // I don't think we need to do this
        // memset(labels, 0, sizeof(labels));
        // set converged to true, innocent until proven guilty:)
        // has_converged = true;
        for (folan = 0; folan < N; folan++) {
            for (filan = 0; filan < K; filan++) {
                if (distances[folan][filan] < distances[folan][labels[folan]]) {
                    labels[folan] = filan;
                    // if labels changes, that means it has not converged
                    // has_converged = false;
                    // this isn't very smart tho, cause it'll write to this var way too many times...
                    // checking K centroids later is much simpler
                }
            }
        }

        std::cout << "set labels" << std::endl;


        // calc new centroids
        // TODO: DONE
        // make copy of old centroids
        // apparently this also does not work, copies pointer somehow I think...
        // memcpy(old_centroids, centroids, sizeof(double) * K * D);
        for (folan = 0; folan < K; folan++) {
            for (filan = 0; filan < D; filan++) {
                old_centroids[folan][filan] = centroids[folan][filan];
            }
        }
        std::cout << "copied centroids to old centroids" << std::endl;
        // set centroids to 0
        memset(cluster_counts, 0, sizeof(int) * K);
        std::cout << "set cluster counts to 0" << std::endl;
        for (folan = 0; folan < K; folan++) {
            // just testing
            // cluster_counts[folan] = 0;
            for (filan = 0; filan < D; filan++) {
                centroids[folan][filan] = 0.0;
            }
        }
        // This doesn't work on doubles
        // memset(centroids, 0, sizeof(double) * K * D);
        std::cout << "after all the memcpys" << std::endl;

        for (folan = 0; folan < N; folan++) {
            cluster_counts[labels[folan]]++;
            for (filan = 0; filan < D; filan++) {
                centroids[labels[folan]][filan] += data_arr[folan][filan];
            }
        }
        for (folan = 0; folan < K; folan++) {
            // to deal with empty clusters
            // if the cluster is empty, keep the old centroid
            if(cluster_counts[folan] > 0){
                for (filan = 0; filan < D; filan++) {
                    centroids[folan][filan] /= cluster_counts[folan];
                }
            } else{
                std::cout << "EMPTYYYYYYYYYYYYYYYYYYYY SOUUUUULLLLLLLLLLLLL" << std::endl;
                std::cout << "EMPTYYYYYYYYYYYYYYYYYYYY SOUUUUULLLLLLLLLLLLL" << std::endl;
                std::cout << "EMPTYYYYYYYYYYYYYYYYYYYY SOUUUUULLLLLLLLLLLLL" << std::endl;
                for (filan = 0; filan < D; filan++) {
                    centroids[folan][filan] = old_centroids[folan][filan];
                }
            }
        }
        std::cout << "calculated new centroids" << std::endl;
        // just to check
        int sanity_check = 0;
        std::cout << "cluster counts..." << std::endl;
        for (folan = 0; folan < K; folan++) {
            sanity_check += cluster_counts[folan];
            std::cout << cluster_counts[folan] << " ";
        }
        std::cout << sanity_check << std::endl;
        // check convergence
        // TODO: gonna do it in labels assignment, changed my mind will do it here, DONE
        has_converged = true;
        for (folan = 0; folan < K; folan++) {
            for (filan = 0; filan < D; filan++) {
                if (old_centroids[folan][filan] != centroids[folan][filan]) {
                    has_converged = false;
                    break;
                }
            }
        }
        std::cout << "checked convergence" << std::endl;

        // end if converged
        if (has_converged) break;
    }
}


// in this version I will calculatet the dists a bit differently
// ||x-c|| ^2 = ||x||^2 + ||c||^2 - 2X.CT
// the ||x||^2 part doesn't matter in the comparisons I will do, so I'll skip them
// the ||c||^2 I will calculate only once and use all the time
// the x.cT I will first do simply with a for loop
// TESTED
void kmeans_v2() {
    std::cout << "in kmeans_v2 ..." << std::endl;
    int folan = 0;
    int filan = 0;
    int ashghal = 0;// for loop usage
    // set initial centroids
    // TODO
    // FIX: for now let's set them to the first k points --> these are a really bad init case, gonna try the random selection 
    // for(folan = 0; folan < K; folan++){
    //     for(filan = 0; filan < D; filan++){
    //         centroids[folan][filan] = data_arr[folan][filan];
    //     }
    // }
    // TODO : check if I can do this with memcpy, DONE but not tested
    memcpy(centroids, data_arr, sizeof(double) * K * D);
    std::cout << "copied init centroids" << std::endl;
    // srand(time(0)); 
    // int rand_index = 0;
    // for(folan = 0; folan < K; folan++){
    //     // generate random index
    //     rand_index = rand() % N;
    //     for(filan = 0; filan < D; filan++){
    //         centroids[folan][filan] = data_arr[rand_index][filan];
    //     }
    // }
    // std::cout << "filled centroids randomly for init" << std::endl;



    bool has_converged = false;
    int cluster_counts[K];
    // V2 DIFF start
    double centroid_squares[K];
    // V2 DIFF end

    // loop over max_iter
    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
        std::cout << "iteration " << iter << "..." << std::endl;
        // assign points to closest centroid
        // TODO : DONE
        // V2 DIFF start
        // calculate the square sum of centroids

        for (folan = 0; folan < K; folan++) {
            centroid_squares[folan] = 0;
            for (filan = 0; filan < D; filan++) {
                centroid_squares[folan] += (centroids[folan][filan] * centroids[folan][filan]);
            }
        }
        // V2 DIFF end
        for (folan = 0; folan < N; folan++) {
            for (filan = 0; filan < K; filan++) {
                // V2 DIFF START
                distances[folan][filan] = centroid_squares[filan];
                // the dot product
                for (ashghal = 0; ashghal < D; ashghal++) {
                    distances[folan][filan] -= (2 * data_arr[folan][ashghal] * centroids[filan][ashghal]);
                }
                // V2 DIFF end
            }
        }
        std::cout << "calculated distances" << std::endl;

        // I don't think we need to do this
        // memset(labels, 0, sizeof(labels));
        // set converged to true, innocent until proven guilty:)
        // has_converged = true;
        for (folan = 0; folan < N; folan++) {
            for (filan = 0; filan < K; filan++) {
                if (distances[folan][filan] < distances[folan][labels[folan]]) {
                    labels[folan] = filan;
                    // if labels changes, that means it has not converged
                    // has_converged = false;
                    // this isn't very smart tho, cause it'll write to this var way too many times...
                    // checking K centroids later is much simpler
                }
            }
        }

        std::cout << "set labels" << std::endl;


        // calc new centroids
        // TODO: DONE
        // make copy of old centroids
        // apparently this also does not work, copies pointer somehow I think...
        // memcpy(old_centroids, centroids, sizeof(double) * K * D);
        for (folan = 0; folan < K; folan++) {
            for (filan = 0; filan < D; filan++) {
                old_centroids[folan][filan] = centroids[folan][filan];
            }
        }
        std::cout << "copied centroids to old centroids" << std::endl;
        // set centroids to 0
        memset(cluster_counts, 0, sizeof(int) * K);
        std::cout << "set cluster counts to 0" << std::endl;
        for (folan = 0; folan < K; folan++) {
            // just testing
            // cluster_counts[folan] = 0;
            for (filan = 0; filan < D; filan++) {
                centroids[folan][filan] = 0.0;
            }
        }
        // This doesn't work on doubles
        // memset(centroids, 0, sizeof(double) * K * D);
        std::cout << "after all the memcpys" << std::endl;

        for (folan = 0; folan < N; folan++) {
            cluster_counts[labels[folan]]++;
            for (filan = 0; filan < D; filan++) {
                centroids[labels[folan]][filan] += data_arr[folan][filan];
            }
        }
        for (folan = 0; folan < K; folan++) {
            // to deal with empty clusters
            // if the cluster is empty, keep the old centroid
            if(cluster_counts[folan] > 0){
                for (filan = 0; filan < D; filan++) {
                    centroids[folan][filan] /= cluster_counts[folan];
                }
            } else{
                for (filan = 0; filan < D; filan++) {
                    centroids[folan][filan] = old_centroids[folan][filan];
                }
            }
        }
        std::cout << "calculated new centroids" << std::endl;
        // just to check
        int sanity_check = 0;
        std::cout << "cluster counts..." << std::endl;
        for (folan = 0; folan < K; folan++) {
            sanity_check += cluster_counts[folan];
            std::cout << cluster_counts[folan] << " ";
        }
        std::cout << sanity_check << std::endl;
        // check convergence
        // TODO: gonna do it in labels assignment, changed my mind will do it here, DONE
        has_converged = true;
        for (folan = 0; folan < K; folan++) {
            for (filan = 0; filan < D; filan++) {
                if (old_centroids[folan][filan] != centroids[folan][filan]) {
                    has_converged = false;
                    break;
                }
            }
        }
        std::cout << "checked convergence" << std::endl;

        // end if converged
        if (has_converged) break;
    }
}


// STEPWISE EXTRAS
// will fill up the data_arr_ss correctly
// data_ss[i][0] = total_sum_of_sqaures
// data_ss[i][1] = total_sum_of_sqaures - 1x1square_ss
// data_ss[i][2] = total_sum_of_sqaures - 2x2sqaure_ss
// ...
// data_ss[i][-1] = 0
// TODO: done not tested
void calculate_data_square_sums() {
    int ashghal, halghe;
    int d_sqrt = sqrt(D);
    int two_p_level_m1, two_p_level;


    for (int i = 0; i < N; i++) {
        // last one is 0
        data_arr_ss[i][int(log2(sqrt(D))) + 2] = 0;
        // fill from back to front
        // every time, we will just add the extra cells needed to complete it
        // the loop ranges are from the inc_dot calculations, I trusted them, haven't checked again
        for(int j = int(log2(sqrt(D))) + 1; j > 0; j--){

            data_arr_ss[i][j] = data_arr_ss[i][j+1];

            two_p_level_m1 = int(pow(2, j - 1));
            two_p_level = std::min(int(pow(2, j)), d_sqrt);

            for(ashghal = 0; ashghal < two_p_level_m1; ashghal++){
                // amoodi ha
                for(halghe = ashghal * d_sqrt + two_p_level_m1; halghe < ashghal * d_sqrt + two_p_level; halghe++){
                    data_arr_ss[i][j] += (data_arr[i][halghe] * data_arr[i][halghe]);
                }
            }
            for(ashghal = 0; ashghal < std::min(two_p_level_m1, d_sqrt - two_p_level_m1); ashghal++){
                // ofoghi ha
                for(halghe = d_sqrt * (two_p_level_m1 + ashghal); halghe < d_sqrt * (two_p_level_m1 + ashghal) + two_p_level; halghe++){
                    data_arr_ss[i][j] += (data_arr[i][halghe] * data_arr[i][halghe]);
                }
            }
        }
        data_arr_ss[i][0] = data_arr_ss[i][1] +  (data_arr[i][0] * data_arr[i][0]);

    }

}

// will fill up the centroids_ss correctly
// TODO: done not tested
void calculate_centroids_square_sums() {
    int ashghal, halghe;
    int d_sqrt = sqrt(D);
    int two_p_level_m1, two_p_level;


    for (int i = 0; i < K; i++) {
        // last one is 0
        centroids_ss[i][int(log2(sqrt(D))) + 2] = 0;
        // fill from back to front
        // every time, we will just add the extra cells needed to complete it
        // the loop ranges are from the inc_dot calculations, I trusted them, haven't checked again
        for(int j = int(log2(sqrt(D))) + 1; j > 0; j--){

            centroids_ss[i][j] = centroids_ss[i][j+1];

            two_p_level_m1 = int(pow(2, j - 1));
            two_p_level = std::min(int(pow(2, j)), d_sqrt);

            for(ashghal = 0; ashghal < two_p_level_m1; ashghal++){
                // amoodi ha
                for(halghe = ashghal * d_sqrt + two_p_level_m1; halghe < ashghal * d_sqrt + two_p_level; halghe++){
                    centroids_ss[i][j] += (centroids[i][halghe] * centroids[i][halghe]);
                }
            }
            for(ashghal = 0; ashghal < std::min(two_p_level_m1, d_sqrt - two_p_level_m1); ashghal++){
                // ofoghi ha
                for(halghe = d_sqrt * (two_p_level_m1 + ashghal); halghe < d_sqrt * (two_p_level_m1 + ashghal) + two_p_level; halghe++){
                    centroids_ss[i][j] += (centroids[i][halghe] * centroids[i][halghe]);
                }
            }
        }
        centroids_ss[i][0] = centroids_ss[i][1] + (centroids[i][0] * centroids[i][0]);

    }

}


// for all the data_points
// TODO: done, not tested
void calculate_distances_till_level(int level) {
    int folan, filan, ashghal, halghe;
    int d_sqrt = sqrt(D);
    int two_p_level_m1 = int(pow(2, level - 1)), two_p_level = std::min(int(pow(2, level)), d_sqrt);
    double this_dot = 0.0;

    for (folan = 0; folan < N; folan++) {

        if (labels[folan] > 0) continue;

        for (filan = 0; filan < K; filan++) {

            if (is_candidate[folan][filan] == false) continue;

            // known parts

            upper_bounds[folan][filan] = centroids_ss[filan][0];
            lower_bounds[folan][filan] = centroids_ss[filan][0];

            // debug
            if (level == 0) {
                incremental_dots[folan][filan] = data_arr[folan][0] * centroids[filan][0];
            } else {
                for(ashghal = 0; ashghal < two_p_level_m1; ashghal++){
                    // amoodi ha
                    for(halghe = ashghal * d_sqrt + two_p_level_m1; halghe < ashghal * d_sqrt + two_p_level; halghe++){
                        incremental_dots[folan][filan] += (data_arr[folan][halghe] * centroids[filan][halghe]);
                    }
                }
                for(ashghal = 0; ashghal < std::min(two_p_level_m1, d_sqrt - two_p_level_m1); ashghal++){
                    // ofoghi ha
                    for(halghe = d_sqrt * (two_p_level_m1 + ashghal); halghe < d_sqrt * (two_p_level_m1 + ashghal) + two_p_level; halghe++){
                        incremental_dots[folan][filan] += (data_arr[folan][halghe] * centroids[filan][halghe]);
                    }
                }
            }
            upper_bounds[folan][filan] -= 2 * incremental_dots[folan][filan];
            lower_bounds[folan][filan] -= 2 * incremental_dots[folan][filan];

            // this_dot = 0.0;
            // for(int i = 0; i < two_p_level; i++){
            //     for(int j = 0; j < two_p_level; j++){
            //         this_dot += (data_arr[folan][d_sqrt*i + j] * centroids[filan][d_sqrt*i + j]);
            //     }
            // }
            // upper_bounds[folan][filan] -= (2 * this_dot);
            // lower_bounds[folan][filan] -= (2 * this_dot);

            // debug



            // unknown parts
            this_dot = sqrt(data_arr_ss[folan][level + 1] * centroids_ss[filan][level + 1]);
            upper_bounds[folan][filan] += (2 * this_dot);
            lower_bounds[folan][filan] -= (2 * this_dot);


        }


    }

}


// copied from calculate_distances_till_level
// just has the x^2 and is sqrt of whatever it was
void calculate_sqrt_distances_till_level(int level) {
    int folan, filan, ashghal, halghe;
    int d_sqrt = sqrt(D);
    int two_p_level_m1 = int(pow(2, level - 1)), two_p_level = std::min(int(pow(2, level)), d_sqrt);
    double this_dot = 0.0;
    double tmp_ub, tmp_lb;

    for (folan = 0; folan < N; folan++) {

        if (labels[folan] > 0) continue;

        for (filan = 0; filan < K; filan++) {

            if (is_candidate[folan][filan] == false) continue;

            // known parts

            tmp_ub = centroids_ss[filan][0] + data_arr_ss[folan][0];
            tmp_lb = centroids_ss[filan][0] + data_arr_ss[folan][0];

            // debug
            if (level == 0) {
                incremental_dots[folan][filan] = data_arr[folan][0] * centroids[filan][0];
            } else {
                for(ashghal = 0; ashghal < two_p_level_m1; ashghal++){
                    // amoodi ha
                    for(halghe = ashghal * d_sqrt + two_p_level_m1; halghe < ashghal * d_sqrt + two_p_level; halghe++){
                        incremental_dots[folan][filan] += (data_arr[folan][halghe] * centroids[filan][halghe]);
                    }
                }
                for(ashghal = 0; ashghal < std::min(two_p_level_m1, d_sqrt - two_p_level_m1); ashghal++){
                    // ofoghi ha
                    for(halghe = d_sqrt * (two_p_level_m1 + ashghal); halghe < d_sqrt * (two_p_level_m1 + ashghal) + two_p_level; halghe++){
                        incremental_dots[folan][filan] += (data_arr[folan][halghe] * centroids[filan][halghe]);
                    }
                }
            }
            tmp_ub -= 2 * incremental_dots[folan][filan];
            tmp_lb -= 2 * incremental_dots[folan][filan];

            // this_dot = 0.0;
            // for(int i = 0; i < two_p_level; i++){
            //     for(int j = 0; j < two_p_level; j++){
            //         this_dot += (data_arr[folan][d_sqrt*i + j] * centroids[filan][d_sqrt*i + j]);
            //     }
            // }
            // upper_bounds[folan][filan] -= (2 * this_dot);
            // lower_bounds[folan][filan] -= (2 * this_dot);

            // debug



            // unknown parts
            this_dot = sqrt(data_arr_ss[folan][level + 1] * centroids_ss[filan][level + 1]);
            tmp_ub += (2 * this_dot);
            tmp_lb -= (2 * this_dot);

            // sqrt the whole thing
            if(tmp_lb < 0.0) tmp_lb = 0.0;
            if(tmp_ub < 0.0) tmp_ub = 0.0;
            upper_bounds[folan][filan] = sqrt(tmp_ub);
            lower_bounds[folan][filan] = sqrt(tmp_lb);

        }


    }

}


// copied from calculate_sqrt_distances_till_level
// just checks assigned instead of labels
void calculate_sqrt_distances_till_level_with_assigned(int level) {
    int folan, filan, ashghal, halghe;
    int d_sqrt = sqrt(D);
    int two_p_level_m1 = int(pow(2, level - 1)), two_p_level = std::min(int(pow(2, level)), d_sqrt);
    double this_dot = 0.0;
    double tmp_ub, tmp_lb;

    for (folan = 0; folan < N; folan++) {

        if (assigned[folan] > 0) continue;

        for (filan = 0; filan < K; filan++) {

            if (is_candidate[folan][filan] == false) continue;

            // known parts

            tmp_ub = centroids_ss[filan][0] + data_arr_ss[folan][0];
            tmp_lb = centroids_ss[filan][0] + data_arr_ss[folan][0];

            // debug
            if (level == 0) {
                incremental_dots[folan][filan] = data_arr[folan][0] * centroids[filan][0];
            } else {
                for(ashghal = 0; ashghal < two_p_level_m1; ashghal++){
                    // amoodi ha
                    for(halghe = ashghal * d_sqrt + two_p_level_m1; halghe < ashghal * d_sqrt + two_p_level; halghe++){
                        incremental_dots[folan][filan] += (data_arr[folan][halghe] * centroids[filan][halghe]);
                    }
                }
                for(ashghal = 0; ashghal < std::min(two_p_level_m1, d_sqrt - two_p_level_m1); ashghal++){
                    // ofoghi ha
                    for(halghe = d_sqrt * (two_p_level_m1 + ashghal); halghe < d_sqrt * (two_p_level_m1 + ashghal) + two_p_level; halghe++){
                        incremental_dots[folan][filan] += (data_arr[folan][halghe] * centroids[filan][halghe]);
                    }
                }
            }
            tmp_ub -= 2 * incremental_dots[folan][filan];
            tmp_lb -= 2 * incremental_dots[folan][filan];

            // this_dot = 0.0;
            // for(int i = 0; i < two_p_level; i++){
            //     for(int j = 0; j < two_p_level; j++){
            //         this_dot += (data_arr[folan][d_sqrt*i + j] * centroids[filan][d_sqrt*i + j]);
            //     }
            // }
            // upper_bounds[folan][filan] -= (2 * this_dot);
            // lower_bounds[folan][filan] -= (2 * this_dot);

            // debug



            // unknown parts
            this_dot = sqrt(data_arr_ss[folan][level + 1] * centroids_ss[filan][level + 1]);
            tmp_ub += (2 * this_dot);
            tmp_lb -= (2 * this_dot);

            // sqrt the whole thing
            if(tmp_lb < 0.0) tmp_lb = 0.0;
            if(tmp_ub < 0.0) tmp_ub = 0.0;
            upper_bounds[folan][filan] = sqrt(tmp_ub);
            lower_bounds[folan][filan] = sqrt(tmp_lb);

        }


    }

}


// for one point in particular
// TODO
void calculate_distances_till_level(int point_index, int level) {}

// has the while(level < log(D))
// starts with level = 0
// TODO: done not tested
void calculate_labels() {

    int folan, filan, ashghal, alaki;
    // I really thought this would work:(, sadly it doesn't
    // memset(labels, 1, sizeof(int) * N);
    // sanity check
    // for(folan = 0; folan < N; folan++){
    //     if(labels[folan] != -1){
    //         std::cout << "TERROR! DISASTER! WE WERE DECIEVED:(" << std::endl;
    //     }
    // }
    // end of sanity check
    for (folan = 0; folan < N; folan++) labels[folan] = -1;

    std::cout << "after memset labels" << std::endl;


    int level = 0;
    // int level = int(log2(int(sqrt(D))) +  1);
    // I forgot that you can't do this to doubles
    // will do it in main
    // memset(incremental_dots, 0, sizeof(double) * K * N);

    int *smallest_ub = (int *) malloc(N * sizeof(int));

    std::cout << "after init labels and inc_dots and smallest_ub" << std::endl;

    while (level < int(log2(int(sqrt(D))) + 2)) {
        std::cout << "in level while loop, level = " << level << std::endl;
        calculate_distances_till_level(level);
        std::cout << "after calculate dist till level" << std::endl;



        // find the smallest upperbound per point
        memset(smallest_ub, 0, sizeof(int) * N);
        for (folan = 0; folan < N; folan++) {
            for (filan = 0; filan < K; filan++) {
                if (upper_bounds[folan][filan] < upper_bounds[folan][smallest_ub[folan]]) smallest_ub[folan] = filan;
            }
        }

        for (folan = 0; folan < N; folan++) {
            if (labels[folan] > 0) continue;
            if (level == int(log2(int(sqrt(D))) + 1)) labels[folan] = smallest_ub[folan];
            else {
                for (filan = 0; filan < K; filan++) {
                    if (filan == smallest_ub[folan]) continue;
                    if (!is_candidate[folan][filan]) continue;
                    if (lower_bounds[folan][filan] >= upper_bounds[folan][smallest_ub[folan]])
                        is_candidate[folan][filan] = false;
                }

                // check if only one is true
                alaki = 0;
                for (filan = 0; filan < K; filan++) {
                    if (is_candidate[folan][filan]) alaki++;
                }
                // then the only one left is the one with the smallest_ub
                if (alaki == 1) labels[folan] = smallest_ub[folan];
            }
        }

        level++;

    }


}


// copied from calculate_labels
// just calls the sqrt_dist function
void calculate_labels_with_sqrt() {

    int folan, filan, ashghal, alaki;
    // I really thought this would work:(, sadly it doesn't
    // memset(labels, 1, sizeof(int) * N);
    // sanity check
    // for(folan = 0; folan < N; folan++){
    //     if(labels[folan] != -1){
    //         std::cout << "TERROR! DISASTER! WE WERE DECIEVED:(" << std::endl;
    //     }
    // }
    // end of sanity check
    for (folan = 0; folan < N; folan++) labels[folan] = -1;

    std::cout << "after memset labels" << std::endl;


    int level = 0;
    // int level = int(log2(int(sqrt(D))) +  1);
    // I forgot that you can't do this to doubles
    // will do it in main
    // memset(incremental_dots, 0, sizeof(double) * K * N);

    int *smallest_ub = (int *) malloc(N * sizeof(int));

    std::cout << "after init labels and inc_dots and smallest_ub" << std::endl;

    while (level < int(log2(int(sqrt(D))) + 2)) {
        std::cout << "in level while loop, level = " << level << std::endl;
        calculate_sqrt_distances_till_level(level);
        std::cout << "after calculate dist till level" << std::endl;



        // find the smallest upperbound per point
        memset(smallest_ub, 0, sizeof(int) * N);
        for (folan = 0; folan < N; folan++) {
            for (filan = 0; filan < K; filan++) {
                if (upper_bounds[folan][filan] < upper_bounds[folan][smallest_ub[folan]]) smallest_ub[folan] = filan;
            }
        }

        for (folan = 0; folan < N; folan++) {
            if (labels[folan] > 0) continue;
            if (level == int(log2(int(sqrt(D))) + 1)) labels[folan] = smallest_ub[folan];
            else {
                for (filan = 0; filan < K; filan++) {
                    if (filan == smallest_ub[folan]) continue;
                    if (!is_candidate[folan][filan]) continue;
                    if (lower_bounds[folan][filan] >= upper_bounds[folan][smallest_ub[folan]])
                        is_candidate[folan][filan] = false;
                }

                // check if only one is true
                alaki = 0;
                for (filan = 0; filan < K; filan++) {
                    if (is_candidate[folan][filan]) alaki++;
                }
                // then the only one left is the one with the smallest_ub
                if (alaki == 1) labels[folan] = smallest_ub[folan];
            }
        }

        level++;

    }


}


// copied from calculate_labels
// just calls the sqrt_dist function
// fills up the hamerly and elkan ub and lbs
// IM NOT SURE WHAT THIS IS, IT DOESN'T LOOK RIGTH ANYWAY!!
void calculate_labels_with_sqrt_integrated() {

    int folan, filan, ashghal, alaki;
    // I really thought this would work:(, sadly it doesn't
    // memset(labels, 1, sizeof(int) * N);
    // sanity check
    // for(folan = 0; folan < N; folan++){
    //     if(labels[folan] != -1){
    //         std::cout << "TERROR! DISASTER! WE WERE DECIEVED:(" << std::endl;
    //     }
    // }
    // end of sanity check
    for (folan = 0; folan < N; folan++) labels[folan] = -1;

    std::cout << "after memset labels" << std::endl;


    int level = 0;
    // int level = int(log2(int(sqrt(D))) +  1);
    // I forgot that you can't do this to doubles
    // will do it in main
    // memset(incremental_dots, 0, sizeof(double) * K * N);

    int *smallest_ub = (int *) malloc(N * sizeof(int));

    std::cout << "after init labels and inc_dots and smallest_ub" << std::endl;

    while (level < int(log2(int(sqrt(D))) + 2)) {
        std::cout << "in level while loop, level = " << level << std::endl;
        calculate_sqrt_distances_till_level(level);
        std::cout << "after calculate dist till level" << std::endl;



        // find the smallest upperbound per point
        memset(smallest_ub, 0, sizeof(int) * N);
        for (folan = 0; folan < N; folan++) {
            for (filan = 0; filan < K; filan++) {
                if (upper_bounds[folan][filan] < upper_bounds[folan][smallest_ub[folan]]) smallest_ub[folan] = filan;
            }
        }

        for (folan = 0; folan < N; folan++) {
            if (labels[folan] > 0) continue;
            if (level == int(log2(int(sqrt(D))) + 1)) {
                labels[folan] = smallest_ub[folan];
                // V9
                hamerly_upper_bounds[folan] = upper_bounds[folan][smallest_ub[folan]];
                // v9
            } else {
                for (filan = 0; filan < K; filan++) {
                    if (filan == smallest_ub[folan]) continue;
                    if (!is_candidate[folan][filan]) continue;
                    if (lower_bounds[folan][filan] >= upper_bounds[folan][smallest_ub[folan]]) {
                        is_candidate[folan][filan] = false;
                        // V9
                        // this means we will no longer calculate the rest of the levels for this centroid
                        // so we can set the elkan_lb as the last lb calculated
                        elkan_lower_bounds[folan][filan] = lower_bounds[folan][filan];
                        // we should also try to find the second closest centroid
                        // TODO: im not sure if i should clear the old ham_lbs from last round
                        // what if it's masked completely, then we would never get here, so the old one should still be valid.
                        // nvm, i think im right, i won't clear it
                        if (lower_bounds[folan][filan] < hamerly_lower_bounds[folan]) {
                            hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                        }
                        // V9
                    }
                }

                // check if only one is true
                alaki = 0;
                for (filan = 0; filan < K; filan++) {
                    if (is_candidate[folan][filan]) alaki++;
                }
                // then the only one left is the one with the smallest_ub
                if (alaki == 1) {
                    labels[folan] = smallest_ub[folan];
                    // V9
                    hamerly_upper_bounds[folan] = upper_bounds[folan][smallest_ub[folan]];
                    // v9
                }
            }
        }

        level++;

    }


}


// copied from calculate_labels
// just calls the sqrt_dist function
// fills up the hamerly and elkan ub and lbs
void calculate_labels_with_sqrt_hamerly_integrated() {
    std::cout << "in calc_labels_ham_integrated" << std::endl;
    int folan, filan, ashghal, alaki;
    // I really thought this would work:(, sadly it doesn't
    // memset(labels, 1, sizeof(int) * N);
    // sanity check
    // for(folan = 0; folan < N; folan++){
    //     if(labels[folan] != -1){
    //         std::cout << "TERROR! DISASTER! WE WERE DECIEVED:(" << std::endl;
    //     }
    // }
    // end of sanity check
    // V9
    // for(folan = 0; folan < N; folan++) labels[folan] = -1;
    // V9

    // std::cout << "after memset labels" << std::endl;


    int level = 0;
    // int level = int(log2(int(sqrt(D))) +  1);
    // I forgot that you can't do this to doubles
    // will do it in main
    // memset(incremental_dots, 0, sizeof(double) * K * N);

    int *smallest_ub = (int *) malloc(N * sizeof(int));


    // V9
    // I can't set the labels to -1, cause they have to be the last one
    // so to still have that continue condition
    // im gonna make a "assigned" list that gets filled every time
    // int *assigned = (int*)calloc(N,sizeof(int));
    memset(assigned, 0, sizeof(int) * N);
    // V9

    std::cout << "after init assigned and inc_dots and smallest_ub" << std::endl;

    while (level < int(log2(int(sqrt(D))) + 2)) {
        std::cout << "in level while loop, level = " << level << std::endl;
        calculate_sqrt_distances_till_level_with_assigned(level);
        std::cout << "after calculate dist till level" << std::endl;



        // find the smallest upperbound per point
        memset(smallest_ub, 0, sizeof(int) * N);
        for (folan = 0; folan < N; folan++) {
            for (filan = 0; filan < K; filan++) {
                if (upper_bounds[folan][filan] < upper_bounds[folan][smallest_ub[folan]]) smallest_ub[folan] = filan;
            }
        }

        for (folan = 0; folan < N; folan++) {
            // V9
            // if(labels[folan] > 0) continue;
            if (assigned[folan] > 0) continue;
            // V9
            if (level == int(log2(int(sqrt(D))) + 1)) {
                labels[folan] = smallest_ub[folan];
                // V9
                hamerly_upper_bounds[folan] = upper_bounds[folan][smallest_ub[folan]];
                assigned[folan] = 1;
                // v9
            } else {
                for (filan = 0; filan < K; filan++) {
                    if (filan == smallest_ub[folan]) continue;
                    if (!is_candidate[folan][filan]) continue;
                    if (lower_bounds[folan][filan] >= upper_bounds[folan][smallest_ub[folan]]) {
                        is_candidate[folan][filan] = false;
                        // V9
                        // we should also try to find the second closest centroid
                        // TODO: im not sure if i should clear the old ham_lbs from last round
                        // what if it's masked completely, then we would never get here, so the old one should still be valid.
                        // nvm, i think im right, i won't clear it
                        if (lower_bounds[folan][filan] < hamerly_lower_bounds[folan]) {
                            hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                        }
                        // V9
                    }
                }

                // check if only one is true
                alaki = 0;
                for (filan = 0; filan < K; filan++) {
                    if (is_candidate[folan][filan]) alaki++;
                }
                // then the only one left is the one with the smallest_ub
                if (alaki == 1) {
                    labels[folan] = smallest_ub[folan];
                    // V9
                    hamerly_upper_bounds[folan] = upper_bounds[folan][smallest_ub[folan]];
                    assigned[folan] = 1;
                    // v9
                }
            }
        }

        level++;

    }


}



// STEPWISE EXTRAS end


// step wise distance calculation
// copied from v2
// TODO: DONE
void kmeans_v4() {
    std::cout << "in kmeans_v4 ..." << std::endl;
    int folan = 0;
    int filan = 0;
    int ashghal = 0;// for loop usage



    // set initial centroids
    // TODO
    // FIX: for now let's set them to the first k points --> these are a really bad init case, gonna try the random selection 
    for(folan = 0; folan < K; folan++){
        for(filan = 0; filan < D; filan++){
            centroids[folan][filan] = data_arr[folan][filan];
        }
    }
    // TODO : check if I can do this with memcpy, DONE but not tested
    // memcpy(centroids, data_arr, sizeof(double) * K * D);
    calculate_centroids_square_sums();
    std::cout << "copied init centroids, and ss-ed them" << std::endl;


    // srand(time(0)); 
    // int rand_index = 0;
    // for(folan = 0; folan < K; folan++){
    //     // generate random index
    //     rand_index = rand() % N;
    //     for(filan = 0; filan < D; filan++){
    //         centroids[folan][filan] = data_arr[rand_index][filan];
    //     }
    // }
    // std::cout << "filled centroids randomly for init" << std::endl;



    bool has_converged = false;
    int cluster_counts[K];
    // V2 DIFF start
    // double centroid_squares[K];
    // V2 DIFF end

    // loop over max_iter
    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
        std::cout << "iteration " << iter << "..." << std::endl;

        // assign points to closest centroid
        // TODO 
        // V4 DIFF
        // reset the is_candidates
        for (folan = 0; folan < N; folan++) {
            for (filan = 0; filan < K; filan++) {
                is_candidate[folan][filan] = true;
            }
        }
        // calculate_labels();
        // calculate_labels_with_sqrt();
        calculate_labels_with_sqrt_hamerly_integrated();

        // V4 DIFF end
        std::cout << "set labels" << std::endl;


        // calc new centroids
        // TODO: DONE
        // make copy of old centroids
        // apparently this also does not work, copies pointer somehow I think...
        // memcpy(old_centroids, centroids, sizeof(double) * K * D);
        for (folan = 0; folan < K; folan++) {
            for (filan = 0; filan < D; filan++) {
                old_centroids[folan][filan] = centroids[folan][filan];
            }
        }
        std::cout << "copied centroids to old centroids" << std::endl;
        // set centroids to 0
        memset(cluster_counts, 0, sizeof(int) * K);
        std::cout << "set cluster counts to 0" << std::endl;
        for (folan = 0; folan < K; folan++) {
            // just testing
            // cluster_counts[folan] = 0;
            for (filan = 0; filan < D; filan++) {
                centroids[folan][filan] = 0.0;
            }
        }
        // This doesn't work on doubles
        // memset(centroids, 0, sizeof(double) * K * D);
        std::cout << "after all the memcpys" << std::endl;

        for (folan = 0; folan < N; folan++) {
            // std::cout << "debug checkpoint 1, folan = " << folan << std::endl;
            // std::cout << "labels[folan] " << labels[folan] << std::endl;
            cluster_counts[labels[folan]]++;
            // std::cout << "debug checkpoint 2, folan = " << folan << std::endl;
            for (filan = 0; filan < D; filan++) {
                // std::cout << "debug checkpoint 3, folan = " << folan << " filan = " << filan << std::endl;
                centroids[labels[folan]][filan] += data_arr[folan][filan];
                // std::cout << "debug checkpoint 4, folan = " << folan << " filan = " << filan << std::endl;
            }
            // std::cout << "debug checkpoint 5, folan = " << folan << std::endl;
        }
        // std::cout << "debug checkpoint 6" << std::endl;

        for (folan = 0; folan < K; folan++) {
            // to deal with empty clusters
            // if the cluster is empty, keep the old centroid
            if(cluster_counts[folan] > 0){
                for (filan = 0; filan < D; filan++) {
                    centroids[folan][filan] /= cluster_counts[folan];
                }
            } else{
                std::cout << "EMPTYYYYYYYYYYYYYYYYYYYY SOUUUUULLLLLLLLLLLLL" << std::endl;
                std::cout << "EMPTYYYYYYYYYYYYYYYYYYYY SOUUUUULLLLLLLLLLLLL" << std::endl;
                std::cout << "EMPTYYYYYYYYYYYYYYYYYYYY SOUUUUULLLLLLLLLLLLL" << std::endl;
                for (filan = 0; filan < D; filan++) {
                    centroids[folan][filan] = old_centroids[folan][filan];
                }
            }
        }
        std::cout << "calculated new centroids" << std::endl;
        calculate_centroids_square_sums();
        std::cout << "and ss-ed them" << std::endl;
        // just to check
        int sanity_check = 0;
        std::cout << "cluster counts..." << std::endl;
        for (folan = 0; folan < K; folan++) {
            sanity_check += cluster_counts[folan];
            std::cout << cluster_counts[folan] << " ";
        }
        std::cout << sanity_check << std::endl;
        // check convergence
        // TODO: gonna do it in labels assignment, changed my mind will do it here, DONE
        has_converged = true;
        for (folan = 0; folan < K; folan++) {
            for (filan = 0; filan < D; filan++) {
                if (old_centroids[folan][filan] != centroids[folan][filan]) {
                    has_converged = false;
                    break;
                }
            }
        }
        std::cout << "checked convergence" << std::endl;

        // end if converged
        if (has_converged) break;
    }
}


// hamerly
// copied from v2
// TODO: STILL BUGGY, WAITING FOR MOHAMMAD TO FIX HIS CODE SO I CAN CHECK THE UB AND LB
// UPDATE: I think it's the x^2 and sqrt that I removed then the TE doesn't make sense anymore so I have to add them back 
// DONE: holy shit, that was it:)) I just added the sqrt and x^2 everywhere and now it's correct:)
void kmeans_v5() {
    std::cout << "in kmeans_v5 ..." << std::endl;
    int folan = 0;
    int filan = 0;
    int ashghal = 0;// for loop usage
    int furthest_moving_centroid, second_furthest_moving_centroid;
    double smallest, second_smallest;
    double tmp, hamerly_bound;
    // set initial centroids

    // FATEMEH: kasper says the memcpy ruins things:-?

    for(folan = 0; folan < K; folan++){
        for(filan = 0; filan < D; filan++){
            centroids[folan][filan] = data_arr[folan][filan];
        }
    }

    // memcpy(centroids, data_arr, sizeof(double) * K * D);
    std::cout << "copied init centroids" << std::endl;


    bool has_converged = false;
    int cluster_counts[K];
    // V2 DIFF start
    double centroid_squares[K];
    // V2 DIFF end

    // loop over max_iter
    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
        std::cout << "iteration " << iter << "..." << std::endl;

        // calculate the square sum of centroids
        for (folan = 0; folan < K; folan++) {
            centroid_squares[folan] = 0.0;
            for (filan = 0; filan < D; filan++) {
                centroid_squares[folan] += (centroids[folan][filan] * centroids[folan][filan]);
            }
        }

        // calculate closest centroid to each centroid
        for (folan = 0; folan < K; folan++) {
            smallest = DBL_MAX;
            for (filan = 0; filan < K; filan++) {
                tmp = centroid_squares[filan] + centroid_squares[folan];
                for (ashghal = 0; ashghal < D; ashghal++) {
                    tmp -= (2 * centroids[folan][ashghal] * centroids[filan][ashghal]);
                }
                if (tmp < smallest) smallest = tmp;
            }
            if(smallest < 0.0) smallest = 0.0;
            closest_centroid_distance[folan] = sqrt(smallest);
        }
        std::cout << "found closest distance to each centroid" << std::endl;


        if (iter == 0) {
            // assign points to closest centroid
            for (folan = 0; folan < N; folan++) {
                smallest = DBL_MAX;
                second_smallest = DBL_MAX;
                for (filan = 0; filan < K; filan++) {
                    // TODO: for now I am still filling up the ss anyway so ill use but later ill add a data_squared arr
                    tmp = centroid_squares[filan] + data_arr_ss[folan][0];
                    // the dot product
                    for (ashghal = 0; ashghal < D; ashghal++) {
                        tmp -= (2 * data_arr[folan][ashghal] * centroids[filan][ashghal]);
                    }
                    if(folan == 0 && filan == 0) std::cout << "dist tmp " << tmp << std::endl;
                    if(tmp < 0.0) tmp = 0.0;
                    distances[folan][filan] = sqrt(tmp);
                    if (distances[folan][filan] < smallest) {
                        labels[folan] = filan;
                        smallest = distances[folan][filan];
                    }
                }
                // find second smallest
                for (filan = 0; filan < K; filan++) {
                    if (distances[folan][filan] < second_smallest) {
                        if (distances[folan][filan] == smallest) continue;
                        second_smallest = distances[folan][filan];
                    }
                }

                if(folan == 0) std::cout << "smallest " << smallest << " sec_smallest " << second_smallest << std::endl;

                hamerly_upper_bounds[folan] = smallest;
                hamerly_lower_bounds[folan] = second_smallest;
            }
            std::cout << "calculated distances" << std::endl;
            // bad implementation for now
            // TODO: make it more efficient
            // if first iter, fill out the hamer_ub_lb
            // for(folan = 0; folan < N; folan++){
            //     smallest = DBL_MAX; second_smallest = DBL_MAX;
            //     // find smallest
            //     for(filan = 0; filan < K; filan++){
            //         if(distances[folan][filan] < smallest) smallest = distances[folan][filan];
            //     }
            //     // find second smallest
            //     for(filan = 0; filan < K; filan++){
            //         if(distances[folan][filan] < second_smallest){
            //             if(distances[folan][filan] == smallest) continue;
            //             second_smallest = distances[folan][filan];
            //         } 
            //     }
            //     hamerly_upper_bounds[folan] = smallest;
            //     hamerly_lower_bounds[folan] = second_smallest;
            // }
            std::cout << "filled out hamerly ub and lb for first time" << std::endl;

        } else {
            int hamerly_count = 0;
            for (folan = 0; folan < N; folan++) {
                // hamerly check
                hamerly_bound = ((0.5 * closest_centroid_distance[labels[folan]]) > hamerly_lower_bounds[folan]) ? (
                        0.5 * closest_centroid_distance[labels[folan]]) : hamerly_lower_bounds[folan];
                // hamerly_bound = max(0.5 * closest_centroid_distance[labels[folan]], hamerly_lower_bounds[folan]);
                if (hamerly_bound < hamerly_upper_bounds[folan]) {
                    // calculate correct distance to labels[folan] and update ham_upper
                    tmp = centroid_squares[labels[folan]] + data_arr_ss[folan][0];
                    for (ashghal = 0; ashghal < D; ashghal++) {
                            tmp -= (2 * data_arr[folan][ashghal] * centroids[labels[folan]][ashghal]);
                    }
                    if(tmp < 0.0) tmp = 0.0;
                    distances[folan][labels[folan]] = sqrt(tmp);
                    hamerly_upper_bounds[folan] = distances[folan][labels[folan]]; 
                    // this should tighten the bounds
                    // check again
                    if (hamerly_bound < hamerly_upper_bounds[folan]){
                        for (filan = 0; filan < K; filan++) {

                            tmp = centroid_squares[filan] + data_arr_ss[folan][0];
                            // the dot product
                            for (ashghal = 0; ashghal < D; ashghal++) {
                                tmp -= (2 * data_arr[folan][ashghal] * centroids[filan][ashghal]);
                            }
                            if(tmp < 0.0) tmp = 0.0;
                            distances[folan][filan] = sqrt(tmp);
                            if (distances[folan][filan] < distances[folan][labels[folan]]) {
                                // keep the second smallest
                                hamerly_lower_bounds[folan] = distances[folan][labels[folan]];
                                labels[folan] = filan;
                                // we can just do this once after the for, ghamerly does this too
                                // hamerly_upper_bounds[folan] = distances[folan][labels[folan]];
                            } else if (hamerly_lower_bounds[folan] > distances[folan][filan]) {
                                hamerly_lower_bounds[folan] = distances[folan][filan];
                            }
                        }
                        hamerly_upper_bounds[folan] = distances[folan][labels[folan]];
                    }
                    else{
                        // otherwise we skip this distance calculation and the label remains the same
                        hamerly_count += (K - 1);
                    }
                } else {
                    // otherwise we skip this distance calculation and the label remains the same
                    hamerly_count += K;
                }
            }

            std::cout << "hamerly pruned " << hamerly_count << std::endl;
        }

        // I do this in the loops that calculate the distances now, no need to do it here:)

        // for(folan = 0; folan < N; folan++){
        //     for(filan = 0; filan < K; filan++){
        //         if(distances[folan][filan] < distances[folan][labels[folan]]) {
        //             labels[folan] = filan;
        //         }
        //     }
        // }

        std::cout << "set labels" << std::endl;

        // FATEMEH DEBUG
        std::cout << "data_arr_ss[0][0] " << data_arr_ss[0][0] << std::endl;
        std::cout << "hamerly_ub[0] " << hamerly_upper_bounds[0] << " hamerly_lb[0] " << hamerly_lower_bounds[0] << std::endl;
        // FATEMEH DEBUG


        // calc new centroids
        // make copy of old centroids
        // apparently this also does not work, copies pointer somehow I think...
        // memcpy(old_centroids, centroids, sizeof(double) * K * D);
        for (folan = 0; folan < K; folan++) {
            for (filan = 0; filan < D; filan++) {
                old_centroids[folan][filan] = centroids[folan][filan];
            }
        }
        std::cout << "copied centroids to old centroids" << std::endl;
        // set centroids to 0
        memset(cluster_counts, 0, sizeof(int) * K);
        std::cout << "set cluster counts to 0" << std::endl;
        for (folan = 0; folan < K; folan++) {
            // just testing
            // cluster_counts[folan] = 0;
            for (filan = 0; filan < D; filan++) {
                centroids[folan][filan] = 0.0;
            }
        }
        // This doesn't work on doubles
        // memset(centroids, 0, sizeof(double) * K * D);
        std::cout << "after all the memcpys" << std::endl;

        for (folan = 0; folan < N; folan++) {
            cluster_counts[labels[folan]]++;
            for (filan = 0; filan < D; filan++) {
                centroids[labels[folan]][filan] += data_arr[folan][filan];
            }
        }
        for (folan = 0; folan < K; folan++) {
            // to deal with empty clusters
            // if the cluster is empty, keep the old centroid
            if(cluster_counts[folan] > 0){
                for (filan = 0; filan < D; filan++) {
                    centroids[folan][filan] /= cluster_counts[folan];
                }
            } else{
                for (filan = 0; filan < D; filan++) {
                    centroids[folan][filan] = old_centroids[folan][filan];
                }
            }
        }
        std::cout << "calculated new centroids" << std::endl;
        // just to check
        int sanity_check = 0;
        std::cout << "cluster counts..." << std::endl;
        for (folan = 0; folan < K; folan++) {
            sanity_check += cluster_counts[folan];
            std::cout << cluster_counts[folan] << " ";
        }
        std::cout << sanity_check << std::endl;

        // calculating the movement of new to old cluster centers
        furthest_moving_centroid = 0;
        second_furthest_moving_centroid = 1;
        if(centroid_movement[second_furthest_moving_centroid] > centroid_movement[furthest_moving_centroid]){
            std::swap(furthest_moving_centroid, second_furthest_moving_centroid);
        }
        for (folan = 0; folan < K; folan++) {
            tmp = 0.0;
            for (filan = 0; filan < D; filan++) {
                tmp += ((centroids[folan][filan] - old_centroids[folan][filan]) *
                        (centroids[folan][filan] - old_centroids[folan][filan]));
            }
            if(tmp < 0.0) tmp = 0.0;
            centroid_movement[folan] = sqrt(tmp);
            if (centroid_movement[folan] > centroid_movement[furthest_moving_centroid]){
                second_furthest_moving_centroid = furthest_moving_centroid;
                furthest_moving_centroid = folan;
            }
            else if (centroid_movement[folan] >
                     centroid_movement[second_furthest_moving_centroid])
                second_furthest_moving_centroid = folan;
        }
        std::cout << "calculated centroid movements" << std::endl;
        // std::cout << "centr_movements: " << std::endl;
        // for(folan = 0; folan < K; folan++){
        //     std::cout << centroid_movement[folan] << " ";
        // }
        // std::cout << std::endl;

        // FATEMEH DEBUG
        std::cout << "before updating the bounds, hamerly_ub[0] " << hamerly_upper_bounds[0] << " hamerly_lb[0] "
             << hamerly_lower_bounds[0] << std::endl;
        // FATEMEH DEBUG

        // update upper and lower hamerly bounds based on centroid movements
        for (folan = 0; folan < N; folan++) {
            if (folan == 0) std::cout << "moving ub " << centroid_movement[labels[folan]] << std::endl;
            hamerly_upper_bounds[folan] += centroid_movement[labels[folan]];
            if (labels[folan] == furthest_moving_centroid) {
                if (folan == 0) std::cout << "moving lb " << centroid_movement[second_furthest_moving_centroid] << std::endl;
                hamerly_lower_bounds[folan] -= centroid_movement[second_furthest_moving_centroid];
            } else {
                if (folan == 0) std::cout << "moving lb " << centroid_movement[furthest_moving_centroid] << std::endl;
                hamerly_lower_bounds[folan] -= centroid_movement[furthest_moving_centroid];
            }
        }


        // FATEMEH DEBUG
        std::cout << "after updating the bounds, hamerly_ub[0] " << hamerly_upper_bounds[0] << " hamerly_lb[0] "
             << hamerly_lower_bounds[0] << std::endl;
        // FATEMEH DEBUG


        // check convergence
        // TODO: gonna do it in labels assignment, changed my mind will do it here, DONE
        // has_converged = true;
        // for (folan = 0; folan < K; folan++) {
        //     for (filan = 0; filan < D; filan++) {
        //         if (old_centroids[folan][filan] != centroids[folan][filan]) {
        //             has_converged = false;
        //             break;
        //         }
        //     }
        // }
        // learnt from ghamerly
        has_converged = (0.0 == centroid_movement[furthest_moving_centroid]);
        std::cout << "checked convergence" << std::endl;

        // end if converged
        if (has_converged) break;
    }
}


// elkan
// copied from v5
// NOT TESTED!!!!!
// TODO: also buggy, still waiting for Mo
// DONE: just had to change them to the sqrt thing
// the numbers don't exactly match but that's because of fp
// I changed all the floats to doubles, and that helps a bit, but still a bit different
// I DIDN'T CHECK THE HAMERLY UB CHANGES IN DETAIL, MAYBE LATER TODO
void kmeans_v6() {
    std::cout << "in kmeans_v6 ..." << std::endl;
    int folan = 0;
    int filan = 0;
    int ashghal = 0;// for loop usage
    int furthest_moving_centroid, second_furthest_moving_centroid;
    double smallest, second_smallest;
    double tmp, hamerly_bound;
    bool r;
    // set initial centroids

    for(folan = 0; folan < K; folan++){
        for(filan = 0; filan < D; filan++){
            centroids[folan][filan] = data_arr[folan][filan];
        }
    }

    // memcpy(centroids, data_arr, sizeof(double) * K * D);
    std::cout << "copied init centroids" << std::endl;


    bool has_converged = false;
    int cluster_counts[K];
    // V2 DIFF start
    double centroid_squares[K];
    // V2 DIFF end

    // loop over max_iter
    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
        std::cout << "iteration " << iter << "..." << std::endl;

        // calculate the square sum of centroids
        for (folan = 0; folan < K; folan++) {
            centroid_squares[folan] = 0;
            for (filan = 0; filan < D; filan++) {
                centroid_squares[folan] += (centroids[folan][filan] * centroids[folan][filan]);
            }
        }

        // calculate closest centroid to each centroid
        // TODO: in matris motegharene, ziadi dari hesab mikoni, lazem nis
        // TODO: mitoonim filan=folan shoroo konim k nesfe matris ro por konim
        // vali mitarsam badan estefade azash bad baseh
        for (folan = 0; folan < K; folan++) {
            smallest = DBL_MAX;
            for (filan = folan; filan < K; filan++) {
                tmp = centroid_squares[filan] + centroid_squares[folan];
                for (ashghal = 0; ashghal < D; ashghal++) {
                    tmp -= (2 * centroids[folan][ashghal] * centroids[filan][ashghal]);
                }
                if(tmp < 0.0) tmp = 0.0;
                centroid_to_centroid_distances[folan][filan] = sqrt(tmp);
                // THEY'RE THE SAME
                centroid_to_centroid_distances[filan][folan] = centroid_to_centroid_distances[folan][filan];
                if (centroid_to_centroid_distances[folan][filan] < smallest)
                    smallest = centroid_to_centroid_distances[folan][filan];
            }
            closest_centroid_distance[folan] = smallest;
        }
        std::cout << "found closest distance to each centroid, and the centr-centr distances" << std::endl;


        if (iter == 0) {
            // assign points to closest centroid
            for (folan = 0; folan < N; folan++) {
                smallest = DBL_MAX;
                second_smallest = DBL_MAX;
                for (filan = 0; filan < K; filan++) {
                    tmp = centroid_squares[filan] + data_arr_ss[folan][0];
                    // the dot product
                    for (ashghal = 0; ashghal < D; ashghal++) {
                        tmp -= (2 * data_arr[folan][ashghal] * centroids[filan][ashghal]);
                    }
                    if(tmp < 0.0) tmp = 0.0;
                    distances[folan][filan] = sqrt(tmp);

                    elkan_lower_bounds[folan][filan] = distances[folan][filan];
                    if (distances[folan][filan] < smallest) {
                        labels[folan] = filan;
                        smallest = distances[folan][filan];
                    }
                }
                // find second smallest
                for (filan = 0; filan < K; filan++) {
                    if (distances[folan][filan] < second_smallest) {
                        if (distances[folan][filan] == smallest) continue;
                        second_smallest = distances[folan][filan];
                    }
                }
                hamerly_upper_bounds[folan] = smallest;
                // hamerly_lower_bounds[folan] = second_smallest;
            }
            std::cout << "calculated distances" << std::endl;
            // bad implementation for now
            // TODO: make it more efficient
            // if first iter, fill out the hamer_ub_lb
            // for(folan = 0; folan < N; folan++){
            //     smallest = DBL_MAX; second_smallest = DBL_MAX;
            //     // find smallest
            //     for(filan = 0; filan < K; filan++){
            //         if(distances[folan][filan] < smallest) smallest = distances[folan][filan];
            //     }
            //     // find second smallest
            //     for(filan = 0; filan < K; filan++){
            //         if(distances[folan][filan] < second_smallest){
            //             if(distances[folan][filan] == smallest) continue;
            //             second_smallest = distances[folan][filan];
            //         } 
            //     }
            //     hamerly_upper_bounds[folan] = smallest;
            //     hamerly_lower_bounds[folan] = second_smallest;
            // }
            std::cout << "filled out hamerly ub and elkan lb for first time" << std::endl;

        } else {
            for (folan = 0; folan < N; folan++) {
                r = true;
                // elkan lemma 1
                if ((0.5 * closest_centroid_distance[labels[folan]]) < hamerly_upper_bounds[folan]) {
                    for (filan = 0; filan < K; filan++) {
                        if (filan == labels[folan]) continue;
                        if (hamerly_upper_bounds[folan] > elkan_lower_bounds[folan][filan]
                            && hamerly_upper_bounds[folan] >
                               (0.5 * centroid_to_centroid_distances[filan][labels[folan]])) {

                            // elkan lemma 3a
                            if (r) {
                                // update the upper bound
                                // FATEMEH: I'll update the distances too just in case
                                tmp = centroid_squares[labels[folan]] + data_arr_ss[folan][0];
                                for (ashghal = 0; ashghal < D; ashghal++) {
                                    tmp -= (2 * data_arr[folan][ashghal] * centroids[labels[folan]][ashghal]);
                                }
                                if(tmp < 0.0) tmp = 0.0;
                                distances[folan][labels[folan]] = sqrt(tmp);
                                // I'll do it once after the for
                                // but I'm changing the elkan_lb too, so i'll do it here too
                                hamerly_upper_bounds[folan] = distances[folan][labels[folan]];
                                elkan_lower_bounds[folan][labels[folan]] = hamerly_upper_bounds[folan];
                                r = false;
                                if (folan == 0) {
                                    std::cout
                                            << "ALSO CHANGED THE HAMERLY_UB to same thing, WHICH I HAD THOUGHT WAS HARMLESS, BUT MAYBE I SHOULDN'T "
                                            << std::endl;
                                    std::cout << "changing elkan_lb[0][" << labels[folan] << "] to "
                                         << hamerly_upper_bounds[folan] << std::endl;
                                }
                            }

                            tmp = centroid_squares[filan] + data_arr_ss[folan][0];
                            for (ashghal = 0; ashghal < D; ashghal++) {
                                tmp -= (2 * data_arr[folan][ashghal] * centroids[filan][ashghal]);
                            }
                            if(tmp < 0.0) tmp = 0.0;
                            distances[folan][filan] = sqrt(tmp);
                            elkan_lower_bounds[folan][filan] = distances[folan][filan];
                            if (folan == 0) {
                                std::cout << "changing elkan_lb[0][" << filan << "] to " << distances[folan][filan] << std::endl;
                            }

                            if (distances[folan][filan] < distances[folan][labels[folan]]) {
                                // keep the second smallest
                                // hamerly_lower_bounds[folan] = distances[folan][labels[folan]];
                                labels[folan] = filan;
                                // i am doing this under duress...
                                hamerly_upper_bounds[folan] = distances[folan][labels[folan]];
                                if (folan == 0) {
                                    std::cout << "changing hamerly_ub[0] to " << distances[folan][labels[folan]] << std::endl;
                                }
                            }
                            // else if(hamerly_lower_bounds[folan] > distances[folan][filan]){
                            //     hamerly_lower_bounds[folan] = distances[folan][filan];
                            // }

                        }
                    }
                    // TODO: this is from the hamerly imp, im not sure if i should keep it or not, for now let's cmnt it
                    // DONE: I remembered why, because I don't want to assign it
                    // but the cpp code does it too, so for now ill do it mult times till i debug properly
                    // hamerly_upper_bounds[folan] = distances[folan][labels[folan]];
                }
                // otherwise we skip this distance calculation and the label remains the same
            }
        }

        // I do this in the loops that calculate the distances now, no need to do it here:)

        // for(folan = 0; folan < N; folan++){
        //     for(filan = 0; filan < K; filan++){
        //         if(distances[folan][filan] < distances[folan][labels[folan]]) {
        //             labels[folan] = filan;
        //         }
        //     }
        // }

        std::cout << "set labels" << std::endl;


        // calc new centroids
        // make copy of old centroids
        // apparently this also does not work, copies pointer somehow I think...
        // memcpy(old_centroids, centroids, sizeof(double) * K * D);
        for (folan = 0; folan < K; folan++) {
            for (filan = 0; filan < D; filan++) {
                old_centroids[folan][filan] = centroids[folan][filan];
            }
        }
        std::cout << "copied centroids to old centroids" << std::endl;
        // set centroids to 0
        memset(cluster_counts, 0, sizeof(int) * K);
        std::cout << "set cluster counts to 0" << std::endl;
        for (folan = 0; folan < K; folan++) {
            // just testing
            // cluster_counts[folan] = 0;
            for (filan = 0; filan < D; filan++) {
                centroids[folan][filan] = 0.0;
            }
        }
        // This doesn't work on doubles
        // memset(centroids, 0, sizeof(double) * K * D);
        std::cout << "after all the memcpys" << std::endl;

        for (folan = 0; folan < N; folan++) {
            cluster_counts[labels[folan]]++;
            for (filan = 0; filan < D; filan++) {
                centroids[labels[folan]][filan] += data_arr[folan][filan];
            }
        }
        for (folan = 0; folan < K; folan++) {
            // to deal with empty clusters
            // if the cluster is empty, keep the old centroid
            if(cluster_counts[folan] > 0){
                for (filan = 0; filan < D; filan++) {
                    centroids[folan][filan] /= cluster_counts[folan];
                }
            } else{
                for (filan = 0; filan < D; filan++) {
                    centroids[folan][filan] = old_centroids[folan][filan];
                }
            }
        }
        std::cout << "calculated new centroids" << std::endl;


        std::cout << "centroid 0 ************** " << std::endl;
        for (filan = 120; filan < 140; filan++) {
            std::cout << centroids[0][filan] << " ";
        }
        tmp = 0.0;
        for (filan = 0; filan < D; filan++) tmp += (centroids[0][filan] * centroids[0][filan]);
        std::cout << std::endl;

        std::cout << "sum_ss_centroids[0] " << tmp << std::endl;

        // just to check
        int sanity_check = 0;
        std::cout << "cluster counts..." << std::endl;
        for (folan = 0; folan < K; folan++) {
            sanity_check += cluster_counts[folan];
            std::cout << cluster_counts[folan] << " ";
        }
        std::cout << sanity_check << std::endl;

        // calculating the movement of new to old cluster centers
        furthest_moving_centroid = 0;
        second_furthest_moving_centroid = 1;
        if(centroid_movement[second_furthest_moving_centroid] > centroid_movement[furthest_moving_centroid]){
            std::swap(furthest_moving_centroid, second_furthest_moving_centroid);
        }
        for (folan = 0; folan < K; folan++) {
            tmp = 0.0;
            for (filan = 0; filan < D; filan++) {
                tmp += ((centroids[folan][filan] - old_centroids[folan][filan]) *
                        (centroids[folan][filan] - old_centroids[folan][filan]));
            }
            if(tmp < 0.0) tmp = 0.0;
            centroid_movement[folan] = sqrt(tmp);
            if (centroid_movement[folan] > centroid_movement[furthest_moving_centroid]){
                second_furthest_moving_centroid = furthest_moving_centroid;
                furthest_moving_centroid = folan;
            }
            else if (centroid_movement[folan] >
                     centroid_movement[second_furthest_moving_centroid])
                second_furthest_moving_centroid = folan;
        }
        std::cout << "calculated centroid movements" << std::endl;

        std::cout << "before updating bounds ----------" << std::endl;
        std::cout << "elkan lb[0][0] " << elkan_lower_bounds[0][0] << " hamerly ub[0] " << hamerly_upper_bounds[0] << std::endl;

        // update upper and lower elkan bounds based on centroid movements
        for (folan = 0; folan < N; folan++) {
            hamerly_upper_bounds[folan] += centroid_movement[labels[folan]];
            for (filan = 0; filan < K; filan++) {
                elkan_lower_bounds[folan][filan] -= centroid_movement[filan];
            }
        }

        std::cout << "after updating bounds ----------" << std::endl;
        std::cout << "elkan lb[0][0] " << elkan_lower_bounds[0][0] << " hamerly ub[0] " << hamerly_upper_bounds[0] << std::endl;



        // check convergence
        // TODO: gonna do it in labels assignment, changed my mind will do it here, DONE
        // has_converged = true;
        // for (folan = 0; folan < K; folan++) {
        //     for (filan = 0; filan < D; filan++) {
        //         if (old_centroids[folan][filan] != centroids[folan][filan]) {
        //             has_converged = false;
        //             break;
        //         }
        //     }
        // }
        has_converged = (0.0 == centroid_movement[furthest_moving_centroid]);
        std::cout << "checked convergence" << std::endl;

        // end if converged
        if (has_converged) break;
    }
}


// elkan + new_lemma
// copied from v6
// NOT TESTED!!!!! this is terrible idea since v6 was buggy already:(
// TODO: also buggy, still waiting for Mo
// DONE: changed all the sqrts, the counts of clusters are not right, but im not it might just be the fp error
void kmeans_v7() {
    std::cout << "in kmeans_v7 ..." << std::endl;
    int folan = 0;
    int filan = 0;
    int ashghal = 0;// for loop usage
    int furthest_moving_centroid, second_furthest_moving_centroid;
    double smallest, second_smallest;
    double tmp, hamerly_bound;
    bool r;
    // set initial centroids

    for(folan = 0; folan < K; folan++){
        for(filan = 0; filan < D; filan++){
            centroids[folan][filan] = data_arr[folan][filan];
        }
    }

    // memcpy(centroids, data_arr, sizeof(double) * K * D);
    std::cout << "copied init centroids" << std::endl;


    bool has_converged = false;
    int cluster_counts[K];
    // V2 DIFF start
    double centroid_squares[K];
    // V2 DIFF end

    // loop over max_iter
    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
        std::cout << "iteration " << iter << "..." << std::endl;

        // calculate the square sum of centroids
        for (folan = 0; folan < K; folan++) {
            centroid_squares[folan] = 0;
            for (filan = 0; filan < D; filan++) {
                centroid_squares[folan] += (centroids[folan][filan] * centroids[folan][filan]);
            }
        }

        // calculate closest centroid to each centroid
        // TODO: in matris motegharene, ziadi dari hesab mikoni, lazem nis
        // TODO: mitoonim filan=folan shoroo konim k nesfe matris ro por konim
        // vali mitarsam badan estefade azash bad baseh
        for (folan = 0; folan < K; folan++) {
            smallest = DBL_MAX;
            for (filan = folan; filan < K; filan++) {
                tmp = centroid_squares[filan] + centroid_squares[folan];
                for (ashghal = 0; ashghal < D; ashghal++) {
                    tmp -= (2 * centroids[folan][ashghal] * centroids[filan][ashghal]);
                }
                if(tmp < 0.0) tmp = 0.0;
                centroid_to_centroid_distances[folan][filan] = sqrt(tmp);
                // THEY'RE THE SAME
                centroid_to_centroid_distances[filan][folan] = centroid_to_centroid_distances[folan][filan];
                if (centroid_to_centroid_distances[folan][filan] < smallest)
                    smallest = centroid_to_centroid_distances[folan][filan];
            }
            closest_centroid_distance[folan] = smallest;
        }
        std::cout << "found closest distance to each centroid, and the centr-centr distances" << std::endl;


        if (iter == 0) {
            // assign points to closest centroid
            for (folan = 0; folan < N; folan++) {
                smallest = DBL_MAX;
                second_smallest = DBL_MAX;
                for (filan = 0; filan < K; filan++) {
                    tmp = centroid_squares[filan] + data_arr_ss[folan][0];
                    // the dot product
                    for (ashghal = 0; ashghal < D; ashghal++) {
                        tmp -= (2 * data_arr[folan][ashghal] * centroids[filan][ashghal]);
                    }
                    if(tmp < 0.0) tmp = 0.0;
                    distances[folan][filan] = sqrt(tmp);
                    elkan_lower_bounds[folan][filan] = distances[folan][filan];
                    if (distances[folan][filan] < smallest) {
                        labels[folan] = filan;
                        smallest = distances[folan][filan];
                    }
                }
                // find second smallest
                for (filan = 0; filan < K; filan++) {
                    if (distances[folan][filan] < second_smallest) {
                        if (distances[folan][filan] == smallest) continue;
                        second_smallest = distances[folan][filan];
                    }
                }
                hamerly_upper_bounds[folan] = smallest;
                // hamerly_lower_bounds[folan] = second_smallest;
            }
            std::cout << "calculated distances" << std::endl;
            // bad implementation for now
            // TODO: make it more efficient
            // if first iter, fill out the hamer_ub_lb
            // for(folan = 0; folan < N; folan++){
            //     smallest = DBL_MAX; second_smallest = DBL_MAX;
            //     // find smallest
            //     for(filan = 0; filan < K; filan++){
            //         if(distances[folan][filan] < smallest) smallest = distances[folan][filan];
            //     }
            //     // find second smallest
            //     for(filan = 0; filan < K; filan++){
            //         if(distances[folan][filan] < second_smallest){
            //             if(distances[folan][filan] == smallest) continue;
            //             second_smallest = distances[folan][filan];
            //         } 
            //     }
            //     hamerly_upper_bounds[folan] = smallest;
            //     hamerly_lower_bounds[folan] = second_smallest;
            // }
            std::cout << "filled out hamerly ub and elkan lb for first time" << std::endl;

        } else {
            for (folan = 0; folan < N; folan++) {
                r = true;
                // elkan lemma 1
                if ((0.5 * closest_centroid_distance[labels[folan]]) < hamerly_upper_bounds[folan]) {
                    for (filan = 0; filan < K; filan++) {
                        if (filan == labels[folan]) continue;
                        // new lemma
                        if (elkan_lower_bounds[folan][filan] >
                            hamerly_upper_bounds[folan] + centroid_movement[labels[folan]] + centroid_movement[filan]) {
                            // TODO: i think this should be a break instead of a continue, but im waiting for Mo's response
                            // DONE: I was wrong, it is just continue, cause it just means that filan will not be the correct cluster for this guy
                            continue;
                        }
                        if (hamerly_upper_bounds[folan] > elkan_lower_bounds[folan][filan]
                            && hamerly_upper_bounds[folan] >
                               (0.5 * centroid_to_centroid_distances[filan][labels[folan]])) {

                            // elkan lemma 3a
                            if (r) {
                                // update the upper bound
                                // FATEMEH: I'll update the distances too just in case
                                tmp = centroid_squares[labels[folan]] + data_arr_ss[folan][0];
                                for (ashghal = 0; ashghal < D; ashghal++) {
                                    tmp -= (2 * data_arr[folan][ashghal] * centroids[labels[folan]][ashghal]);
                                }
                                if(tmp < 0.0) tmp = 0.0;
                                distances[folan][labels[folan]] = sqrt(tmp);
                                // I'll do it once after the for
                                // but I'm changing the elkan_lb too, so i'll do it here too
                                hamerly_upper_bounds[folan] = distances[folan][labels[folan]];
                                elkan_lower_bounds[folan][labels[folan]] = hamerly_upper_bounds[folan];
                                r = false;
                            }

                            tmp = centroid_squares[filan] + data_arr_ss[folan][0];
                            for (ashghal = 0; ashghal < D; ashghal++) {
                                tmp -= (2 * data_arr[folan][ashghal] * centroids[filan][ashghal]);
                            }
                            if(tmp < 0.0) tmp = 0.0;
                            distances[folan][filan] = sqrt(tmp);
                            elkan_lower_bounds[folan][filan] = distances[folan][filan];

                            if (distances[folan][filan] < distances[folan][labels[folan]]) {
                                // keep the second smallest
                                // hamerly_lower_bounds[folan] = distances[folan][labels[folan]];
                                labels[folan] = filan;
                                // i am doing this under duress...
                                hamerly_upper_bounds[folan] = distances[folan][labels[folan]];
                            }
                            // else if(hamerly_lower_bounds[folan] > distances[folan][filan]){
                            //     hamerly_lower_bounds[folan] = distances[folan][filan];
                            // }

                        }
                    }
                    // TODO: this is from the hamerly imp, im not sure if i should keep it or not, for now let's cmnt it
                    // DONE: I remembered why, because I don't want to assign it
                    // but the cpp code does it too, so for now ill do it mult times till i debug properly
                    // hamerly_upper_bounds[folan] = distances[folan][labels[folan]];
                }
                // otherwise we skip this distance calculation and the label remains the same
            }
        }

        // I do this in the loops that calculate the distances now, no need to do it here:)

        // for(folan = 0; folan < N; folan++){
        //     for(filan = 0; filan < K; filan++){
        //         if(distances[folan][filan] < distances[folan][labels[folan]]) {
        //             labels[folan] = filan;
        //         }
        //     }
        // }

        std::cout << "set labels" << std::endl;


        // calc new centroids
        // make copy of old centroids
        // apparently this also does not work, copies pointer somehow I think...
        // memcpy(old_centroids, centroids, sizeof(double) * K * D);
        for (folan = 0; folan < K; folan++) {
            for (filan = 0; filan < D; filan++) {
                old_centroids[folan][filan] = centroids[folan][filan];
            }
        }
        std::cout << "copied centroids to old centroids" << std::endl;
        // set centroids to 0
        memset(cluster_counts, 0, sizeof(int) * K);
        std::cout << "set cluster counts to 0" << std::endl;
        for (folan = 0; folan < K; folan++) {
            // just testing
            // cluster_counts[folan] = 0;
            for (filan = 0; filan < D; filan++) {
                centroids[folan][filan] = 0.0;
            }
        }
        // This doesn't work on doubles
        // memset(centroids, 0, sizeof(double) * K * D);
        std::cout << "after all the memcpys" << std::endl;

        for (folan = 0; folan < N; folan++) {
            cluster_counts[labels[folan]]++;
            for (filan = 0; filan < D; filan++) {
                centroids[labels[folan]][filan] += data_arr[folan][filan];
            }
        }
        for (folan = 0; folan < K; folan++) {
            // to deal with empty clusters
            // if the cluster is empty, keep the old centroid
            if(cluster_counts[folan] > 0){
                for (filan = 0; filan < D; filan++) {
                    centroids[folan][filan] /= cluster_counts[folan];
                }
            } else{
                for (filan = 0; filan < D; filan++) {
                    centroids[folan][filan] = old_centroids[folan][filan];
                }
            }
        }
        std::cout << "calculated new centroids" << std::endl;
        // just to check
        int sanity_check = 0;
        std::cout << "cluster counts..." << std::endl;
        for (folan = 0; folan < K; folan++) {
            sanity_check += cluster_counts[folan];
            std::cout << cluster_counts[folan] << " ";
        }
        std::cout << sanity_check << std::endl;

        // calculating the movement of new to old cluster centers
        furthest_moving_centroid = 0;
        second_furthest_moving_centroid = 1;
        if(centroid_movement[second_furthest_moving_centroid] > centroid_movement[furthest_moving_centroid]){
            std::swap(furthest_moving_centroid, second_furthest_moving_centroid);
        }
        for (folan = 0; folan < K; folan++) {
            tmp = 0.0;
            for (filan = 0; filan < D; filan++) {
                tmp += ((centroids[folan][filan] - old_centroids[folan][filan]) *
                        (centroids[folan][filan] - old_centroids[folan][filan]));
            }
            if(tmp < 0.0) tmp = 0.0;
            centroid_movement[folan] = sqrt(tmp);
            if (centroid_movement[folan] > centroid_movement[furthest_moving_centroid]){
                second_furthest_moving_centroid = furthest_moving_centroid;
                furthest_moving_centroid = folan;
            }
            else if (centroid_movement[folan] >
                     centroid_movement[second_furthest_moving_centroid])
                second_furthest_moving_centroid = folan;
        }
        std::cout << "calculated centroid movements" << std::endl;

        // update upper and lower elkan bounds based on centroid movements
        for (folan = 0; folan < N; folan++) {
            hamerly_upper_bounds[folan] += centroid_movement[labels[folan]];
            for (filan = 0; filan < K; filan++) {
                elkan_lower_bounds[folan][filan] -= centroid_movement[filan];
            }
        }


        // check convergence
        // TODO: gonna do it in labels assignment, changed my mind will do it here, DONE
        // has_converged = true;
        // for (folan = 0; folan < K; folan++) {
        //     for (filan = 0; filan < D; filan++) {
        //         if (old_centroids[folan][filan] != centroids[folan][filan]) {
        //             has_converged = false;
        //             break;
        //         }
        //     }
        // }
        has_converged = (0.0 == centroid_movement[furthest_moving_centroid]);
        std::cout << "checked convergence" << std::endl;

        // end if converged
        if (has_converged) break;
    }
}


// elkan + new_lemma + hamerly
// copied from v7
// TODO: not even compiled, but i think I have filled out the correct parts
// DONE: compiled, runs, the cluster counts till iter 2 look ok
void kmeans_v8() {
    std::cout << "in kmeans_v8 ..." << std::endl;
    int folan = 0;
    int filan = 0;
    int ashghal = 0;// for loop usage
    int furthest_moving_centroid, second_furthest_moving_centroid;
    double smallest, second_smallest;
    double tmp, hamerly_bound;
    bool r;
    // set initial centroids
    for(folan = 0; folan < K; folan++){
        for(filan = 0; filan < D; filan++){
            centroids[folan][filan] = data_arr[folan][filan];
        }
    }
    // memcpy(centroids, data_arr, sizeof(double) * K * D);
    std::cout << "copied init centroids" << std::endl;


    bool has_converged = false;
    int cluster_counts[K];
    double centroid_squares[K];

    // loop over max_iter
    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
        std::cout << "iteration " << iter << "..." << std::endl;

        // calculate the square sum of centroids
        for (folan = 0; folan < K; folan++) {
            centroid_squares[folan] = 0;
            for (filan = 0; filan < D; filan++) {
                centroid_squares[folan] += (centroids[folan][filan] * centroids[folan][filan]);
            }
        }

        // calculate closest centroid to each centroid
        // TODO: in matris motegharene, ziadi dari hesab mikoni, lazem nis
        // TODO: mitoonim filan=folan shoroo konim k nesfe matris ro por konim
        // vali mitarsam badan estefade azash bad baseh
        for (folan = 0; folan < K; folan++) {
            smallest = DBL_MAX;
            for (filan = folan; filan < K; filan++) {
                tmp = centroid_squares[filan] + centroid_squares[folan];
                for (ashghal = 0; ashghal < D; ashghal++) {
                    tmp -= (2 * centroids[folan][ashghal] * centroids[filan][ashghal]);
                }
                if(tmp < 0.0) tmp = 0.0;
                centroid_to_centroid_distances[folan][filan] = sqrt(tmp);
                // THEY'RE THE SAME
                centroid_to_centroid_distances[filan][folan] = centroid_to_centroid_distances[folan][filan];
                if (centroid_to_centroid_distances[folan][filan] < smallest)
                    smallest = centroid_to_centroid_distances[folan][filan];
            }
            closest_centroid_distance[folan] = smallest;
        }
        std::cout << "found closest distance to each centroid, and the centr-centr distances" << std::endl;


        if (iter == 0) {
            // assign points to closest centroid
            for (folan = 0; folan < N; folan++) {
                smallest = DBL_MAX;
                second_smallest = DBL_MAX;
                for (filan = 0; filan < K; filan++) {
                    tmp = centroid_squares[filan] + data_arr_ss[folan][0];
                    // the dot product
                    for (ashghal = 0; ashghal < D; ashghal++) {
                        tmp -= (2 * data_arr[folan][ashghal] * centroids[filan][ashghal]);
                    }
                    if(tmp < 0.0) tmp = 0.0;
                    distances[folan][filan] = sqrt(tmp);
                    elkan_lower_bounds[folan][filan] = distances[folan][filan];
                    if (distances[folan][filan] < smallest) {
                        labels[folan] = filan;
                        smallest = distances[folan][filan];
                    }
                }
                // find second smallest
                for (filan = 0; filan < K; filan++) {
                    if (distances[folan][filan] < second_smallest) {
                        if (distances[folan][filan] == smallest) continue;
                        second_smallest = distances[folan][filan];
                    }
                }
                hamerly_upper_bounds[folan] = smallest;
                hamerly_lower_bounds[folan] = second_smallest;
            }
            std::cout << "calculated distances" << std::endl;

            std::cout << "filled out hamerly ub and lb and elkan lb for first time" << std::endl;

        } else {
            for (folan = 0; folan < N; folan++) {
                r = true;
                hamerly_bound = ((0.5 * closest_centroid_distance[labels[folan]]) > hamerly_lower_bounds[folan]) ? (
                        0.5 * closest_centroid_distance[labels[folan]]) : hamerly_lower_bounds[folan];
                // elkan lemma 1 and hamerly
                if (hamerly_bound < hamerly_upper_bounds[folan]) {
                    // if((0.5 * closest_centroid_distance[labels[folan]]) < hamerly_upper_bounds[folan]){
                    for (filan = 0; filan < K; filan++) {
                        if (filan == labels[folan]) continue;
                        // new lemma
                        if (elkan_lower_bounds[folan][filan] >
                            hamerly_upper_bounds[folan] + centroid_movement[labels[folan]] + centroid_movement[filan]) {
                            // TODO: i think this should be a break instead of a continue, but im waiting for Mo's response
                            // DONE: I was wrong, it is just continue, cause it just means that filan will not be the correct cluster for this guy
                            continue;
                        }
                        if (hamerly_upper_bounds[folan] > elkan_lower_bounds[folan][filan]
                            && hamerly_upper_bounds[folan] >
                               (0.5 * centroid_to_centroid_distances[filan][labels[folan]])) {

                            // elkan lemma 3a
                            if (r) {
                                // update the upper bound
                                // FATEMEH: I'll update the distances too just in case
                                tmp = centroid_squares[labels[folan]] + data_arr_ss[folan][0];
                                for (ashghal = 0; ashghal < D; ashghal++) {
                                    tmp -= (2 * data_arr[folan][ashghal] * centroids[labels[folan]][ashghal]);
                                }
                                if(tmp < 0.0) tmp = 0.0;
                                distances[folan][labels[folan]] = sqrt(tmp);
                                // I'll do it once after the for
                                // but I'm changing the elkan_lb too, so i'll do it here too
                                hamerly_upper_bounds[folan] = distances[folan][labels[folan]];
                                elkan_lower_bounds[folan][labels[folan]] = hamerly_upper_bounds[folan];
                                r = false;
                            }

                            tmp = centroid_squares[filan] + data_arr_ss[folan][0];
                            for (ashghal = 0; ashghal < D; ashghal++) {
                                tmp -= (2 * data_arr[folan][ashghal] * centroids[filan][ashghal]);
                            }
                            if(tmp < 0.0) tmp = 0.0;
                            distances[folan][filan] = sqrt(tmp);
                            elkan_lower_bounds[folan][filan] = distances[folan][filan];

                            if (distances[folan][filan] < distances[folan][labels[folan]]) {
                                // keep the second smallest
                                hamerly_lower_bounds[folan] = distances[folan][labels[folan]];
                                labels[folan] = filan;
                                // i am doing this under duress...
                                hamerly_upper_bounds[folan] = distances[folan][labels[folan]];
                            } else if (hamerly_lower_bounds[folan] > distances[folan][filan]) {
                                hamerly_lower_bounds[folan] = distances[folan][filan];
                            }

                        }
                    }
                    // TODO: this is from the hamerly imp, im not sure if i should keep it or not, for now let's cmnt it
                    // DONE: I remembered why, because I don't want to assign it
                    // but the cpp code does it too, so for now ill do it mult times till i debug properly
                    // hamerly_upper_bounds[folan] = distances[folan][labels[folan]];
                }
                // otherwise we skip this distance calculation and the label remains the same
            }
        }

        // I do this in the loops that calculate the distances now, no need to do it here:)

        // for(folan = 0; folan < N; folan++){
        //     for(filan = 0; filan < K; filan++){
        //         if(distances[folan][filan] < distances[folan][labels[folan]]) {
        //             labels[folan] = filan;
        //         }
        //     }
        // }

        std::cout << "set labels" << std::endl;


        // calc new centroids
        // make copy of old centroids
        // apparently this also does not work, copies pointer somehow I think...
        // memcpy(old_centroids, centroids, sizeof(double) * K * D);
        for (folan = 0; folan < K; folan++) {
            for (filan = 0; filan < D; filan++) {
                old_centroids[folan][filan] = centroids[folan][filan];
            }
        }
        std::cout << "copied centroids to old centroids" << std::endl;
        // set centroids to 0
        memset(cluster_counts, 0, sizeof(int) * K);
        std::cout << "set cluster counts to 0" << std::endl;
        for (folan = 0; folan < K; folan++) {
            // just testing
            // cluster_counts[folan] = 0;
            for (filan = 0; filan < D; filan++) {
                centroids[folan][filan] = 0.0;
            }
        }
        // This doesn't work on doubles
        // memset(centroids, 0, sizeof(double) * K * D);
        std::cout << "after all the memcpys" << std::endl;

        for (folan = 0; folan < N; folan++) {
            cluster_counts[labels[folan]]++;
            for (filan = 0; filan < D; filan++) {
                centroids[labels[folan]][filan] += data_arr[folan][filan];
            }
        }
        for (folan = 0; folan < K; folan++) {
            // to deal with empty clusters
            // if the cluster is empty, keep the old centroid
            if(cluster_counts[folan] > 0){
                for (filan = 0; filan < D; filan++) {
                    centroids[folan][filan] /= cluster_counts[folan];
                }
            } else{
                for (filan = 0; filan < D; filan++) {
                    centroids[folan][filan] = old_centroids[folan][filan];
                }
            }
        }
        std::cout << "calculated new centroids" << std::endl;
        // just to check
        int sanity_check = 0;
        std::cout << "cluster counts..." << std::endl;
        for (folan = 0; folan < K; folan++) {
            sanity_check += cluster_counts[folan];
            std::cout << cluster_counts[folan] << " ";
        }
        std::cout << sanity_check << std::endl;

        // calculating the movement of new to old cluster centers
        furthest_moving_centroid = 0;
        second_furthest_moving_centroid = 1;
        if(centroid_movement[second_furthest_moving_centroid] > centroid_movement[furthest_moving_centroid]){
            std::swap(furthest_moving_centroid, second_furthest_moving_centroid);
        }
        for (folan = 0; folan < K; folan++) {
            tmp = 0.0;
            for (filan = 0; filan < D; filan++) {
                tmp += ((centroids[folan][filan] - old_centroids[folan][filan]) *
                        (centroids[folan][filan] - old_centroids[folan][filan]));
            }
            if(tmp < 0.0) tmp = 0.0;
            centroid_movement[folan] = sqrt(tmp);
            if (centroid_movement[folan] > centroid_movement[furthest_moving_centroid]){
                second_furthest_moving_centroid = furthest_moving_centroid;
                furthest_moving_centroid = folan;
            }
            else if (centroid_movement[folan] >
                     centroid_movement[second_furthest_moving_centroid])
                second_furthest_moving_centroid = folan;
        }
        std::cout << "calculated centroid movements" << std::endl;

        // update upper and lower elkan bounds based on centroid movements
        for (folan = 0; folan < N; folan++) {
            hamerly_upper_bounds[folan] += centroid_movement[labels[folan]];
            if (labels[folan] == furthest_moving_centroid) {
                hamerly_lower_bounds[folan] -= centroid_movement[second_furthest_moving_centroid];
            } else {
                hamerly_lower_bounds[folan] -= centroid_movement[furthest_moving_centroid];
            }

            for (filan = 0; filan < K; filan++) {
                elkan_lower_bounds[folan][filan] -= centroid_movement[filan];
            }
        }


        // check convergence
        // TODO: gonna do it in labels assignment, changed my mind will do it here, DONE
        // has_converged = true;
        // for (folan = 0; folan < K; folan++) {
        //     for (filan = 0; filan < D; filan++) {
        //         if (old_centroids[folan][filan] != centroids[folan][filan]) {
        //             has_converged = false;
        //             break;
        //         }
        //     }
        // }
        has_converged = (0.0 == centroid_movement[furthest_moving_centroid]);
        std::cout << "checked convergence" << std::endl;

        // end if converged
        if (has_converged) break;
    }
}


// hamerly + step-wise
// copied from v5
// TODO
void kmeans_v9() {
    std::cout << "in kmeans_v9 ..." << std::endl;
    int folan = 0;
    int filan = 0;
    int ashghal = 0;// for loop usage
    int furthest_moving_centroid, second_furthest_moving_centroid;
    double smallest, second_smallest;
    double tmp, hamerly_bound;
    int hamerly_count = 0;
    // set initial centroids

    for(folan = 0; folan < K; folan++){
        for(filan = 0; filan < D; filan++){
            centroids[folan][filan] = data_arr[folan][filan];
        }
    }

    // memcpy(centroids, data_arr, sizeof(double) * K * D);
    calculate_centroids_square_sums();
    std::cout << "copied init centroids" << std::endl;


    bool has_converged = false;
    int cluster_counts[K];
    // V2 DIFF start
    double centroid_squares[K];
    // V2 DIFF end

    // loop over max_iter
    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
        std::cout << "iteration " << iter << "..." << std::endl;

        // calculate the square sum of centroids
        for (folan = 0; folan < K; folan++) {
            centroid_squares[folan] = 0;
            for (filan = 0; filan < D; filan++) {
                centroid_squares[folan] += (centroids[folan][filan] * centroids[folan][filan]);
            }
        }

        // calculate closest centroid to each centroid
        for (folan = 0; folan < K; folan++) {
            smallest = DBL_MAX;
            for (filan = 0; filan < K; filan++) {
                tmp = centroid_squares[filan] + centroid_squares[folan];
                for (ashghal = 0; ashghal < D; ashghal++) {
                    tmp -= (2 * centroids[folan][ashghal] * centroids[filan][ashghal]);
                }
                if (tmp < smallest) smallest = tmp;
            }
            if(smallest < 0.0) smallest = 0.0;
            closest_centroid_distance[folan] = sqrt(smallest);
        }
        std::cout << "found closest distance to each centroid" << std::endl;


        if (iter == 0) {
            std::cout << "in if iter == 0" << std::endl;
            // assign points to closest centroid
            // V9
            // for(folan = 0; folan < N; folan++){
            //     smallest = DBL_MAX; second_smallest = DBL_MAX;
            //     for(filan =0; filan < K; filan++){
            //         // TODO: for now I am still filling up the ss anyway so ill use but later ill add a data_squared arr
            //         tmp = centroid_squares[filan] + data_arr_ss[folan][0];
            //         // the dot product
            //         for(ashghal = 0; ashghal < D; ashghal++){
            //             tmp -= (2 * data_arr[folan][ashghal] * centroids[filan][ashghal]);
            //         }
            //         distances[folan][filan] =  sqrt(tmp);
            //         if(distances[folan][filan] < smallest){ 
            //             labels[folan] = filan;
            //             smallest = distances[folan][filan];
            //         }
            //     }
            //     // find second smallest
            //     for(filan = 0; filan < K; filan++){
            //         if(distances[folan][filan] < second_smallest){
            //             if(distances[folan][filan] == smallest) continue;
            //             second_smallest = distances[folan][filan];
            //         } 
            //     }

            //     hamerly_upper_bounds[folan] = smallest;
            //     hamerly_lower_bounds[folan] = second_smallest;
            // }
            // the is_cand is all one from cause it's the first time from main
            calculate_labels_with_sqrt_hamerly_integrated();
            // V9
            std::cout << "calculated distances" << std::endl;
            // bad implementation for now
            // TODO: make it more efficient
            // if first iter, fill out the hamer_ub_lb
            // for(folan = 0; folan < N; folan++){
            //     smallest = DBL_MAX; second_smallest = DBL_MAX;
            //     // find smallest
            //     for(filan = 0; filan < K; filan++){
            //         if(distances[folan][filan] < smallest) smallest = distances[folan][filan];
            //     }
            //     // find second smallest
            //     for(filan = 0; filan < K; filan++){
            //         if(distances[folan][filan] < second_smallest){
            //             if(distances[folan][filan] == smallest) continue;
            //             second_smallest = distances[folan][filan];
            //         } 
            //     }
            //     hamerly_upper_bounds[folan] = smallest;
            //     hamerly_lower_bounds[folan] = second_smallest;
            // }
            std::cout << "filled out hamerly ub and lb for first time" << std::endl;

        } else {
            std::cout << "iter was not 0" << std::endl;
            // V9
            // fill out mask
            // set all to 1 first
            for (folan = 0; folan < N; folan++) {
                for (filan = 0; filan < K; filan++) {
                    is_candidate[folan][filan] = true;
                }
            }
            // now let's do pruning based on TI
            for (folan = 0; folan < N; folan++) {
                // hamerly check
                hamerly_bound = ((0.5 * closest_centroid_distance[labels[folan]]) > hamerly_lower_bounds[folan]) ? (
                        0.5 * closest_centroid_distance[labels[folan]]) : hamerly_lower_bounds[folan];
                if (hamerly_bound >= hamerly_upper_bounds[folan]) {
                    hamerly_count += K;
                    // this one will not move, so none of them are candidates
                    for (filan = 0; filan < K; filan++) {
                        is_candidate[folan][filan] = false;
                    }
                    // TODO: we may have to leave the labels[folan] = true, but i don't think so
                    // I will leave it as false for now and change it if necessary
                }
            }
            std::cout << "hamerly pruned " << hamerly_count << std::endl;
            calculate_labels_with_sqrt_hamerly_integrated();
            // we do all the next lines in the calc_labels so no need
            // for(folan = 0; folan < N; folan++){    
            //     if(hamerly_bound < hamerly_upper_bounds[folan]){
            //         for(filan = 0; filan < K; filan++){
            //                 tmp = centroid_squares[filan] + data_arr_ss[folan][0];
            //                 // the dot product
            //                 for(ashghal = 0; ashghal < D; ashghal++){
            //                     tmp -= (2 * data_arr[folan][ashghal] * centroids[filan][ashghal]);
            //                 }
            //                 distances[folan][filan] = sqrt(tmp);
            //                 if(distances[folan][filan] < distances[folan][labels[folan]]){ 
            //                     // keep the second smallest
            //                     hamerly_lower_bounds[folan] = distances[folan][labels[folan]];
            //                     labels[folan] = filan;
            //                     hamerly_upper_bounds[folan] = distances[folan][labels[folan]];
            //                 }
            //                 else if(hamerly_lower_bounds[folan] > distances[folan][filan]){
            //                     hamerly_lower_bounds[folan] = distances[folan][filan];
            //                 }
            //         }
            //         hamerly_upper_bounds[folan] = distances[folan][labels[folan]];
            //     }
            // V9
            // otherwise we skip this distance calculation and the label remains the same

        }

        // I do this in the loops that calculate the distances now, no need to do it here:)

        // for(folan = 0; folan < N; folan++){
        //     for(filan = 0; filan < K; filan++){
        //         if(distances[folan][filan] < distances[folan][labels[folan]]) {
        //             labels[folan] = filan;
        //         }
        //     }
        // }

        std::cout << "set labels" << std::endl;




        // calc new centroids
        // make copy of old centroids
        // apparently this also does not work, copies pointer somehow I think...
        // memcpy(old_centroids, centroids, sizeof(double) * K * D);
        for (folan = 0; folan < K; folan++) {
            for (filan = 0; filan < D; filan++) {
                old_centroids[folan][filan] = centroids[folan][filan];
            }
        }
        std::cout << "copied centroids to old centroids" << std::endl;
        // set centroids to 0
        memset(cluster_counts, 0, sizeof(int) * K);
        std::cout << "set cluster counts to 0" << std::endl;
        for (folan = 0; folan < K; folan++) {
            // just testing
            // cluster_counts[folan] = 0;
            for (filan = 0; filan < D; filan++) {
                centroids[folan][filan] = 0.0;
            }
        }
        // This doesn't work on doubles
        // memset(centroids, 0, sizeof(double) * K * D);
        std::cout << "after all the memcpys" << std::endl;

        for (folan = 0; folan < N; folan++) {
            cluster_counts[labels[folan]]++;
            for (filan = 0; filan < D; filan++) {
                centroids[labels[folan]][filan] += data_arr[folan][filan];
            }
        }
        for (folan = 0; folan < K; folan++) {
            // to deal with empty clusters
            // if the cluster is empty, keep the old centroid
            if(cluster_counts[folan] > 0){
                for (filan = 0; filan < D; filan++) {
                    centroids[folan][filan] /= cluster_counts[folan];
                }
            } else{
                for (filan = 0; filan < D; filan++) {
                    centroids[folan][filan] = old_centroids[folan][filan];
                }
            }
        }
        std::cout << "calculated new centroids" << std::endl;
        calculate_centroids_square_sums();

        // just to check
        int sanity_check = 0;
        std::cout << "cluster counts..." << std::endl;
        for (folan = 0; folan < K; folan++) {
            sanity_check += cluster_counts[folan];
            std::cout << cluster_counts[folan] << " ";
        }
        std::cout << sanity_check << std::endl;

        // calculating the movement of new to old cluster centers
        furthest_moving_centroid = 0;
        second_furthest_moving_centroid = 1;
        if(centroid_movement[second_furthest_moving_centroid] > centroid_movement[furthest_moving_centroid]){
            std::swap(furthest_moving_centroid, second_furthest_moving_centroid);
        }
        for (folan = 0; folan < K; folan++) {
            tmp = 0.0;
            for (filan = 0; filan < D; filan++) {
                tmp += ((centroids[folan][filan] - old_centroids[folan][filan]) *
                        (centroids[folan][filan] - old_centroids[folan][filan]));
            }
            if(tmp < 0.0) tmp = 0.0;
            centroid_movement[folan] = sqrt(tmp);
            if (centroid_movement[folan] > centroid_movement[furthest_moving_centroid]){
                second_furthest_moving_centroid = furthest_moving_centroid;
                furthest_moving_centroid = folan;
            }
            else if (centroid_movement[folan] >
                     centroid_movement[second_furthest_moving_centroid])
                second_furthest_moving_centroid = folan;
        }
        std::cout << "calculated centroid movements" << std::endl;

        // update upper and lower hamerly bounds based on centroid movements
        for (folan = 0; folan < N; folan++) {
            if (folan == 0) std::cout << "moving ub " << centroid_movement[labels[folan]] << std::endl;
            hamerly_upper_bounds[folan] += centroid_movement[labels[folan]];
            if (labels[folan] == furthest_moving_centroid) {
                if (folan == 0) std::cout << "moving lb " << centroid_movement[second_furthest_moving_centroid] << std::endl;
                hamerly_lower_bounds[folan] -= centroid_movement[second_furthest_moving_centroid];
            } else {
                if (folan == 0) std::cout << "moving lb " << centroid_movement[furthest_moving_centroid] << std::endl;
                hamerly_lower_bounds[folan] -= centroid_movement[furthest_moving_centroid];
            }
        }


        // FATEMEH DEBUG
        std::cout << "after updating the bounds, hamerly_ub[0] " << hamerly_upper_bounds[0] << " hamerly_lb[0] "
             << hamerly_lower_bounds[0] << std::endl;
        // FATEMEH DEBUG


        // check convergence
        // TODO: gonna do it in labels assignment, changed my mind will do it here, DONE
        // has_converged = true;
        // for (folan = 0; folan < K; folan++) {
        //     for (filan = 0; filan < D; filan++) {
        //         if (old_centroids[folan][filan] != centroids[folan][filan]) {
        //             has_converged = false;
        //             break;
        //         }
        //     }
        // }
        has_converged = (0.0 == centroid_movement[furthest_moving_centroid]);
        std::cout << "checked convergence" << std::endl;

        // end if converged
        if (has_converged) break;
    }

    std::cout << "hamerly pruned " << hamerly_count << std::endl;
}


int main(int argc, char **argv) {
    // allocate all arrays
    labels = (int *) calloc(N, sizeof(int));

    data_arr = (double **) malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        data_arr[i] = (double *) malloc(D * sizeof(double));
    }

    centroids = (double **) malloc(K * sizeof(double *));
    for (int i = 0; i < K; i++) {
        centroids[i] = (double *) malloc(D * sizeof(double));
    }

    old_centroids = (double **) malloc(K * sizeof(double *));
    for (int i = 0; i < K; i++) {
        old_centroids[i] = (double *) malloc(D * sizeof(double));
    }

    distances = (double **) malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        distances[i] = (double *) malloc(K * sizeof(double));
    }


    std::cout << "after allocation" << std::endl;

    // read data from file
    // put into "data_arr" array
    // TODO: DONE
    std::string data_file_name = argv[1];
    std::string label_file_name = argv[2];

    std::cout << "data file name: " << data_file_name << std::endl;
    std::cout << "label file name: " << label_file_name << std::endl;


    std::ifstream data_file(data_file_name.c_str());
    int data_count = std::count(std::istreambuf_iterator<char>(data_file),
                                std::istreambuf_iterator<char>(), '\n');
    if (data_count < N) {
        std::cout << "NOT ENOUGH DATA!!\n";
        exit(3);
    }

    std::cout << "count of data in file: " << data_count << " N: " << N << std::endl;

    data_file.clear();
    data_file.seekg(0, std::ios::beg);

    double tmp;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < D; j++) {
            // data_file >> tmp;
            // data_arr[i][j] = round(tmp / 100);
            data_file >> data_arr[i][j];
        }
    }
    data_file.close();

    std::cout << "read data..." << std::endl;

    // STEPWISE EXTRAS
    is_candidate = (bool **) malloc(N * sizeof(bool *));
    for (int i = 0; i < N; i++) {
        is_candidate[i] = (bool *) malloc(K * sizeof(bool));
        // TODO: check the memset trick instead of this
        for (int j = 0; j < K; j++) {
            is_candidate[i][j] = true;
        }
    }
    upper_bounds = (double **) malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        upper_bounds[i] = (double *) malloc(K * sizeof(double));
    }
    lower_bounds = (double **) malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        lower_bounds[i] = (double *) malloc(K * sizeof(double));
    }
    incremental_dots = (double **) malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        incremental_dots[i] = (double *) malloc(K * sizeof(double));
        // can't memset doubles
        for (int j = 0; j < K; j++) {
            incremental_dots[i][j] = 0.0;
        }
    }

    data_arr_ss = (double **) malloc(N * sizeof(double *));

    for (int i = 0; i < N; i++) {
        (data_arr_ss)[i] = (double *) malloc((int(log2(int(sqrt(D))) + 3) * sizeof(double)));
    }

    centroids_ss = (double **) malloc(K * sizeof(double *));
    for (int i = 0; i < K; i++) {
        (centroids_ss)[i] = (double *) malloc((int(log2(int(sqrt(D))) + 3) * sizeof(double)));
    }
    calculate_data_square_sums();
    // let's check if these are correct

    // STEPWISE EXTRAS end

    // HAMERLY EXTRAS
    hamerly_lower_bounds = (double *) calloc(N, sizeof(double));
    hamerly_upper_bounds = (double *) calloc(N, sizeof(double));
    closest_centroid_distance = (double *) calloc(K, sizeof(double));
    centroid_movement = (double *) calloc(K, sizeof(double));
    // HAMERLY EXTRAS end

    // ELKAN EXTRAS
    elkan_lower_bounds = (double **) malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        elkan_lower_bounds[i] = (double *) malloc(K * sizeof(double));
    }
    centroid_to_centroid_distances = (double **) malloc(K * sizeof(double *));
    for (int i = 0; i < K; i++) {
        centroid_to_centroid_distances[i] = (double *) malloc(K * sizeof(double));
    }
    // ELKAN EXTRAS end

    // just helper
    assigned = (int *) malloc(N * sizeof(int));
    // just helper

    std::cout << "filled out the ss arr..." << std::endl;


    // set labels to 0 init state
    memset(labels, 0, sizeof(int) * N);

    std::cout << "set labels to 0, calling kmeans..." << std::endl;

    // do the clustering
    kmeans_v4();

    // write labels to somewhere I guess...
    // TODO
    std::ofstream label_file;
    label_file.open(label_file_name.c_str());

    for (int i = 0; i < N; i++) {
        label_file << labels[i] << "\n";
    }
    label_file.close();


    return 0;
}

