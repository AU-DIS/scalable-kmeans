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
//#include <format>
#include <float.h>
#include <cblas.h>
#include "benchmark/benchmark.h"

char **argv_;
int argc_;

    #ifdef DEBUG
    #define DEBUG_TEST 1
    #else
    #define DEBUG_TEST 0
    #endif
    #define DEBUGPRINT(fmt, ...) \
                do { if (DEBUG_TEST) fprintf(stderr, fmt, __VA_ARGS__); } while (0)


    //#ifndef D 
    //#define D 16384// data dimensionality
    //#endif

    /*#ifndef K
    #define K 40 // k variable in kmeans
    #endif
    #ifndef N
    #define N 168 // count of data points
    #endif*/
    #ifndef HYBRID_SWITCH_THRESHOLD
    #define HYBRID_SWITCH_THRESHOLD 1024
    #endif
    #ifndef MAX_ITERATIONS
    #define MAX_ITERATIONS 100 // maximum number of iterations to do before giving up convergence
    #endif

double **data_arr;

class Kmeans_bench {

    // double data_arr[N][D]; // the data itself 
    // int labels[N]; // where the labels of the algorithm will be eventually
    // double centroids[K][D]; // centroid of clusters
    // double distances[N][K]; // distance of data[i] to centroid[j]
    private:
    
    
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

    // MARIGOLD EXTRAS
    int *smallest_ub;
    int **last_level_calculated;
    // MARIGOLD EXTRAS end

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

    public:
    int *labels;
    const int D; // = 16384;
    const int N;
    const int K;
    Kmeans_bench(int D, int N, int K) : D(D), N(N), K(K) {
        //std::cout << D << " " << N << " " << K;
    };

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
        ////std::cout << "in kmeans..." << std::endl;
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
        //std::cout << "copied init centroids" << std::endl;
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
        //std::cout << sizeof(cluster_counts) << "\n";
        //std::cout << K << "\n";

        // loop over max_iter
        for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
            ////std::cout << "iteration " << iter << "..." << std::endl;
            // assign points to closest centroid
            // TODO : DONE
            for (folan = 0; folan < N; folan++) {
                for (filan = 0; filan < K; filan++) {
                    distances[folan][filan] = euclidean_distance(folan, filan);
                }
            }
            ////std::cout << "calculated distances" << std::endl;
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

            ////std::cout << "set labels" << std::endl;


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
            //std::cout << "copied centroids to old centroids" << std::endl;
            // set centroids to 0
            memset(cluster_counts, 0, sizeof(int) * K);
            //std::cout << "set cluster counts to 0" << std::endl;
            for (folan = 0; folan < K; folan++) {
                // just testing
                // cluster_counts[folan] = 0;
                for (filan = 0; filan < D; filan++) {
                    centroids[folan][filan] = 0.0;
                }
            }
            // This doesn't work on doubles
            // memset(centroids, 0, sizeof(double) * K * D);
            //std::cout << "after all the memcpys" << std::endl;

            for (folan = 0; folan < N; folan++) {
                cluster_counts[labels[folan]]++;
                for (filan = 0; filan < D; filan++) {
                    centroids[labels[folan]][filan] += data_arr[folan][filan];
                }
            }

            /*for(int i = 0; i < K; i++) {
                for (int j = 0; j < 10; j++) {
                    std::cout << centroids[i][j] << " ";
                }
                std::cout << "\n";
            }
            std::cout << std::endl;*/


            for (folan = 0; folan < K; folan++) {
                //std::cout << cluster_counts[folan] << std::endl;
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

            

            //std::cout << "calculated new centroids" << std::endl;
            // just to check
            /*int sanity_check = 0;
            std::cout << "cluster counts..." << std::endl;
            for (folan = 0; folan < K; folan++) {
                sanity_check += cluster_counts[folan];
                std::cout << cluster_counts[folan] << " ";
            }
            std::cout << sanity_check << std::endl;*/
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
            //std::cout << "checked convergence" << std::endl;

            // end if converged
            if (has_converged) {
                if (D >= 8*8) {
                    int sanity_check = 0;
                    for (folan = 0; folan < K; folan++) {
                        sanity_check += cluster_counts[folan];
                        std::cout << cluster_counts[folan] << " ";
                    }
                    std::cout << sanity_check << std::endl;
                    std::cout << "Final iter " << iter << std::endl;
                }
                break;
            }
        }
    }


    // in this version I will calculatet the dists a bit differently
    // ||x-c|| ^2 = ||x||^2 + ||c||^2 - 2X.CT
    // the ||x||^2 part doesn't matter in the comparisons I will do, so I'll skip them
    // the ||c||^2 I will calculate only once and use all the time
    // the x.cT I will first do simply with a for loop
    // TESTED
    void kmeans_v2() {
        //std::cout << "in kmeans_v2 ..." << std::endl;
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
        //std::cout << "copied init centroids" << std::endl;
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
            //std::cout << "iteration " << iter << "..." << std::endl;
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
            //std::cout << "calculated distances" << std::endl;

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

            //std::cout << "set labels" << std::endl;


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
            //std::cout << "copied centroids to old centroids" << std::endl;
            // set centroids to 0
            memset(cluster_counts, 0, sizeof(int) * K);
            //std::cout << "set cluster counts to 0" << std::endl;
            for (folan = 0; folan < K; folan++) {
                // just testing
                // cluster_counts[folan] = 0;
                for (filan = 0; filan < D; filan++) {
                    centroids[folan][filan] = 0.0;
                }
            }
            // This doesn't work on doubles
            // memset(centroids, 0, sizeof(double) * K * D);
            //std::cout << "after all the memcpys" << std::endl;

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
            //std::cout << "calculated new centroids" << std::endl;
            // just to check
            int sanity_check = 0;
            //std::cout << "cluster counts..." << std::endl;
            for (folan = 0; folan < K; folan++) {
                sanity_check += cluster_counts[folan];
                //std::cout << cluster_counts[folan] << " ";
            }
            //std::cout << sanity_check << std::endl;
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
            //std::cout << "checked convergence" << std::endl;

            // end if converged
            if (has_converged) {
                DEBUGPRINT("final iter: %d", iter);
                break;
            }
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

        //std::cout << "after memset labels" << std::endl;


        int level = 0;
        // int level = int(log2(int(sqrt(D))) +  1);
        // I forgot that you can't do this to doubles
        // will do it in main
        // memset(incremental_dots, 0, sizeof(double) * K * N);

        //int *smallest_ub = (int *) malloc(N * sizeof(int));

        //std::cout << "after init labels and inc_dots and smallest_ub" << std::endl;

        while (level < int(log2(int(sqrt(D))) + 2)) {
            //std::cout << "in level while loop, level = " << level << std::endl;
            calculate_distances_till_level(level);
            //std::cout << "after calculate dist till level" << std::endl;



            // find the smallest upperbound per point
            memset(smallest_ub, 0, sizeof(int) * N);
            for (folan = 0; folan < N; folan++) {
                for (filan = 0; filan < K; filan++) {
                    if(!is_candidate[folan][filan]) continue;
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
        //         //std::cout << "TERROR! DISASTER! WE WERE DECIEVED:(" << std::endl;
        //     }
        // }
        // end of sanity check
        for (folan = 0; folan < N; folan++) labels[folan] = -1;

        //std::cout << "after memset labels" << std::endl;


        int level = 0;
        // int level = int(log2(int(sqrt(D))) +  1);
        // I forgot that you can't do this to doubles
        // will do it in main
        // memset(incremental_dots, 0, sizeof(double) * K * N);

        //int *smallest_ub = (int *) malloc(N * sizeof(int));

        //std::cout << "after init labels and inc_dots and smallest_ub" << std::endl;

        while (level < int(log2(int(sqrt(D))) + 2)) {
            //std::cout << "in level while loop, level = " << level << std::endl;
            calculate_sqrt_distances_till_level(level);
            //std::cout << "after calculate dist till level" << std::endl;



            // find the smallest upperbound per point
            memset(smallest_ub, 0, sizeof(int) * N);
            for (folan = 0; folan < N; folan++) {
                for (filan = 0; filan < K; filan++) {
                    if(!is_candidate[folan][filan]) continue;
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
        //         //std::cout << "TERROR! DISASTER! WE WERE DECIEVED:(" << std::endl;
        //     }
        // }
        // end of sanity check
        for (folan = 0; folan < N; folan++) labels[folan] = -1;

        //std::cout << "after memset labels" << std::endl;


        int level = 0;
        // int level = int(log2(int(sqrt(D))) +  1);
        // I forgot that you can't do this to doubles
        // will do it in main
        // memset(incremental_dots, 0, sizeof(double) * K * N);

        //int *smallest_ub = (int *) malloc(N * sizeof(int));

        //std::cout << "after init labels and inc_dots and smallest_ub" << std::endl;

        while (level < int(log2(int(sqrt(D))) + 2)) {
            //std::cout << "in level while loop, level = " << level << std::endl;
            calculate_sqrt_distances_till_level(level);
            //std::cout << "after calculate dist till level" << std::endl;



            // find the smallest upperbound per point
            memset(smallest_ub, 0, sizeof(int) * N);
            for (folan = 0; folan < N; folan++) {
                for (filan = 0; filan < K; filan++) {
                    if(!is_candidate[folan][filan]) continue;
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
    // fills up the hamerly ub and lbs
    void calculate_labels_with_sqrt_hamerly_integrated() {
    //std::cout << "in calc_labels_ham_integrated" << std::endl;
    int folan, filan, ashghal, alaki;
    bool candidates_exist;
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
    // memset(assigned, 0, sizeof(int) * N);
    // V9
    //std::cout << "after init assigned and inc_dots and smallest_ub" << std::endl;
    while (level < int(log2(int(sqrt(D))) + 2)) {
        //std::cout << "in level while loop, level = " << level << std::endl;
        calculate_sqrt_distances_till_level_with_assigned(level);
        //std::cout << "after calculate dist till level" << std::endl;
        // find the smallest upperbound per point
        memset(smallest_ub, 0, sizeof(int) * N);
        for (folan = 0; folan < N; folan++) {
            for (filan = 0; filan < K; filan++) {
                if(!is_candidate[folan][filan]) continue;
                if (upper_bounds[folan][filan] < upper_bounds[folan][smallest_ub[folan]]) smallest_ub[folan] = filan;
            }
        }

        for (folan = 0; folan < N; folan++) {
            // V9
            // if(labels[folan] > 0) continue;
            if (assigned[folan] > 0) continue;
            // V9
            candidates_exist = false;
            if (level == int(log2(int(sqrt(D))) + 1)) {
                labels[folan] = smallest_ub[folan];
                // V9
                hamerly_upper_bounds[folan] = upper_bounds[folan][smallest_ub[folan]];
                assigned[folan] = 1;
                // v9

                // but we still need to update the ham_lb
                for (filan = 0; filan < K; filan++) {
                    if (filan == smallest_ub[folan]) continue;
                    if (!is_candidate[folan][filan]) continue;
                    if(!candidates_exist){
                        candidates_exist = true;
                        // fake_smallest_lb = DBL_MAX;
                        hamerly_lower_bounds[folan] = DBL_MAX;
                    }

                    // TODO: so im changing the definition of ham_lb
                    // let ham_lb be the smallest lb betweent the last batch of candidates
                    // if there are other candidates at THIS point, then the ham_lb should just be fake_smallest_lb here
                    // but what if there are no other candidates?
                    // then we would have set fake_smallest_lb to DBL_MAX and it stays that way...
                    // I fixed it with the exists_candidate
                    // it's not the label == smallest_ub
                    // so it's the smallest lb to get pruned
                    if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                        hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                    }

                }
            } else {
                for (filan = 0; filan < K; filan++) {
                    if (filan == smallest_ub[folan]) continue;
                    if (!is_candidate[folan][filan]) continue;
                    if (lower_bounds[folan][filan] >= upper_bounds[folan][smallest_ub[folan]]) {
                        is_candidate[folan][filan] = false;
                        if(!candidates_exist){
                            candidates_exist = true;
                            // fake_smallest_lb = DBL_MAX;
                            hamerly_lower_bounds[folan] = DBL_MAX;
                        }
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


    // copied from calculate_labels_with_sqrt_hamerly_integrated
// just calls the sqrt_dist function
// fills up the hamerly and elkan ub and lbs
void calculate_labels_with_sqrt_hamerly_elkan_integrated() {
    std::cout << "in calc_labels_ham_elk_integrated" << std::endl;
    int folan, filan, ashghal, alaki;
    bool candidates_exist;
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
    // memset(assigned, 0, sizeof(int) * N);
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
                if(!is_candidate[folan][filan]) continue;
                elkan_lower_bounds[folan][filan] = lower_bounds[folan][filan];
                if (upper_bounds[folan][filan] < upper_bounds[folan][smallest_ub[folan]]) smallest_ub[folan] = filan;
            }
        }
        for (folan = 0; folan < N; folan++) {
            // V9
            // if(labels[folan] > 0) continue;
            if (assigned[folan] > 0) continue;
            // V9
            candidates_exist = false;
            if (level == int(log2(int(sqrt(D))) + 1)) {
                labels[folan] = smallest_ub[folan];
                // V9
                hamerly_upper_bounds[folan] = upper_bounds[folan][smallest_ub[folan]];
                assigned[folan] = 1;
                // v9
                // but we still need to update the ham_lb
                for (filan = 0; filan < K; filan++) {
                    if (filan == smallest_ub[folan]) continue;
                    if (!is_candidate[folan][filan]) continue;
                    if(!candidates_exist){
                        candidates_exist = true;
                        // fake_smallest_lb = DBL_MAX;
                        hamerly_lower_bounds[folan] = DBL_MAX;
                    }
                    // TODO: so im changing the definition of ham_lb
                    // let ham_lb be the smallest lb betweent the last batch of candidates
                    // if there are other candidates at THIS point, then the ham_lb should just be fake_smallest_lb here
                    // but what if there are no other candidates?
                    // then we would have set fake_smallest_lb to DBL_MAX and it stays that way...
                    // I fixed it with the exists_candidate
                    // it's not the label == smallest_ub
                    // so it's the smallest lb to get pruned
                    if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                        hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                    }
                }
            } else {
                for (filan = 0; filan < K; filan++) {
                    if (filan == smallest_ub[folan]) continue;
                    if (!is_candidate[folan][filan]) continue;
                    if (lower_bounds[folan][filan] >= upper_bounds[folan][smallest_ub[folan]]) {
                        is_candidate[folan][filan] = false;
                        if(!candidates_exist){
                            candidates_exist = true;
                            // fake_smallest_lb = DBL_MAX;
                            hamerly_lower_bounds[folan] = DBL_MAX;
                        }
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
        //std::cout << "in kmeans_v4 ..." << std::endl;
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
        //memcpy(centroids, data_arr, sizeof(double) * K * D);


        calculate_centroids_square_sums();
        //std::cout << "copied init centroids, and ss-ed them" << std::endl;


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
            //std::cout << "iteration " << iter << "..." << std::endl;

            // assign points to closest centroid
            // TODO 
            // V4 DIFF
            // reset the is_candidates
            for (folan = 0; folan < N; folan++) {
                for (filan = 0; filan < K; filan++) {
                    is_candidate[folan][filan] = true;
                }
            }
             calculate_labels();
            // calculate_labels_with_sqrt();
            //calculate_labels_with_sqrt_hamerly_integrated();

            // V4 DIFF end
            //std::cout << "set labels" << std::endl;


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
            //std::cout << "copied centroids to old centroids" << std::endl;
            // set centroids to 0
            memset(cluster_counts, 0, sizeof(int) * K);
            //std::cout << "set cluster counts to 0" << std::endl;
            for (folan = 0; folan < K; folan++) {
                // just testing
                // cluster_counts[folan] = 0;
                for (filan = 0; filan < D; filan++) {
                    centroids[folan][filan] = 0.0;
                }
            }
            // This doesn't work on doubles
            // memset(centroids, 0, sizeof(double) * K * D);
            //std::cout << "after all the memcpys" << std::endl;

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
            //std::cout << "calculated new centroids" << std::endl;
            calculate_centroids_square_sums();
            //std::cout << "and ss-ed them" << std::endl;
            // just to check
            /*int sanity_check = 0;
            //std::cout << "cluster counts..." << std::endl;
            for (folan = 0; folan < K; folan++) {
                sanity_check += cluster_counts[folan];
                //std::cout << cluster_counts[folan] << " ";
            }*/
            //std::cout << sanity_check << std::endl;
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
            //std::cout << "checked convergence" << std::endl;

            // end if converged
            if (has_converged) {
                if (D >= 256*256) {
                    int sanity_check = 0;
                    for (folan = 0; folan < K; folan++) {
                        sanity_check += cluster_counts[folan];
                        std::cout << cluster_counts[folan] << " ";
                    }
                    std::cout << sanity_check << std::endl;
                    std::cout << "Final iter " << iter << std::endl;
                }
                break;
            }
        }
    }


    // hamerly
    // copied from v2
    // TODO: STILL BUGGY, WAITING FOR MOHAMMAD TO FIX HIS CODE SO I CAN CHECK THE UB AND LB
    // UPDATE: I think it's the x^2 and sqrt that I removed then the TE doesn't make sense anymore so I have to add them back 
    // DONE: holy shit, that was it:)) I just added the sqrt and x^2 everywhere and now it's correct:)
    void kmeans_v5() {
        //std::cout << "in kmeans_v5 ..." << std::endl;
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
        //std::cout << "copied init centroids" << std::endl;


        bool has_converged = false;
        int cluster_counts[K];
        // V2 DIFF start
        double centroid_squares[K];
        // V2 DIFF end

        // loop over max_iter
        for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
            //std::cout << "iteration " << iter << "..." << std::endl;

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
            //std::cout << "found closest distance to each centroid" << std::endl;


            if (iter == 0) {
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
                    // COUNT
                    //features_accessed += D;
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

                //if(folan == 0) std::cout << "smallest " << smallest << " sec_smallest " << second_smallest << std::endl;

                hamerly_upper_bounds[folan] = smallest;
                hamerly_lower_bounds[folan] = second_smallest;
                }
                //std::cout << "calculated distances" << std::endl;
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
                //std::cout << "filled out hamerly ub and lb for first time" << std::endl;

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

                //std::cout << "hamerly pruned " << hamerly_count << std::endl;
            }

            // I do this in the loops that calculate the distances now, no need to do it here:)

            // for(folan = 0; folan < N; folan++){
            //     for(filan = 0; filan < K; filan++){
            //         if(distances[folan][filan] < distances[folan][labels[folan]]) {
            //             labels[folan] = filan;
            //         }
            //     }
            // }

            //std::cout << "set labels" << std::endl;

            // FATEMEH DEBUG
            //std::cout << "data_arr_ss[0][0] " << data_arr_ss[0][0] << std::endl;
            //std::cout << "hamerly_ub[0] " << hamerly_upper_bounds[0] << " hamerly_lb[0] " << hamerly_lower_bounds[0] << std::endl;
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
            //std::cout << "copied centroids to old centroids" << std::endl;
            // set centroids to 0
            memset(cluster_counts, 0, sizeof(int) * K);
            //std::cout << "set cluster counts to 0" << std::endl;
            for (folan = 0; folan < K; folan++) {
                // just testing
                // cluster_counts[folan] = 0;
                for (filan = 0; filan < D; filan++) {
                    centroids[folan][filan] = 0.0;
                }
            }
            // This doesn't work on doubles
            // memset(centroids, 0, sizeof(double) * K * D);
            //std::cout << "after all the memcpys" << std::endl;

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
            //std::cout << "calculated new centroids" << std::endl;
            // just to check
            int sanity_check = 0;
            //std::cout << "cluster counts..." << std::endl;
            for (folan = 0; folan < K; folan++) {
                sanity_check += cluster_counts[folan];
                //std::cout << cluster_counts[folan] << " ";
            }
            //std::cout << sanity_check << std::endl;

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
            //std::cout << "calculated centroid movements" << std::endl;
            // std::cout << "centr_movements: " << std::endl;
            // for(folan = 0; folan < K; folan++){
            //     std::cout << centroid_movement[folan] << " ";
            // }
            // std::cout << std::endl;

            // FATEMEH DEBUG
            //std::cout << "before updating the bounds, hamerly_ub[0] " << hamerly_upper_bounds[0] << " hamerly_lb[0] "
            //     << hamerly_lower_bounds[0] << std::endl;
            // FATEMEH DEBUG

            // update upper and lower hamerly bounds based on centroid movements
            for (folan = 0; folan < N; folan++) {
                //if (folan == 0) std::cout << "moving ub " << centroid_movement[labels[folan]] << std::endl;
                hamerly_upper_bounds[folan] += centroid_movement[labels[folan]];
                if (labels[folan] == furthest_moving_centroid) {
                    //if (folan == 0) std::cout << "moving lb " << centroid_movement[second_furthest_moving_centroid] << std::endl;
                    hamerly_lower_bounds[folan] -= centroid_movement[second_furthest_moving_centroid];
                } else {
                    //if (folan == 0) std::cout << "moving lb " << centroid_movement[furthest_moving_centroid] << std::endl;
                    hamerly_lower_bounds[folan] -= centroid_movement[furthest_moving_centroid];
                }
            }


            // FATEMEH DEBUG
            //std::cout << "after updating the bounds, hamerly_ub[0] " << hamerly_upper_bounds[0] << " hamerly_lb[0] "
            //     << hamerly_lower_bounds[0] << std::endl;
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
            //std::cout << "checked convergence" << std::endl;

            if (has_converged) {
                if (D >= 256*256) {
                    int sanity_check = 0;
                    for (folan = 0; folan < K; folan++) {
                        sanity_check += cluster_counts[folan];
                        std::cout << cluster_counts[folan] << " ";
                    }
                    std::cout << sanity_check << std::endl;
                }
                break;
            }
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
        //std::cout << "in kmeans_v6 ..." << std::endl;
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
        //std::cout << "copied init centroids" << std::endl;


        bool has_converged = false;
        int cluster_counts[K];
        // V2 DIFF start
        double centroid_squares[K];
        // V2 DIFF end

        // loop over max_iter
        for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
            //std::cout << "iteration " << iter << "..." << std::endl;

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
            //std::cout << "found closest distance to each centroid, and the centr-centr distances" << std::endl;


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
                        if(elkan_lower_bounds[folan][filan] < second_smallest && filan != labels[folan]){
                            second_smallest = elkan_lower_bounds[folan][filan];
                        }
                    }
                    hamerly_upper_bounds[folan] = smallest;
                    // hamerly_lower_bounds[folan] = second_smallest;
                }
                //std::cout << "calculated distances" << std::endl;
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
                //std::cout << "filled out hamerly ub and elkan lb for first time" << std::endl;

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
                                        //std::cout << "ALSO CHANGED THE HAMERLY_UB to same thing, WHICH I HAD THOUGHT WAS HARMLESS, BUT MAYBE I SHOULDN'T " << std::endl;
                                        //std::cout << "changing elkan_lb[0][" << labels[folan] << "] to "<< hamerly_upper_bounds[folan] << std::endl;
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
                                    //std::cout << "changing elkan_lb[0][" << filan << "] to " << distances[folan][filan] << std::endl;
                                }

                                if (distances[folan][filan] < distances[folan][labels[folan]]) {
                                    // keep the second smallest
                                    // hamerly_lower_bounds[folan] = distances[folan][labels[folan]];
                                    labels[folan] = filan;
                                    // i am doing this under duress...
                                    hamerly_upper_bounds[folan] = distances[folan][labels[folan]];
                                    if (folan == 0) {
                                        //std::cout << "changing hamerly_ub[0] to " << distances[folan][labels[folan]] << std::endl;
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

            //std::cout << "set labels" << std::endl;


            // calc new centroids
            // make copy of old centroids
            // apparently this also does not work, copies pointer somehow I think...
            // memcpy(old_centroids, centroids, sizeof(double) * K * D);
            for (folan = 0; folan < K; folan++) {
                for (filan = 0; filan < D; filan++) {
                    old_centroids[folan][filan] = centroids[folan][filan];
                }
            }
            //std::cout << "copied centroids to old centroids" << std::endl;
            // set centroids to 0
            memset(cluster_counts, 0, sizeof(int) * K);
            //std::cout << "set cluster counts to 0" << std::endl;
            for (folan = 0; folan < K; folan++) {
                // just testing
                // cluster_counts[folan] = 0;
                for (filan = 0; filan < D; filan++) {
                    centroids[folan][filan] = 0.0;
                }
            }
            // This doesn't work on doubles
            // memset(centroids, 0, sizeof(double) * K * D);
            //std::cout << "after all the memcpys" << std::endl;

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
            //std::cout << "calculated new centroids" << std::endl;


            //std::cout << "centroid 0 ************** " << std::endl;
            for (filan = 120; filan < 140; filan++) {
                //std::cout << centroids[0][filan] << " ";
            }
            tmp = 0.0;
            for (filan = 0; filan < D; filan++) tmp += (centroids[0][filan] * centroids[0][filan]);
            //std::cout << std::endl;

            //std::cout << "sum_ss_centroids[0] " << tmp << std::endl;

            // just to check
            int sanity_check = 0;
            //std::cout << "cluster counts..." << std::endl;
            for (folan = 0; folan < K; folan++) {
                sanity_check += cluster_counts[folan];
                //std::cout << cluster_counts[folan] << " ";
            }
            //std::cout << sanity_check << std::endl;

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
            //std::cout << "calculated centroid movements" << std::endl;

            //std::cout << "before updating bounds ----------" << std::endl;
            //std::cout << "elkan lb[0][0] " << elkan_lower_bounds[0][0] << " hamerly ub[0] " << hamerly_upper_bounds[0] << std::endl;

            // update upper and lower elkan bounds based on centroid movements
            for (folan = 0; folan < N; folan++) {
                hamerly_upper_bounds[folan] += centroid_movement[labels[folan]];
                for (filan = 0; filan < K; filan++) {
                    elkan_lower_bounds[folan][filan] -= centroid_movement[filan];
                }
            }

            //std::cout << "after updating bounds ----------" << std::endl;
            //std::cout << "elkan lb[0][0] " << elkan_lower_bounds[0][0] << " hamerly ub[0] " << hamerly_upper_bounds[0] << std::endl;



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
            //std::cout << "checked convergence" << std::endl;

            // end if converged
            if (has_converged) {
                if (D >= 256*256) {
                    int sanity_check = 0;
                    for (folan = 0; folan < K; folan++) {
                        sanity_check += cluster_counts[folan];
                        std::cout << cluster_counts[folan] << " ";
                    }
                    std::cout << sanity_check << std::endl;
                }
                break;
            }
        }
    }


    // elkan + new_lemma
    // copied from v6
    // NOT TESTED!!!!! this is terrible idea since v6 was buggy already:(
    // TODO: also buggy, still waiting for Mo
    // DONE: changed all the sqrts, the counts of clusters are not right, but im not it might just be the fp error
    void kmeans_v7() {
        //std::cout << "in kmeans_v7 ..." << std::endl;
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
        //std::cout << "copied init centroids" << std::endl;


        bool has_converged = false;
        int cluster_counts[K];
        // V2 DIFF start
        double centroid_squares[K];
        // V2 DIFF end

        // loop over max_iter
        for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
            //std::cout << "iteration " << iter << "..." << std::endl;

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
            //std::cout << "found closest distance to each centroid, and the centr-centr distances" << std::endl;


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
                        if(elkan_lower_bounds[folan][filan] < second_smallest && filan != labels[folan]){
                            second_smallest = elkan_lower_bounds[folan][filan];
                        }
                    }
                    hamerly_upper_bounds[folan] = smallest;
                    // hamerly_lower_bounds[folan] = second_smallest;
                }
                //std::cout << "calculated distances" << std::endl;
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
                //std::cout << "filled out hamerly ub and elkan lb for first time" << std::endl;

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

            //std::cout << "set labels" << std::endl;


            // calc new centroids
            // make copy of old centroids
            // apparently this also does not work, copies pointer somehow I think...
            // memcpy(old_centroids, centroids, sizeof(double) * K * D);
            for (folan = 0; folan < K; folan++) {
                for (filan = 0; filan < D; filan++) {
                    old_centroids[folan][filan] = centroids[folan][filan];
                }
            }
            //std::cout << "copied centroids to old centroids" << std::endl;
            // set centroids to 0
            memset(cluster_counts, 0, sizeof(int) * K);
            //std::cout << "set cluster counts to 0" << std::endl;
            for (folan = 0; folan < K; folan++) {
                // just testing
                // cluster_counts[folan] = 0;
                for (filan = 0; filan < D; filan++) {
                    centroids[folan][filan] = 0.0;
                }
            }
            // This doesn't work on doubles
            // memset(centroids, 0, sizeof(double) * K * D);
            //std::cout << "after all the memcpys" << std::endl;

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
            //std::cout << "calculated new centroids" << std::endl;
            // just to check
            int sanity_check = 0;
            //std::cout << "cluster counts..." << std::endl;
            for (folan = 0; folan < K; folan++) {
                sanity_check += cluster_counts[folan];
                //std::cout << cluster_counts[folan] << " ";
            }
            //std::cout << sanity_check << std::endl;

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
            //std::cout << "calculated centroid movements" << std::endl;

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
            //std::cout << "checked convergence" << std::endl;

            // end if converged
            if (has_converged) {
                DEBUGPRINT("final iter: %d", iter);
                break;
            }
        }
    }

    //Hammerly + New Lemma
    void kmeans_v7_5() {
        //std::cout << "in kmeans_v8 ..." << std::endl;
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
        //std::cout << "copied init centroids" << std::endl;


        bool has_converged = false;
        int cluster_counts[K];
        double centroid_squares[K];

        // loop over max_iter
        for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
            //std::cout << "iteration " << iter << "..." << std::endl;

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
            //std::cout << "found closest distance to each centroid, and the centr-centr distances" << std::endl;


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
                        if(elkan_lower_bounds[folan][filan] < second_smallest && filan != labels[folan]){
                            second_smallest = elkan_lower_bounds[folan][filan];
                        }
                    }
                    hamerly_upper_bounds[folan] = smallest;
                    hamerly_lower_bounds[folan] = second_smallest;
                }
                //std::cout << "calculated distances" << std::endl;

                //std::cout << "filled out hamerly ub and lb and elkan lb for first time" << std::endl;

            } else {
                for (folan = 0; folan < N; folan++) {
                    r = true;
                    hamerly_bound = ((0.5 * closest_centroid_distance[labels[folan]]) > hamerly_lower_bounds[folan]) ? (
                            0.5 * closest_centroid_distance[labels[folan]]) : hamerly_lower_bounds[folan];
                    // elkan lemma 1 and hamerly
                    if (hamerly_bound < hamerly_upper_bounds[folan]) {
                        // if((0.5 * closest_centroid_distance[labels[folan]]) < hamerly_upper_bounds[folan]){
                        //SAME TO HERE 1326
                        for (filan = 0; filan < K; filan++) {
                            if (filan == labels[folan]) continue;
                            // new lemma
                            if (elkan_lower_bounds[folan][filan] >
                                hamerly_upper_bounds[folan] + centroid_movement[labels[folan]] + centroid_movement[filan]) {
                                // TODO: i think this should be a break instead of a continue, but im waiting for Mo's response
                                // DONE: I was wrong, it is just continue, cause it just means that filan will not be the correct cluster for this guy
                                continue;
                            }
                            //WOOP STEP
                            if (hamerly_upper_bounds[folan] >
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
                                    //hamerly_upper_bounds[folan] = distances[folan][labels[folan]];
                                    //elkan_lower_bounds[folan][labels[folan]] = hamerly_upper_bounds[folan];
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

                            } //WOOPSTEP
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

            //std::cout << "set labels" << std::endl;


            // calc new centroids
            // make copy of old centroids
            // apparently this also does not work, copies pointer somehow I think...
            // memcpy(old_centroids, centroids, sizeof(double) * K * D);
            for (folan = 0; folan < K; folan++) {
                for (filan = 0; filan < D; filan++) {
                    old_centroids[folan][filan] = centroids[folan][filan];
                }
            }
            //std::cout << "copied centroids to old centroids" << std::endl;
            // set centroids to 0
            memset(cluster_counts, 0, sizeof(int) * K);
            //std::cout << "set cluster counts to 0" << std::endl;
            for (folan = 0; folan < K; folan++) {
                // just testing
                // cluster_counts[folan] = 0;
                for (filan = 0; filan < D; filan++) {
                    centroids[folan][filan] = 0.0;
                }
            }
            // This doesn't work on doubles
            // memset(centroids, 0, sizeof(double) * K * D);
            //std::cout << "after all the memcpys" << std::endl;

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
            //std::cout << "calculated new centroids" << std::endl;
            // just to check
            int sanity_check = 0;
            //std::cout << "cluster counts..." << std::endl;
            for (folan = 0; folan < K; folan++) {
                sanity_check += cluster_counts[folan];
                //std::cout << cluster_counts[folan] << " ";
            }
            //std::cout << sanity_check << std::endl;

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
            //std::cout << "calculated centroid movements" << std::endl;

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
            //std::cout << "checked convergence" << std::endl;

            // end if converged
            if (has_converged) {
                DEBUGPRINT("final iter: %d", iter);
                break;
            }
        }
    }


    // elkan + new_lemma + hamerly
    // copied from v7
    // TODO: not even compiled, but i think I have filled out the correct parts
    // DONE: compiled, runs, the cluster counts till iter 2 look ok
    void kmeans_v8() {
        //std::cout << "in kmeans_v8 ..." << std::endl;
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
        //std::cout << "copied init centroids" << std::endl;


        bool has_converged = false;
        int cluster_counts[K];
        double centroid_squares[K];

        // loop over max_iter
        for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
            //std::cout << "iteration " << iter << "..." << std::endl;

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
            //std::cout << "found closest distance to each centroid, and the centr-centr distances" << std::endl;


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
                        if(elkan_lower_bounds[folan][filan] < second_smallest && filan != labels[folan]){
                            second_smallest = elkan_lower_bounds[folan][filan];
                        }
                    }
                    hamerly_upper_bounds[folan] = smallest;
                    hamerly_lower_bounds[folan] = second_smallest;
                }
                //std::cout << "calculated distances" << std::endl;

                //std::cout << "filled out hamerly ub and lb and elkan lb for first time" << std::endl;

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

            //std::cout << "set labels" << std::endl;


            // calc new centroids
            // make copy of old centroids
            // apparently this also does not work, copies pointer somehow I think...
            // memcpy(old_centroids, centroids, sizeof(double) * K * D);
            for (folan = 0; folan < K; folan++) {
                for (filan = 0; filan < D; filan++) {
                    old_centroids[folan][filan] = centroids[folan][filan];
                }
            }
            //std::cout << "copied centroids to old centroids" << std::endl;
            // set centroids to 0
            memset(cluster_counts, 0, sizeof(int) * K);
            //std::cout << "set cluster counts to 0" << std::endl;
            for (folan = 0; folan < K; folan++) {
                // just testing
                // cluster_counts[folan] = 0;
                for (filan = 0; filan < D; filan++) {
                    centroids[folan][filan] = 0.0;
                }
            }
            // This doesn't work on doubles
            // memset(centroids, 0, sizeof(double) * K * D);
            //std::cout << "after all the memcpys" << std::endl;

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
            //std::cout << "calculated new centroids" << std::endl;
            // just to check
            int sanity_check = 0;
            //std::cout << "cluster counts..." << std::endl;
            for (folan = 0; folan < K; folan++) {
                sanity_check += cluster_counts[folan];
                //std::cout << cluster_counts[folan] << " ";
            }
            //std::cout << sanity_check << std::endl;

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
            //std::cout << "calculated centroid movements" << std::endl;

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
            //std::cout << "checked convergence" << std::endl;

            // end if converged
            if (has_converged) {
                DEBUGPRINT("final iter: %d", iter);
                break;
            }
        }
    }

// elkan + hamerly
// copied from v8
// deleted redundant new lemma
// and the redundant r variable
void kmeans_v85() {
    //std::cout << "in kmeans_v85 ..." << std::endl;
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
    //std::cout << "copied init centroids" << std::endl;


    bool has_converged = false;
    int cluster_counts[K];
    double centroid_squares[K];

    // loop over max_iter
    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
        //std::cout << "iteration " << iter << "..." << std::endl;

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
        //std::cout << "found closest distance to each centroid, and the centr-centr distances" << std::endl;


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
                    if(elkan_lower_bounds[folan][filan] < second_smallest && filan != labels[folan]){
                        second_smallest = elkan_lower_bounds[folan][filan];
                    }
                }
                hamerly_upper_bounds[folan] = smallest;
                hamerly_lower_bounds[folan] = second_smallest;
            }
            //std::cout << "calculated distances" << std::endl;

            //std::cout << "filled out hamerly ub and lb and elkan lb for first time" << std::endl;

        } else {
            for (folan = 0; folan < N; folan++) {
                r = true;
                hamerly_bound = ((0.5 * closest_centroid_distance[labels[folan]]) > hamerly_lower_bounds[folan]) ? (
                        0.5 * closest_centroid_distance[labels[folan]]) : hamerly_lower_bounds[folan];
                // elkan lemma 1 and hamerly
                if (hamerly_bound < hamerly_upper_bounds[folan]) {
                    /*if(folan == 1){
                        std::cout << "not pruned by ham_bound\n";
                    }*/

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
                    /*if(folan == 1){
                        std::cout << "1 changed ham_ub[" << folan << "] to " << hamerly_upper_bounds[folan];
                        std::cout << " and elkan_lb[" << folan << "][" << labels[folan] << "] to " << elkan_lower_bounds[folan][labels[folan]] << std::endl;
                    }*/
                    // if((0.5 * closest_centroid_distance[labels[folan]]) < hamerly_upper_bounds[folan]){
                    for (filan = 0; filan < K; filan++) {
                        if (filan == labels[folan]) continue;
                        
                        if (hamerly_upper_bounds[folan] > elkan_lower_bounds[folan][filan]
                            && hamerly_upper_bounds[folan] >
                               (0.5 * centroid_to_centroid_distances[filan][labels[folan]])) {

                            tmp = centroid_squares[filan] + data_arr_ss[folan][0];
                            for (ashghal = 0; ashghal < D; ashghal++) {
                                tmp -= (2 * data_arr[folan][ashghal] * centroids[filan][ashghal]);
                            }
                            if(tmp < 0.0) tmp = 0.0;
                            distances[folan][filan] = sqrt(tmp);
                            elkan_lower_bounds[folan][filan] = distances[folan][filan];
                            /*if(folan == 1){
                                std::cout << "2 changed elkan_lb[" << folan << "][" << filan << "] to " << elkan_lower_bounds[folan][filan] << std::endl;
                            }*/


                            if (distances[folan][filan] < distances[folan][labels[folan]]) {
                                // keep the second smallest
                                hamerly_lower_bounds[folan] = distances[folan][labels[folan]];
                                labels[folan] = filan;
                                // i am doing this under duress...
                                hamerly_upper_bounds[folan] = distances[folan][labels[folan]];
                                /*if(folan == 1){
                                    std::cout << "1 changed ham_ub[" << folan << "] to " << hamerly_upper_bounds[folan];
                                    std::cout << " and ham_lb[" << folan << "] to " << hamerly_lower_bounds[folan];
                                    std::cout << "and labels[" << folan << "] to " << labels[folan] << std::endl;
                                }*/
                            } else if (hamerly_lower_bounds[folan] > distances[folan][filan]) {
                                hamerly_lower_bounds[folan] = distances[folan][filan];
                                /*if(folan == 1){
                                    std::cout << "1 changed ham_lb[" << folan << "] to " << hamerly_lower_bounds[folan] << std::endl;
                                }*/

                            }

                        }
                        /*else{
                            if(folan == 1){
                                std::cout << "pruned elkan " << filan << "\n";
                            }
                        }*/
                    }
                    // TODO: this is from the hamerly imp, im not sure if i should keep it or not, for now let's cmnt it
                    // DONE: I remembered why, because I don't want to assign it
                    // but the cpp code does it too, so for now ill do it mult times till i debug properly
                    // hamerly_upper_bounds[folan] = distances[folan][labels[folan]];
                }
                /*else{
                    if(folan == 1){
                        std::cout << "pruned by ham_bound \n";
                    }
                }*/
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

        //std::cout << "set labels" << std::endl;




        // calc new centroids
        // make copy of old centroids
        // apparently this also does not work, copies pointer somehow I think...
        // memcpy(old_centroids, centroids, sizeof(double) * K * D);
        for (folan = 0; folan < K; folan++) {
            for (filan = 0; filan < D; filan++) {
                old_centroids[folan][filan] = centroids[folan][filan];
            }
        }
        //std::cout << "copied centroids to old centroids" << std::endl;
        // set centroids to 0
        memset(cluster_counts, 0, sizeof(int) * K);
        //std::cout << "set cluster counts to 0" << std::endl;
        for (folan = 0; folan < K; folan++) {
            // just testing
            // cluster_counts[folan] = 0;
            for (filan = 0; filan < D; filan++) {
                centroids[folan][filan] = 0.0;
            }
        }
        // This doesn't work on doubles
        // memset(centroids, 0, sizeof(double) * K * D);
        //std::cout << "after all the memcpys" << std::endl;

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
        //std::cout << "calculated new centroids" << std::endl;
        // just to check
        /*int sanity_check = 0;
        std::cout << "cluster counts..." << std::endl;
        for (folan = 0; folan < K; folan++) {
            sanity_check += cluster_counts[folan];
            std::cout << cluster_counts[folan] << " ";
        }
        std::cout << sanity_check << std::endl;*/

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
        //std::cout << "calculated centroid movements" << std::endl;

        
        /*std::cout << "ham_ub[1] " << hamerly_upper_bounds[1] << " ham_lb[1] " << hamerly_lower_bounds[1] << std::endl;
        for(filan = 0; filan < K; filan++){
            std::cout << elkan_lower_bounds[1][filan] << " ";
        }
        std::cout << std::endl;*/

        //std::cout << "labels[1] " << labels[1] << std::endl;

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

        /*std::cout << "ham_ub[1] " << hamerly_upper_bounds[1] << " ham_lb[1] " << hamerly_lower_bounds[1] << std::endl;
        for(filan = 0; filan < K; filan++){
            std::cout << elkan_lower_bounds[1][filan] << " ";
        }
        std::cout << std::endl;*/



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
        //std::cout << "checked convergence" << std::endl;

        // end if converged
        if (has_converged) {
            if (D >= 256*256) {
                int sanity_check = 0;
                for (folan = 0; folan < K; folan++) {
                    sanity_check += cluster_counts[folan];
                    std::cout << cluster_counts[folan] << " ";
                }
                std::cout << sanity_check << std::endl;
            }
            break;
        }
    }
}


   // hamerly + step-wise
// copied from v5
// TODO: done
void kmeans_v9() {
    //std::cout << "in kmeans_v9 ..." << std::endl;
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
    //std::cout << "copied init centroids" << std::endl;


    bool has_converged = false;
    int cluster_counts[K];
    // V2 DIFF start
    double centroid_squares[K];
    // V2 DIFF end

    // loop over max_iter
    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
        //std::cout << "iteration " << iter << "..." << std::endl;

        // calculate the square sum of centroids
        // for (folan = 0; folan < K; folan++) {
        //     centroid_squares[folan] = 0;
        //     for (filan = 0; filan < D; filan++) {
        //         centroid_squares[folan] += (centroids[folan][filan] * centroids[folan][filan]);
        //     }
        // }

        // calculate closest centroid to each centroid
        for (folan = 0; folan < K; folan++) {
            smallest = DBL_MAX;
            for (filan = 0; filan < K; filan++) {
                tmp = centroids_ss[filan][0] + centroids_ss[folan][0];
                for (ashghal = 0; ashghal < D; ashghal++) {
                    tmp -= (2 * centroids[folan][ashghal] * centroids[filan][ashghal]);
                }
                if (tmp < smallest) smallest = tmp;
            }
            if(smallest < 0.0) smallest = 0.0;
            closest_centroid_distance[folan] = sqrt(smallest);
        }
        //std::cout << "found closest distance to each centroid" << std::endl;


        if (iter == 0) {
            //std::cout << "in if iter == 0" << std::endl;
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
            memset(assigned, 0, sizeof(int) * N);
            calculate_labels_with_sqrt_hamerly_integrated();
            // V9
            //std::cout << "calculated distances" << std::endl;
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
            //std::cout << "filled out hamerly ub and lb for first time" << std::endl;

        } else {
            //std::cout << "iter was not 0" << std::endl;
            // V9
            // fill out mask
            // set all to 1 first
            for (folan = 0; folan < N; folan++) {
                for (filan = 0; filan < K; filan++) {
                    is_candidate[folan][filan] = true;
                }
            }

            memset(assigned, 0, sizeof(int) * N);
            // now let's do pruning based on TI
            for (folan = 0; folan < N; folan++) {
                // hamerly check
                hamerly_bound = ((0.5 * closest_centroid_distance[labels[folan]]) > hamerly_lower_bounds[folan]) ? (
                        0.5 * closest_centroid_distance[labels[folan]]) : hamerly_lower_bounds[folan];
                if (hamerly_bound >= hamerly_upper_bounds[folan]) {
                    hamerly_count += K;
                    // this one will not move, so none of them are candidates
                    assigned[folan] = 1;
                    // for (filan = 0; filan < K; filan++) {
                    //     is_candidate[folan][filan] = false;
                    // }
                    // TODO: we may have to leave the labels[folan] = true, but i don't think so
                    // I will leave it as false for now and change it if necessary
                }
            }
            //std::cout << "hamerly pruned " << hamerly_count << std::endl;
            
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

        //std::cout << "set labels" << std::endl;




        // calc new centroids
        // make copy of old centroids
        // apparently this also does not work, copies pointer somehow I think...
        // memcpy(old_centroids, centroids, sizeof(double) * K * D);
        for (folan = 0; folan < K; folan++) {
            for (filan = 0; filan < D; filan++) {
                old_centroids[folan][filan] = centroids[folan][filan];
            }
        }
        //std::cout << "copied centroids to old centroids" << std::endl;
        // set centroids to 0
        memset(cluster_counts, 0, sizeof(int) * K);
        //std::cout << "set cluster counts to 0" << std::endl;
        for (folan = 0; folan < K; folan++) {
            // just testing
            // cluster_counts[folan] = 0;
            for (filan = 0; filan < D; filan++) {
                centroids[folan][filan] = 0.0;
            }
        }
        // This doesn't work on doubles
        // memset(centroids, 0, sizeof(double) * K * D);
        //std::cout << "after all the memcpys" << std::endl;

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
        //std::cout << "calculated new centroids" << std::endl;
        calculate_centroids_square_sums();

        // just to check
        /*int sanity_check = 0;
        std::cout << "cluster counts..." << std::endl;
        for (folan = 0; folan < K; folan++) {
            sanity_check += cluster_counts[folan];
            std::cout << cluster_counts[folan] << " ";
        }
        std::cout << sanity_check << std::endl;*/

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
        //std::cout << "calculated centroid movements" << std::endl;

        // update upper and lower hamerly bounds based on centroid movements
        for (folan = 0; folan < N; folan++) {
            //if (folan == 0) std::cout << "moving ub " << centroid_movement[labels[folan]] << std::endl;
            hamerly_upper_bounds[folan] += centroid_movement[labels[folan]];
            if (labels[folan] == furthest_moving_centroid) {
                //if (folan == 0) std::cout << "moving lb " << centroid_movement[second_furthest_moving_centroid] << std::endl;
                hamerly_lower_bounds[folan] -= centroid_movement[second_furthest_moving_centroid];
            } else {
                //if (folan == 0) std::cout << "moving lb " << centroid_movement[furthest_moving_centroid] << std::endl;
                hamerly_lower_bounds[folan] -= centroid_movement[furthest_moving_centroid];
            }
        }


        // FATEMEH DEBUG
        //std::cout << "after updating the bounds, hamerly_ub[0] " << hamerly_upper_bounds[0] << " hamerly_lb[0] "
        //     << hamerly_lower_bounds[0] << std::endl;
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
        //std::cout << "checked convergence" << std::endl;

        // end if converged
        if (has_converged) {
                if (D >= 256*256) {
                    int sanity_check = 0;
                    for (folan = 0; folan < K; folan++) {
                        sanity_check += cluster_counts[folan];
                        std::cout << cluster_counts[folan] << " ";
                    }
                    std::cout << sanity_check << std::endl;
                }
            break;
        }
    }

    //std::cout << "hamerly pruned " << hamerly_count << std::endl;
}

    int init(int argc, char **argv) {

        //std::cout << " " << D << " "  << K << " "  << N << std::endl;
        // allocate all arrays
        labels = (int *) malloc(N * sizeof(int));
        for (int i = 0; i < N; i++) {
            labels[i] = 0; // (double *) malloc(D * sizeof(double));
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


        //std::cout << "after allocation" << std::endl;

        

        //std::cout << "read data..." << std::endl;

        smallest_ub = (int *) malloc(N * sizeof(int));

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
        hamerly_lower_bounds = (double *) malloc(N * sizeof(double));
        for (int i = 0; i < N; i++) {
            hamerly_lower_bounds[i] = 0; // (double *) malloc(K * sizeof(double));
        }

        hamerly_upper_bounds = (double *) malloc(N * sizeof(double));
        for (int i = 0; i < N; i++) {
            hamerly_upper_bounds[i] = 0; // (double *) malloc(K * sizeof(double));
        }

        closest_centroid_distance = (double *) malloc(K * sizeof(double));
        for (int i = 0; i < K; i++) {
            closest_centroid_distance[i] = 0; // (double *) malloc(K * sizeof(double));
        }

        centroid_movement = (double *) malloc(K * sizeof(double));
        for (int i = 0; i < K; i++) {
            centroid_movement[i] = 0; // (double *) malloc(K * sizeof(double));
        }
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


        // MARIGOLG EXTRAS
        last_level_calculated = (int **) malloc(N * sizeof(int *));
        for (int i = 0; i < N; i++) {
            last_level_calculated[i] = (int *) malloc(K * sizeof(int));
        }
        // MARIGOLG EXTRAS end

        // just helper
        assigned = (int *) malloc(N * sizeof(int));
        // just helper

        //std::cout << "filled out the ss arr..." << std::endl;


        // set labels to 0 init state
        //memset(labels, 0, sizeof(int) * N);

        /*std::cout << "set labels to 0, calling kmeans..." << std::endl;

        // do the clustering
        //std::time_t t1 = std::time(nullptr);
        auto t1 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
        kmeans_v4();
        auto t2 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
        
        
        std::cout << t2-t1 << " nanoseconds " << (t2-t1)/1000000 << " milliseconds " << (t2-t1)/1000000000 << " seconds\n";
        

        // write labels to somewhere I guess...
        // TODO
        std::ofstream label_file;
        label_file.open(label_file_name.c_str());

        for (int i = 0; i < N; i++) {
            label_file << labels[i] << "\n";
        }
        label_file.close();*/


        return 0;
    }

    
    void clean() {
        
        //std::cout << " " << D << " "  << K << " "  << N << std::endl;
        // allocate all arrays
        free(labels);

        

        
        for (int i = 0; i < K; i++) {
            free(centroids[i]);
        }
        free(centroids);

        
        for (int i = 0; i < K; i++) {
            free(old_centroids[i]);
        }
        free(old_centroids);

        
        for (int i = 0; i < N; i++) {
            free(distances[i]);
        }
        free(distances);



        //std::cout << "after allocation" << std::endl;

        free(smallest_ub);// = (int *) malloc(N * sizeof(int));

        //std::cout << "read data..." << std::endl;

        // STEPWISE EXTRAS
        for (int i = 0; i < N; i++) {
            free(is_candidate[i]);
    
        }
        free(is_candidate);

        for (int i = 0; i < N; i++) {
            free(upper_bounds[i]);
        }
        free(upper_bounds);

        for (int i = 0; i < N; i++) {
            free(lower_bounds[i]);
        }
        free(lower_bounds);

        for (int i = 0; i < N; i++) {
            free(incremental_dots[i]);
        }
        free(incremental_dots);

        for (int i = 0; i < N; i++) {
            free(data_arr_ss[i]);
        }
        free(data_arr_ss); 

        
        for (int i = 0; i < K; i++) {
            free(centroids_ss[i]);
        }
        free(centroids_ss);
        // let's check if these are correct

        // STEPWISE EXTRAS end

        // HAMERLY EXTRAS
        free(hamerly_lower_bounds); 
        free(hamerly_upper_bounds);
        free(closest_centroid_distance);
        free(centroid_movement);
        // HAMERLY EXTRAS end

        // ELKAN EXTRAS
        for (int i = 0; i < N; i++) {
            free(elkan_lower_bounds[i]);
        }
        free(elkan_lower_bounds);

        for (int i = 0; i < K; i++) {
            free(centroid_to_centroid_distances[i]);
        }
        free(centroid_to_centroid_distances);
        // ELKAN EXTRAS end


        // MARIGOLG EXTRAS
        for (int i = 0; i < N; i++) {
            free(last_level_calculated[i]);
        }
        free(last_level_calculated);
        // MARIGOLG EXTRAS end

        // just helper
        free(assigned);
        // just helper

        //std::cout << "filled out the ss arr..." << std::endl;


        // set labels to 0 init state
       // memset(labels, 0, sizeof(int) * N);

        /*std::cout << "set labels to 0, calling kmeans..." << std::endl;

        // do the clustering
        //std::time_t t1 = std::time(nullptr);
        auto t1 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
        kmeans_v4();
        auto t2 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
        
        
        std::cout << t2-t1 << " nanoseconds " << (t2-t1)/1000000 << " milliseconds " << (t2-t1)/1000000000 << " seconds\n";
        

        // write labels to somewhere I guess...
        // TODO
        std::ofstream label_file;
        label_file.open(label_file_name.c_str());

        for (int i = 0; i < N; i++) {
            label_file << labels[i] << "\n";
        }
        label_file.close();*/

    }

// STUFF FOR COMPLETE MARIGOLD

// assume it is candoidate and stuff, and the inc_dot is complete till here
// then calc the ub and lb of distance between data[folan] and cent[filan] till level
void calculate_distance_folan_filan_till_level(int folan, int filan, int level){
    /*if(folan == 4){
        std::cout << "called folan_filan  for " << folan << " " << filan << " " << level << std::endl;
    }*/
    
    // check to not redo anything
    if(last_level_calculated[folan][filan] >= level){ 
        /*if(folan == 4){
            std::cout << "last level was " << last_level_calculated[folan][filan] << " so we are returning" << std::endl;
        }*/
        return;
    }
    if(last_level_calculated[folan][filan] != level - 1 && folan == 4){
        std::cout << "DISTASTER!!!!!!!!!!!!!!!!! " << folan << " " << filan << " " << level << " " << last_level_calculated[folan][filan] << std::endl;
        std:: cout << "fingers crossed..." << std::endl;
        // exit(2);
    }

    if(last_level_calculated[folan][filan] == -1){
        incremental_dots[folan][filan] = 0.0;
    }
    

    int ashghal, halghe;
    int d_sqrt = sqrt(D);
    // int two_p_level_m1 = int(pow(2, level - 1));
    int two_p_level_m1 = int(pow(2, last_level_calculated[folan][filan]));
    int two_p_level = std::min(int(pow(2, level)), d_sqrt);
    double this_dot = 0.0;
    double tmp_ub, tmp_lb;


    tmp_ub = centroids_ss[filan][0] + data_arr_ss[folan][0];
    tmp_lb = centroids_ss[filan][0] + data_arr_ss[folan][0];

    // debug
    // if (level == 0) {
    //     incremental_dots[folan][filan] = data_arr[folan][0] * centroids[filan][0];
    // } else {
    //     for(ashghal = 0; ashghal < two_p_level_m1; ashghal++){
    //         // amoodi ha
    //         for(halghe = ashghal * d_sqrt + two_p_level_m1; halghe < ashghal * d_sqrt + two_p_level; halghe++){
    //             incremental_dots[folan][filan] += (data_arr[folan][halghe] * centroids[filan][halghe]);
    //         }
    //     }
    //     for(ashghal = 0; ashghal < std::min(two_p_level_m1, d_sqrt - two_p_level_m1); ashghal++){
    //         // ofoghi ha
    //         for(halghe = d_sqrt * (two_p_level_m1 + ashghal); halghe < d_sqrt * (two_p_level_m1 + ashghal) + two_p_level; halghe++){
    //             incremental_dots[folan][filan] += (data_arr[folan][halghe] * centroids[filan][halghe]);
    //         }
    //     }
    // }

    // NEW FORMATION to handle starting from 2^last_level
    /*if(folan == 4){
        std::cout << "amoodi ha" << std::endl;
    }*/
    for(ashghal = 0; ashghal < two_p_level_m1; ashghal++){
        // amoodi ha
        for(halghe = ashghal * d_sqrt + two_p_level_m1; halghe < ashghal * d_sqrt + two_p_level; halghe++){
            /*if(folan == 4){
                std::cout << halghe << " ";
            }*/
            incremental_dots[folan][filan] += (data_arr[folan][halghe] * centroids[filan][halghe]);
        }
        /*if(folan == 4){
        std::cout << std::endl;
        }*/
    }
    // for(ashghal = 0; ashghal < std::min(two_p_level_m1, d_sqrt - two_p_level_m1); ashghal++){
    /*if(folan == 4){
        std::cout << "ofoghi ha" << std::endl;
    }*/
    for(ashghal = 0; ashghal < std::min(two_p_level - two_p_level_m1, d_sqrt - two_p_level_m1); ashghal++){
        // ofoghi ha
        for(halghe = d_sqrt * (two_p_level_m1 + ashghal); halghe < d_sqrt * (two_p_level_m1 + ashghal) + two_p_level; halghe++){
            /*if(folan == 4){
                std::cout << halghe << " ";
            }*/
            incremental_dots[folan][filan] += (data_arr[folan][halghe] * centroids[filan][halghe]);
        }
        /*if(folan == 4){
        std::cout << std::endl;
        }*/
    }
    tmp_ub -= 2 * incremental_dots[folan][filan];
    tmp_lb -= 2 * incremental_dots[folan][filan];

    // this_dot = 0.0;
    // for(int i = 0; i < two_p_level; i++){
    //     for(int j = 0; j < two_p_level; j++){
    //         this_dot += (data_arr[folan][d_sqrt*i + j] * centroids[filan][d_sqrt*i + j]);
    //     }
    // }
    // // upper_bounds[folan][filan] -= (2 * this_dot);
    // // lower_bounds[folan][filan] -= (2 * this_dot);
    // tmp_ub -= (2 * this_dot);
    // tmp_lb -= (2 * this_dot);

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

    // otherwise, increment to know we have done it
    last_level_calculated[folan][filan] = level;
    /*if(folan == 4){
        std::cout << "went though successfully, new last-level is " << last_level_calculated[folan][filan] << std::endl;
    }*/
}

// copied from calculate_sqrt_distances_till_level
// just checks assigned instead of labels
void calculate_sqrt_distances_till_level_with_assigned_and_ll(int level) {
    int folan, filan, ashghal, halghe;
    int d_sqrt = sqrt(D);
    // int two_p_level_m1 = int(pow(2, level - 1));
    int two_p_level_m1 = int(pow(2, last_level_calculated[folan][filan]));
    int two_p_level = std::min(int(pow(2, level)), d_sqrt);
    double this_dot = 0.0;
    double tmp_ub, tmp_lb;

    for (folan = 0; folan < N; folan++) {

        if (assigned[folan] > 0) continue;

        for (filan = 0; filan < K; filan++) {

            if (is_candidate[folan][filan] == false) continue;

            two_p_level_m1 = int(pow(2, last_level_calculated[folan][filan]));

            // known parts

            tmp_ub = centroids_ss[filan][0] + data_arr_ss[folan][0];
            tmp_lb = centroids_ss[filan][0] + data_arr_ss[folan][0];

            if(last_level_calculated[folan][filan] == -1){
                incremental_dots[folan][filan] = 0.0;
            }
            if(last_level_calculated[folan][filan] < level){
                // NEW FORMATION to handle starting from 2^last_level
                for(ashghal = 0; ashghal < two_p_level_m1; ashghal++){
                    // amoodi ha
                    for(halghe = ashghal * d_sqrt + two_p_level_m1; halghe < ashghal * d_sqrt + two_p_level; halghe++){
                        incremental_dots[folan][filan] += (data_arr[folan][halghe] * centroids[filan][halghe]);
                    }
                }
                
                for(ashghal = 0; ashghal < std::min(two_p_level - two_p_level_m1, d_sqrt - two_p_level_m1); ashghal++){
                    // ofoghi ha
                    for(halghe = d_sqrt * (two_p_level_m1 + ashghal); halghe < d_sqrt * (two_p_level_m1 + ashghal) + two_p_level; halghe++){
                        incremental_dots[folan][filan] += (data_arr[folan][halghe] * centroids[filan][halghe]);
                    }
                }
                last_level_calculated[folan][filan] = level;
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

// let's try for the whole thing
// hamerly + elkan + new_lemma + step-wise
// TODO
void kmeans_v10(){

    int folan = 0;
    int filan = 0; int alaki;
    int ashghal = 0;// for loop usage
    int furthest_moving_centroid, second_furthest_moving_centroid;
    double smallest, second_smallest;
    double tmp, hamerly_bound;
    int hamerly_count = 0;
    bool r; int level;
    int r_int; double fake_smallest_lb; bool candidates_exist;
    //int *smallest_ub = (int *) malloc(N * sizeof(int));

    // set initial centroids

    for(folan = 0; folan < K; folan++){
        for(filan = 0; filan < D; filan++){
            centroids[folan][filan] = data_arr[folan][filan];
        }
    }

    // memcpy(centroids, data_arr, sizeof(double) * K * D);
    calculate_centroids_square_sums();
    //std::cout << "copied init centroids" << std::endl;


    bool has_converged = false;
    int cluster_counts[K];
    

    // let's do the first iteration out here, so the if inside the loop is resolved
    // TODO
    level = 0;
    // calculate closest centroid to each centroid
    for (folan = 0; folan < K; folan++) {
        smallest = DBL_MAX;
        for (filan = 0; filan < K; filan++) {
            tmp = centroids_ss[filan][0] + centroids_ss[folan][0];
            for (ashghal = 0; ashghal < D; ashghal++) {
                tmp -= (2 * centroids[folan][ashghal] * centroids[filan][ashghal]);
            }
            if (tmp < smallest) smallest = tmp;
        }
        if(smallest < 0.00000) smallest = 0.0;
        closest_centroid_distance[folan] = sqrt(smallest);
    }
    //std::cout << "found closest distance to each centroid" << std::endl;


    // fill out mask
    // set all to 1 first
    for (folan = 0; folan < N; folan++) {
        for (filan = 0; filan < K; filan++) {
            is_candidate[folan][filan] = true;
        }
    }

    //std::cout << "filled out is_candid" << std::endl;

    // set the assigned
    memset(assigned, 0, sizeof(int) * N);
    memset(smallest_ub, 0, sizeof(int) * N);

    //std::cout << "memset out assigned and smallest_ub" << std::endl;

    for(folan = 0; folan < N; folan++){
        for(filan=0; filan <K; filan++){
            last_level_calculated[folan][filan] = -1;
        }
    }
    //std::cout << "set last_level_calc to -1" << std::endl;

    // I think I have to put the step-wise while here

    while (level < int(log2(int(sqrt(D))) + 2)) {
        //std::cout << "iter 0, in while level =  " << level << std::endl;
        for(folan = 0; folan < N; folan++){
            if(assigned[folan] > 0) continue;
            for(filan = 0; filan < K; filan++){
                if(!is_candidate[folan][filan]) continue;
                calculate_distance_folan_filan_till_level(folan, filan, level);
                // fill elkan lb
                elkan_lower_bounds[folan][filan] = lower_bounds[folan][filan];
                // find smallest_ub
                if(upper_bounds[folan][filan] < upper_bounds[folan][smallest_ub[folan]]) smallest_ub[folan] = filan;
            }
            // std::cout << "doing pruning for point " << folan << " ..." << std::endl;
            // TODO: do step-wise pruning here: DONE
            // also fill out the ham_lb: DONE and elkan_lb: DONE somewhere
            // for ham_lb
            // fake_smallest_lb = lower_bounds[folan][smallest_ub[folan]];
            // fake_smallest_lb = DBL_MAX;
            // everybody is guilty until proven innocent...
            candidates_exist = false;
            if (level == int(log2(int(sqrt(D))) + 1)) {
                labels[folan] = smallest_ub[folan];
                // V9
                hamerly_upper_bounds[folan] = upper_bounds[folan][smallest_ub[folan]];
                assigned[folan] = 1;
                // v9
                // but we still need to update the ham_lb
                for (filan = 0; filan < K; filan++) {
                    if (filan == smallest_ub[folan]) continue;
                    if (!is_candidate[folan][filan]) continue;
                    if(!candidates_exist){
                        candidates_exist = true;
                        // fake_smallest_lb = DBL_MAX;
                        hamerly_lower_bounds[folan] = DBL_MAX;
                    }
                    // if(folan == 0){
                    //     std::cout << "before: HERE: filan: " << filan << " lower_bounds[folan][filan] " << lower_bounds[folan][filan] << " = " << elkan_lower_bounds[folan][filan] << std::endl;
                    //     std::cout << "HERE: fake_smallest_ub: " << fake_smallest_lb << " ham_lb[folan] " << hamerly_lower_bounds[folan] << std::endl;
                    //     std::cout << "----------------" << std::endl;
                    // }
                    
                    // find sec smallest lb
                    // if(lower_bounds[folan][filan] < fake_smallest_lb){
                    //     // keep sec smallest
                    //     hamerly_lower_bounds[folan] = fake_smallest_lb;
                    //     fake_smallest_lb = lower_bounds[folan][filan];
                    // } else if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                    //     hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                    // }

                    // TODO: so im changing the definition of ham_lb
                    // let ham_lb be the smallest lb betweent the last batch of candidates
                    // if there are other candidates at THIS point, then the ham_lb should just be fake_smallest_lb here
                    // but what if there are no other candidates?
                    // then we would have set fake_smallest_lb to DBL_MAX and it stays that way...
                    // I fixed it with the exists_candidate
                    // it's not the label == smallest_ub
                    // so it's the smallest lb to get pruned
                    if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                        hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                    }


                    // if(folan == 0){
                    //     std::cout << "after: HERE: filan: " << filan << " lower_bounds[folan][filan] " << lower_bounds[folan][filan] << " = " << elkan_lower_bounds[folan][filan] << std::endl;
                    //     std::cout << "HERE: fake_smallest_ub: " << fake_smallest_lb << " ham_lb[folan] " << hamerly_lower_bounds[folan] << std::endl;
                    //     std::cout << "----------------" << std::endl;
                    // }

                }
            } else {
                for (filan = 0; filan < K; filan++) {
                    if (filan == smallest_ub[folan]) continue;
                    if (!is_candidate[folan][filan]) continue;
                    if(!candidates_exist){
                        candidates_exist = true;
                        // fake_smallest_lb = DBL_MAX;
                        hamerly_lower_bounds[folan] = DBL_MAX;
                    }
                    if (lower_bounds[folan][filan] >= upper_bounds[folan][smallest_ub[folan]]) {
                        is_candidate[folan][filan] = false;
                        // V9
                        // we should also try to find the second closest centroid
                        // TODO: im not sure if i should clear the old ham_lbs from last round
                        // what if it's masked completely, then we would never get here, so the old one should still be valid.
                        // nvm, i think im right, i won't clear it
                        // if (lower_bounds[folan][filan] < hamerly_lower_bounds[folan]) {
                        //     hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                        // }
                        // V9
                        // I'm fairly sure the above is bullshit, can't remember what I was thinking
                        // so ill make another
                        
                        
                    }
                    /*if(folan == 0){
                        std::cout << "before: HERE: filan: " << filan << " lower_bounds[folan][filan] " << lower_bounds[folan][filan] << " = " << elkan_lower_bounds[folan][filan] << std::endl;
                        std::cout << "HERE: ham_lb[folan] " << hamerly_lower_bounds[folan] << std::endl;
                        std::cout << "----------------" << std::endl;
                    }*/
                    
                    // find sec smallest lb
                    // if(lower_bounds[folan][filan] < fake_smallest_lb){
                    //     // keep sec smallest
                    //     hamerly_lower_bounds[folan] = fake_smallest_lb;
                    //     fake_smallest_lb = lower_bounds[folan][filan];
                    // } else if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                    //     hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                    // }
                    
                    // change of ham_lb definition

                    if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                        hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                    }

                    /*if(folan == 0){
                        std::cout << "after: HERE: filan: " << filan << " lower_bounds[folan][filan] " << lower_bounds[folan][filan] << " = " << elkan_lower_bounds[folan][filan] << std::endl;
                        std::cout << "HERE: ham_lb[folan] " << hamerly_lower_bounds[folan] << std::endl;
                        std::cout << "----------------" << std::endl;
                    }*/
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

    //std::cout << "set labels " << std::endl;
        

    // sanity check: is everyone assigned?
    /*for(folan = 0; folan < N; folan++){
        if(assigned[folan] == 0){
            std::cout << folan << " NOT ASSIGNED" << std::endl;
            exit(2);
        }
    }*/

    //std::cout << "after all-assigned sanity check" << std::endl;

    // calc new centroids
    // make copy of old centroids
    // apparently this also does not work, copies pointer somehow I think...
    // memcpy(old_centroids, centroids, sizeof(double) * K * D);
    for (folan = 0; folan < K; folan++) {
        for (filan = 0; filan < D; filan++) {
            old_centroids[folan][filan] = centroids[folan][filan];
        }
    }
    //std::cout << "copied centroids to old centroids" << std::endl;
    // set centroids to 0
    memset(cluster_counts, 0, sizeof(int) * K);
    //std::cout << "set cluster counts to 0" << std::endl;
    for (folan = 0; folan < K; folan++) {
        // just testing
        // cluster_counts[folan] = 0;
        for (filan = 0; filan < D; filan++) {
            centroids[folan][filan] = 0.0;
        }
    }
    // This doesn't work on doubles
    // memset(centroids, 0, sizeof(double) * K * D);
    //std::cout << "after all the memcpys" << std::endl;

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
    //std::cout << "calculated new centroids" << std::endl;
    calculate_centroids_square_sums();
    //std::cout << "ssed them " << std::endl;
    // just to check
    //int sanity_check = 0;
    /*std::cout << "cluster counts..." << std::endl;
    for (folan = 0; folan < K; folan++) {
        sanity_check += cluster_counts[folan];
        std::cout << cluster_counts[folan] << " ";
    }
    std::cout << sanity_check << std::endl;
    */
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
    //std::cout << "calculated centroid movements" << std::endl;

    /*    
    std::cout << "ham_ub[1] " << hamerly_upper_bounds[1] << " ham_lb[1] " << hamerly_lower_bounds[1] <<  std::endl;
    for(filan =0; filan < K; filan++){
        std::cout << elkan_lower_bounds[1][filan] << " ";
    }
    std::cout << std::endl;
    */

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

    /*std::cout << "after updating the bounds...\n";

    std::cout << "ham_ub[1] " << hamerly_upper_bounds[1] << " ham_lb[1] " << hamerly_lower_bounds[1] <<  std::endl;
    for(filan =0; filan < K; filan++){
        std::cout << elkan_lower_bounds[1][filan] << " ";
    }
    std::cout << std::endl;*/




    

    // loop over max_iter
    for (int iter = 1; iter < MAX_ITERATIONS; iter++) {
        level = 0;
        std::cout << "iteration " << iter << "..." << std::endl;

        // calculate closest centroid to each centroid
        for (folan = 0; folan < K; folan++) {
            smallest = DBL_MAX;
            for (filan = 0; filan < K; filan++) {
                tmp = centroids_ss[filan][0] + centroids_ss[folan][0];
                for (ashghal = 0; ashghal < D; ashghal++) {
                    tmp -= (2 * centroids[folan][ashghal] * centroids[filan][ashghal]);
                }
                if (tmp < smallest) smallest = tmp;
            }
            if(smallest < 0.00000) smallest = 0.0;
            closest_centroid_distance[folan] = sqrt(smallest);
        }
        //std::cout << "found closest distance to each centroid" << std::endl;


        // fill out mask
        // set all to 1 first
        for (folan = 0; folan < N; folan++) {
            for (filan = 0; filan < K; filan++) {
                is_candidate[folan][filan] = true;
            }
        }
        //std::cout << "filled out is_candid" << std::endl;


        // set the assigned
        memset(assigned, 0, sizeof(int) * N);
        //std::cout << "memset assigned" << std::endl;

        for(folan = 0; folan < N; folan++){
            for(filan=0; filan <K; filan++){
                last_level_calculated[folan][filan] = -1;
            }
        }
        //std::cout << "set last_level_calc to -1" << std::endl;

        // I think I have to put the step-wise while here

        while (level < int(log2(int(sqrt(D))) + 2)) {
            //std::cout << "in while level " << level << std::endl;
            // set smallest ub
            // let's initiate them to the previous labels so skipping the labels in the for K does not change anything
            // TODO!!!!
            // memset(smallest_ub, 0, sizeof(int) * N);


            for (folan = 0; folan < N; folan++) {
                // TODO: check if sth goes wrong
                smallest_ub[folan] = labels[folan];
                if (assigned[folan] > 0) continue;
                r = true;
                r_int = 0; fake_smallest_lb = DBL_MAX;
                hamerly_bound = ((0.5 * closest_centroid_distance[labels[folan]]) > hamerly_lower_bounds[folan]) ? (
                        0.5 * closest_centroid_distance[labels[folan]]) : hamerly_lower_bounds[folan];
                // elkan lemma 1 and hamerly
                if (hamerly_bound < hamerly_upper_bounds[folan]) {
                    /*if(folan == 4){
                        std::cout << "not pruned by ham_bound..." << std::endl;
                    }*/
                    // if((0.5 * closest_centroid_distance[labels[folan]]) < hamerly_upper_bounds[folan]){
                    for (filan = 0; filan < K; filan++) {
                        // TODO: tentative...
                        if(!is_candidate[folan][filan]) continue;
                        // if(folan == 1){
                        //     std::cout << "am i going crazy? " << filan << " " << labels[folan] << std::endl;
                        // }
                        if (filan == labels[folan]) continue;
                        // new lemma
                        if (elkan_lower_bounds[folan][filan] >
                            hamerly_upper_bounds[folan] + centroid_movement[labels[folan]] + centroid_movement[filan]) {
                                /*if(folan == 4){
                                    std::cout << "skipping for new lemma " << filan << std::endl;
                                }*/
                            // TODO: i think this should be a break instead of a continue, but im waiting for Mo's response
                            // DONE: I was wrong, it is just continue, cause it just means that filan will not be the correct cluster for this guy
                            continue;
                        }
                        /*if(folan == 4){
                            std::cout << "not pruned by new lemma " << filan << std::endl;
                        }*/
                        if (hamerly_upper_bounds[folan] > elkan_lower_bounds[folan][filan]
                            && hamerly_upper_bounds[folan] >
                                (0.5 * centroid_to_centroid_distances[filan][labels[folan]])) {
                                // HERE
                                // I'm thinking the full dist calculation of the prev lables[folan]
                                // we have to do it all in one go
                                // or the r bool has to be activated for all the levels or sth
                                
                                // if(r_int < int(log2(int(sqrt(D))) + 2)){
                                if(r){
                                    // update the upper bound
                                    // im gonna do it till this level only
                                    calculate_distance_folan_filan_till_level(folan, labels[folan], level);
                                    // update the ham_ub and elkan_lb if better
                                    // TODO: not sure, maybe we have to update it anyway
                                    //hamerly_upper_bounds[folan] = std::min(hamerly_upper_bounds[folan], upper_bounds[folan][labels[folan]]);
                                    //elkan_lower_bounds[folan][labels[folan]] = std::max(elkan_lower_bounds[folan][labels[folan]], lower_bounds[folan][labels[folan]]);
                                    // AS DISCUSSED WITH KASPER, lets's try it this way for now TODO: check
                                    hamerly_upper_bounds[folan] = upper_bounds[folan][labels[folan]];
                                    elkan_lower_bounds[folan][labels[folan]] = lower_bounds[folan][labels[folan]];
                                    /*if(folan == 4){
                                        std::cout << "JUST FOR THE SAKE OF FATEMEH'S SANITY: euc_dist(1, 1) " << sqrt(euclidean_distance(1, 1)) << std::endl; 
                                        std::cout << "1: changed ham_ub[" << folan << "] to " << hamerly_upper_bounds[folan];
                                        std::cout << " and elkan_lb[" << folan << "][" << labels[folan] << "] to " << elkan_lower_bounds[folan][labels[folan]] << std::endl;
                                    }*/
                                    // r_int++;
                                    r=false;

                                }
                                // TODO: maybe I should have a if(is_candid) around this
                                // we do the TI pruning with the ifs, then the step-wise prunings with the is_candid
                                // or I could just have a if(!is_candid) continue at the start of for K ?
                                // I think this could work
                                calculate_distance_folan_filan_till_level(folan, filan, level);
                                // i'm not sure about this maxing, maybe I should always set it to lb TODO?
                                elkan_lower_bounds[folan][filan] = lower_bounds[folan][filan];
                                /*if(folan == 4){
                                    std::cout << "2: changed  elkan_lb[" << folan << "][" << filan << "] to " << elkan_lower_bounds[folan][filan] << " ub was " << upper_bounds[folan][filan] << std::endl;
                                } */   
                                // elkan_lower_bounds[folan][filan] = std::max(elkan_lower_bounds[folan][filan], lower_bounds[folan][filan]);

                                // then (in non-step-wise) we check if the dist is smaller than the one for labels[folan]
                                // I'm not sure we can do that since we don't exactly have distances, but ub and lb
                                // we also don't have labels[folan], in the same way as before
                                // we could store a smallest_ub_till_now
                                // then use that as the labels[folan]
                                // we could even set the labels[folan] = smallest_ub_till_now and update as we go
                                // then we'd do the ic pruning here to avoid an extra loop
                                // no, that wouldn't really work. since the smallest_ub_till_now is not necessarily the smallest_ub over all
                                // we don't need a separate variable for smallest_ub, it could just be labels[folan]
                                // then the problem is when do we say it is assigned.
                                // TODO
                                // maybe I should update ham_lb separately
                                // like, if they are not from the same point, it would still be correct
                                // I think...
                                // yeah but then ham_ub==labels[folan] would be the smallest_ub
                                // but ham_lb is supposed to be the second smallest lb:-?
                                // or the second largest lb;-? 
                                // I don't think it should be sec largest lb
                                // i'll go with sec smallest lb, fingers crossed --> for this I will have to keep smallest_lb too
                                if(upper_bounds[folan][filan] < hamerly_upper_bounds[folan]){
                                    // TODO: how the f am i supposed to update ham_lb now?!
                                    // update the ham_lb if necessary
                                    // hamerly_lower_bounds[folan] = std::max(hamerly_lower_bounds[folan], hamerly_upper_bounds[folan]);
                                    // hamerly_lower_bounds[folan] = hamerly_upper_bounds[folan];
                                    labels[folan] = filan;
                                    hamerly_upper_bounds[folan] = upper_bounds[folan][filan];
                                    /*if(folan == 4){
                                        std::cout << "3: changed  ham_ub[" << folan << "] to " << hamerly_upper_bounds[folan];
                                        std::cout << " and labels[" << folan << "]  to " << labels[folan] << std::endl;
                                    }*/
                                }
                                // I have chnaged the def of ham_lb
                                // I will only update the ham_lb in the pruning of stepwise
                                // I only use the ham_lb in the ham_bound if, so I can just update it after the whole for K
                                // will it ever stay unassigned?
                                // after this for K
                                // we go into the step-wise pruning
                                // then, if there are candidates, the smallest lb of the candidates(other than the smallest ub) will be ham_lb
                                // then if there are no candidates, the ham_lb of last level remains
                                // the first level will always have candidates(everybody is a candidate)
                                // so we will always assign it
                                // no im doubting if the def is correct
                                // since the "the smallest lb of the candidates(other than the smallest ub)"
                                // does not really add up to much
                                // idk, let's see what happens:) 
                                // finding sec smallest lb
                                // if(lower_bounds[folan][filan] < fake_smallest_lb){
                                //     // TODO: maybe the max completely ruins the whole thing...
                                //     // i'm not gonna do the max for now
                                //     // hamerly_lower_bounds[folan] = std::max(hamerly_lower_bounds[folan], fake_smallest_lb);
                                //     hamerly_lower_bounds[folan] = fake_smallest_lb;
                                //     fake_smallest_lb = lower_bounds[folan][filan];
                                // } else if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                                //     hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                                // }
                        } 
                        // FOR NOW I WON'T
                        //else{
                            // TODO: I think I need to change the is_candid
                            /*if(folan == 4){
                                std::cout << "3.5 skipping this " << filan << std::endl;
                            }*/
                        //}

                        // im gonna update smallest_ub here
                        // I think it should be fine
                        // cause after all the hassle, the ub and lb here should be valid
                        if(upper_bounds[folan][filan] < upper_bounds[folan][smallest_ub[folan]]){ smallest_ub[folan] = filan; }
                    } // for filan < K
                } // if ham_bound
                else{
                    // TODO: I think I need to change the is_candid
                    //  I think  I should change the assigned...
                    assigned[folan] = 1;
                    // TODO: maybe I should set the ham_lb here too?:-?
                }
                // TODO: now we have to do the pruning based off step-wise... I think...
                // I'll move this to the start of the N loop
                // I added it there, but will keep it here because of the else that can change assigned with ham_boung
                if (assigned[folan] > 0) continue;
                /*if(folan == 4){
                    std::cout << "elkan_lb[4] ";
                    for(filan = 0; filan < K; filan++){
                        std::cout << elkan_lower_bounds[folan][filan] << " ";
                    }
                    std::cout << std::endl;
                }*/

                candidates_exist = false;
                // V9
                if (level == int(log2(int(sqrt(D))) + 1)) {
                    labels[folan] = smallest_ub[folan];
                    // V9
                    hamerly_upper_bounds[folan] = upper_bounds[folan][smallest_ub[folan]];
                    assigned[folan] = 1;
                    // v9

                    // but we still need to update the ham_lb
                    for (filan = 0; filan < K; filan++) {
                        if (filan == smallest_ub[folan]) continue;
                        if (!is_candidate[folan][filan]) continue;
                        if(!candidates_exist){
                            candidates_exist = true;
                            // fake_smallest_lb = DBL_MAX;
                            hamerly_lower_bounds[folan] = DBL_MAX;
                        }

                        // TODO: so im changing the definition of ham_lb
                        // let ham_lb be the smallest lb betweent the last batch of candidates
                        // if there are other candidates at THIS point, then the ham_lb should just be fake_smallest_lb here
                        // but what if there are no other candidates?
                        // then we would have set fake_smallest_lb to DBL_MAX and it stays that way...
                        // I fixed it with the exists_candidate
                        // it's not the label == smallest_ub
                        // so it's the smallest lb to get pruned
                        if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                            hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                        }

                    }

                } else {
                    for (filan = 0; filan < K; filan++) {
                        if (filan == smallest_ub[folan]) continue;
                        if (!is_candidate[folan][filan]) continue;
                        if(!candidates_exist){
                            candidates_exist = true;
                            // fake_smallest_lb = DBL_MAX;
                            hamerly_lower_bounds[folan] = DBL_MAX;
                        }
                        if (lower_bounds[folan][filan] >= upper_bounds[folan][smallest_ub[folan]]) {
                            is_candidate[folan][filan] = false;
                            /*if(folan == 4){
                                std::cout << "4: pruning " << filan << std::endl;
                            }*/
                            // NOT HERE i THINK
                            // V9
                            // we should also try to find the second closest centroid
                            // TODO: im not sure if i should clear the old ham_lbs from last round
                            // what if it's masked completely, then we would never get here, so the old one should still be valid.
                            // nvm, i think im right, i won't clear it
                            //if (lower_bounds[folan][filan] < hamerly_lower_bounds[folan]) {
                            //    hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                            //}
                            // V9
                        }
                        
                        if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                            hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                            /*if(folan == 4){
                                std::cout << "5: chnaged ham_lb[" << folan << "] to " << hamerly_lower_bounds[folan] << std::endl; 
                            }*/
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
                        /*if(folan == 4){
                            std::cout << "6: only one candidate chnaged labels[" << folan << "] to " << labels[folan];
                            std::cout << " and ham_ub[" << folan << "] to " << hamerly_upper_bounds[folan] << std::endl; 
                        }*/
                    }
                }
                
            } // for folan < N
            level++;
        } // step-wise loop
        
        // TODO: sanity check, is everybody assigned?
        /*for(folan = 0; folan < N; folan++){
            if(assigned[folan] == 0){
                std::cout << folan << " NOT ASSIGNED" << std::endl;
                exit(2);
            }
        }*/

        //std::cout << "set labels" << std::endl;


        // calc new centroids
        // make copy of old centroids
        // apparently this also does not work, copies pointer somehow I think...
        // memcpy(old_centroids, centroids, sizeof(double) * K * D);
        for (folan = 0; folan < K; folan++) {
            for (filan = 0; filan < D; filan++) {
                old_centroids[folan][filan] = centroids[folan][filan];
            }
        }
        //std::cout << "copied centroids to old centroids" << std::endl;
        // set centroids to 0
        memset(cluster_counts, 0, sizeof(int) * K);
        //std::cout << "set cluster counts to 0" << std::endl;
        for (folan = 0; folan < K; folan++) {
            // just testing
            // cluster_counts[folan] = 0;
            for (filan = 0; filan < D; filan++) {
                centroids[folan][filan] = 0.0;
            }
        }
        // This doesn't work on doubles
        // memset(centroids, 0, sizeof(double) * K * D);
        //std::cout << "after all the memcpys" << std::endl;

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
        //std::cout << "calculated new centroids" << std::endl;
        calculate_centroids_square_sums();
        //std::cout << "ssed them " << std::endl;
        // just to check
        //int sanity_check = 0;
        /*std::cout << "cluster counts..." << std::endl;
        for (folan = 0; folan < K; folan++) {
            sanity_check += cluster_counts[folan];
            std::cout << cluster_counts[folan] << " ";
        }
        std::cout << sanity_check << std::endl;
        */
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
            if(tmp < 0.0000) tmp = 0.0;
            centroid_movement[folan] = sqrt(tmp);
            if (centroid_movement[folan] > centroid_movement[furthest_moving_centroid]){
                second_furthest_moving_centroid = furthest_moving_centroid;
                furthest_moving_centroid = folan;
            }
            else if (centroid_movement[folan] >
                     centroid_movement[second_furthest_moving_centroid])
                second_furthest_moving_centroid = folan;
        }
        /*std::cout << "calculated centroid movements" << std::endl;

        std::cout << "ham_ub[1] " << hamerly_upper_bounds[1] << " ham_lb[1] " << hamerly_lower_bounds[1] << std::endl;
        for(filan = 0; filan < K; filan++){
            std::cout << elkan_lower_bounds[1][filan] << " " ;
        }
        std::cout << std::endl;*/


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
        //has_converged = (0.0 == centroid_movement[furthest_moving_centroid]);
        has_converged = true;
        for (folan = 0; folan < K; folan++) {
            for (filan = 0; filan < D; filan++) {
                if (old_centroids[folan][filan] != centroids[folan][filan]) {
                    has_converged = false;
                    break;
                }
            }
        }
        //std::cout << "checked convergence" << std::endl;

        // end if converged
        if (has_converged) {
            DEBUGPRINT("final iter: %d", iter);
            break;
        }
        

    } // for iter

}

// let's try for the whole thing
// hamerly + elkan + step-wise
// TODO
// let's try for the whole thing
// hamerly + elkan + step-wise
// TODO: done and tested
void kmeans_v105(){

    int folan = 0;
    int filan = 0; int alaki;
    int ashghal = 0;// for loop usage
    int furthest_moving_centroid, second_furthest_moving_centroid;
    double smallest, second_smallest;
    double tmp, hamerly_bound;
    int hamerly_count = 0;
    bool r; int level;
    int r_int; double fake_smallest_lb; bool candidates_exist;
    int *smallest_ub = (int *) malloc(N * sizeof(int));

    // set initial centroids

    for(folan = 0; folan < K; folan++){
        for(filan = 0; filan < D; filan++){
            centroids[folan][filan] = data_arr[folan][filan];
        }
    }

    // memcpy(centroids, data_arr, sizeof(double) * K * D);
    calculate_centroids_square_sums();
    //std::cout << "copied init centroids" << std::endl;


    bool has_converged = false;
    int cluster_counts[K];


    // let's do the first iteration out here, so the if inside the loop is resolved
    // TODO
    level = 0;
    // calculate closest centroid to each centroid
    for (folan = 0; folan < K; folan++) {
        smallest = DBL_MAX;
        for (filan = folan; filan < K; filan++) {
            tmp = centroids_ss[filan][0] + centroids_ss[folan][0];
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
    //std::cout << "found closest distance to each centroid, and the centr-centr distances" << std::endl;

    // fill out mask
    // set all to 1 first
    for (folan = 0; folan < N; folan++) {
        for (filan = 0; filan < K; filan++) {
            is_candidate[folan][filan] = true;
        }
    }

    //std::cout << "filled out is_candid" << std::endl;

    // set the assigned
    memset(assigned, 0, sizeof(int) * N);
    memset(smallest_ub, 0, sizeof(int) * N);

    //std::cout << "memset out assigned and smallest_ub" << std::endl;

    for(folan = 0; folan < N; folan++){
        for(filan=0; filan <K; filan++){
            last_level_calculated[folan][filan] = -1;
        }
    }
    //std::cout << "set last_level_calc to -1" << std::endl;

    // I think I have to put the step-wise while here
    
    // TODO:
    // for x in X:
    //  check ham_bound: if true, then make every entroid other than labels[x] not candidate


    while (level < int(log2(int(sqrt(D))) + 2)) {
        //std::cout << "iter 0, in while level =  " << level << std::endl;
        for(folan = 0; folan < N; folan++){
            if(assigned[folan] > 0) continue;
            for(filan = 0; filan < K; filan++){
                if(!is_candidate[folan][filan]) continue;
                calculate_distance_folan_filan_till_level(folan, filan, level);
                // fill elkan lb
                elkan_lower_bounds[folan][filan] = lower_bounds[folan][filan];
                // find smallest_ub
                if(upper_bounds[folan][filan] < upper_bounds[folan][smallest_ub[folan]]) smallest_ub[folan] = filan;
            }
            // std::cout << "doing pruning for point " << folan << " ..." << std::endl;
            // TODO: do step-wise pruning here: DONE
            // also fill out the ham_lb: DONE and elkan_lb: DONE somewhere
            // for ham_lb
            // fake_smallest_lb = lower_bounds[folan][smallest_ub[folan]];
            // fake_smallest_lb = DBL_MAX;
            // everybody is guilty until proven innocent...
            candidates_exist = false;
            if (level == int(log2(int(sqrt(D))) + 1)) {
                labels[folan] = smallest_ub[folan];
                // V9
                hamerly_upper_bounds[folan] = upper_bounds[folan][smallest_ub[folan]];
                assigned[folan] = 1;
                // v9
                // but we still need to update the ham_lb
                for (filan = 0; filan < K; filan++) {
                    if (filan == smallest_ub[folan]) continue;
                    if (!is_candidate[folan][filan]) continue;
                    if(!candidates_exist){
                        candidates_exist = true;
                        // fake_smallest_lb = DBL_MAX;
                        hamerly_lower_bounds[folan] = DBL_MAX;
                    }
                    // if(folan == 0){
                    //     std::cout << "before: HERE: filan: " << filan << " lower_bounds[folan][filan] " << lower_bounds[folan][filan] << " = " << elkan_lower_bounds[folan][filan] << std::endl;
                    //     std::cout << "HERE: fake_smallest_ub: " << fake_smallest_lb << " ham_lb[folan] " << hamerly_lower_bounds[folan] << std::endl;
                    //     std::cout << "----------------" << std::endl;
                    // }
                    
                    // find sec smallest lb
                    // if(lower_bounds[folan][filan] < fake_smallest_lb){
                    //     // keep sec smallest
                    //     hamerly_lower_bounds[folan] = fake_smallest_lb;
                    //     fake_smallest_lb = lower_bounds[folan][filan];
                    // } else if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                    //     hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                    // }

                    // TODO: so im changing the definition of ham_lb
                    // let ham_lb be the smallest lb betweent the last batch of candidates
                    // if there are other candidates at THIS point, then the ham_lb should just be fake_smallest_lb here
                    // but what if there are no other candidates?
                    // then we would have set fake_smallest_lb to DBL_MAX and it stays that way...
                    // I fixed it with the exists_candidate
                    // it's not the label == smallest_ub
                    // so it's the smallest lb to get pruned
                    if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                        hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                    }


                    // if(folan == 0){
                    //     std::cout << "after: HERE: filan: " << filan << " lower_bounds[folan][filan] " << lower_bounds[folan][filan] << " = " << elkan_lower_bounds[folan][filan] << std::endl;
                    //     std::cout << "HERE: fake_smallest_ub: " << fake_smallest_lb << " ham_lb[folan] " << hamerly_lower_bounds[folan] << std::endl;
                    //     std::cout << "----------------" << std::endl;
                    // }

                }
            } else {
                for (filan = 0; filan < K; filan++) {
                    if (filan == smallest_ub[folan]) continue;
                    if (!is_candidate[folan][filan]) continue;
                    // if(!candidates_exist){
                    //     candidates_exist = true;
                    //     // fake_smallest_lb = DBL_MAX;
                    //     hamerly_lower_bounds[folan] = DBL_MAX;
                    // }
                    if (lower_bounds[folan][filan] >= upper_bounds[folan][smallest_ub[folan]]) {
                        is_candidate[folan][filan] = false;
                        if(!candidates_exist){
                            candidates_exist = true;
                            // fake_smallest_lb = DBL_MAX;
                            hamerly_lower_bounds[folan] = DBL_MAX;
                        }
                        // V9
                        // we should also try to find the second closest centroid
                        // TODO: im not sure if i should clear the old ham_lbs from last round
                        // what if it's masked completely, then we would never get here, so the old one should still be valid.
                        // nvm, i think im right, i won't clear it
                        // if (lower_bounds[folan][filan] < hamerly_lower_bounds[folan]) {
                        //     hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                        // }
                        // V9
                        // I'm fairly sure the above is bullshit, can't remember what I was thinking
                        // so ill make another
                        if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                            hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                        }
                        
                        
                    }
                    /*if(folan == 0){
                        std::cout << "before: HERE: filan: " << filan << " lower_bounds[folan][filan] " << lower_bounds[folan][filan] << " = " << elkan_lower_bounds[folan][filan] << std::endl;
                        std::cout << "HERE: ham_lb[folan] " << hamerly_lower_bounds[folan] << std::endl;
                        std::cout << "----------------" << std::endl;
                    }*/
                    
                    // find sec smallest lb
                    // if(lower_bounds[folan][filan] < fake_smallest_lb){
                    //     // keep sec smallest
                    //     hamerly_lower_bounds[folan] = fake_smallest_lb;
                    //     fake_smallest_lb = lower_bounds[folan][filan];
                    // } else if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                    //     hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                    // }
                    
                    // change of ham_lb definition

                    // if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                    //     hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                    // }

                    /*if(folan == 0){
                        std::cout << "after: HERE: filan: " << filan << " lower_bounds[folan][filan] " << lower_bounds[folan][filan] << " = " << elkan_lower_bounds[folan][filan] << std::endl;
                        std::cout << "HERE: ham_lb[folan] " << hamerly_lower_bounds[folan] << std::endl;
                        std::cout << "----------------" << std::endl;
                    }*/
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

    //std::cout << "set labels " << std::endl;
        

    // sanity check: is everyone assigned?
    for(folan = 0; folan < N; folan++){
        if(assigned[folan] == 0){
            std::cout << folan << " NOT ASSIGNED" << std::endl;
            exit(2);
        }
    }

    //std::cout << "after all-assigned sanity check" << std::endl;

    // calc new centroids
    // make copy of old centroids
    // apparently this also does not work, copies pointer somehow I think...
    // memcpy(old_centroids, centroids, sizeof(double) * K * D);
    for (folan = 0; folan < K; folan++) {
        for (filan = 0; filan < D; filan++) {
            old_centroids[folan][filan] = centroids[folan][filan];
        }
    }
    //std::cout << "copied centroids to old centroids" << std::endl;
    // set centroids to 0
    memset(cluster_counts, 0, sizeof(int) * K);
    //std::cout << "set cluster counts to 0" << std::endl;
    for (folan = 0; folan < K; folan++) {
        // just testing
        // cluster_counts[folan] = 0;
        for (filan = 0; filan < D; filan++) {
            centroids[folan][filan] = 0.0;
        }
    }
    // This doesn't work on doubles
    // memset(centroids, 0, sizeof(double) * K * D);
    //std::cout << "after all the memcpys" << std::endl;

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
    //std::cout << "calculated new centroids" << std::endl;
    calculate_centroids_square_sums();
    //std::cout << "ssed them " << std::endl;
    // just to check
    int sanity_check = 0;
    //std::cout << "cluster counts..." << std::endl;
    for (folan = 0; folan < K; folan++) {
        sanity_check += cluster_counts[folan];
        //std::cout << cluster_counts[folan] << " ";
    }
    //std::cout << sanity_check << std::endl;

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
    //std::cout << "calculated centroid movements" << std::endl;

    
    /*std::cout << "ham_ub[1] " << hamerly_upper_bounds[1] << " ham_lb[1] " << hamerly_lower_bounds[1] <<  std::endl;
    for(filan =0; filan < K; filan++){
        std::cout << elkan_lower_bounds[1][filan] << " ";
    }
    std::cout << std::endl;*/

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

    /*std::cout << "after updating the bounds...\n";

    std::cout << "ham_ub[1] " << hamerly_upper_bounds[1] << " ham_lb[1] " << hamerly_lower_bounds[1] <<  std::endl;
    for(filan =0; filan < K; filan++){
        std::cout << elkan_lower_bounds[1][filan] << " ";
    }
    std::cout << std::endl;*/




    

    // loop over max_iter
    for (int iter = 1; iter < MAX_ITERATIONS; iter++) {
        level = 0;
        //std::cout << "iteration " << iter << "..." << std::endl;

        // calculate closest centroid to each centroid
        for (folan = 0; folan < K; folan++) {
            smallest = DBL_MAX;
            for (filan = folan; filan < K; filan++) {
                tmp = centroids_ss[filan][0] + centroids_ss[folan][0];
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
        //std::cout << "found closest distance to each centroid, and the centr-centr distances" << std::endl;


        // fill out mask
        // set all to 1 first
        for (folan = 0; folan < N; folan++) {
            for (filan = 0; filan < K; filan++) {
                is_candidate[folan][filan] = true;
            }
        }
        //std::cout << "filled out is_candid" << std::endl;


        // set the assigned
        memset(assigned, 0, sizeof(int) * N);
        //std::cout << "memset assigned" << std::endl;

        for(folan = 0; folan < N; folan++){
            for(filan=0; filan <K; filan++){
                last_level_calculated[folan][filan] = -1;
            }
        }
        //std::cout << "set last_level_calc to -1" << std::endl;

        // I think I have to put the step-wise while here

        // TODO:
        // for x in X:
        //  if ham_bound[x] is true: set is_candid to false for all centroid != labels[x] and set assigned[x] = 1
        for (folan = 0; folan < N; folan++) {
            // hamerly check
            hamerly_bound = ((0.5 * closest_centroid_distance[labels[folan]]) > hamerly_lower_bounds[folan]) ? (
                    0.5 * closest_centroid_distance[labels[folan]]) : hamerly_lower_bounds[folan];
            if (hamerly_bound >= hamerly_upper_bounds[folan]) {
                hamerly_count += K;
                // this one will not move, so none of them are candidates, it is assigned...
                assigned[folan] = 1;
            }
        }

        while (level < int(log2(int(sqrt(D))) + 2)) {
            //std::cout << "in while level " << level << std::endl;
            // set smallest ub
            // let's initiate them to the previous labels so skipping the labels in the for K does not change anything
            // TODO!!!!
            // memset(smallest_ub, 0, sizeof(int) * N);


            for (folan = 0; folan < N; folan++) {
                // TODO: check if sth goes wrong
                //smallest_ub[folan] = labels[folan];
                if (assigned[folan] > 0) continue;
                r = true;
                r_int = 0; fake_smallest_lb = DBL_MAX;
                // hamerly_bound = ((0.5 * closest_centroid_distance[labels[folan]]) > hamerly_lower_bounds[folan]) ? (
                //         0.5 * closest_centroid_distance[labels[folan]]) : hamerly_lower_bounds[folan];
                // elkan lemma 1 and hamerly
                // if (hamerly_bound < hamerly_upper_bounds[folan]) {
                
                
                // if(true){
                /*if(folan == 4){
                    std::cout << "not pruned by ham_bound..." << std::endl;
                }*/
                // if((0.5 * closest_centroid_distance[labels[folan]]) < hamerly_upper_bounds[folan]){
                
                // update the upper bound
                // im gonna do it till this level only
                calculate_distance_folan_filan_till_level(folan, labels[folan], level);
                // update the ham_ub and elkan_lb if better
                // TODO: not sure, maybe we have to update it anyway
                // hamerly_upper_bounds[folan] = std::min(hamerly_upper_bounds[folan], upper_bounds[folan][labels[folan]]);
                // elkan_lower_bounds[folan][labels[folan]] = std::max(elkan_lower_bounds[folan][labels[folan]], lower_bounds[folan][labels[folan]]);
                // AS DISCUSSED WITH KASPER, lets's try it this way for now TODO: check
                hamerly_upper_bounds[folan] = upper_bounds[folan][labels[folan]];
                elkan_lower_bounds[folan][labels[folan]] = lower_bounds[folan][labels[folan]];
                /*if(folan == 4){
                    std::cout << "JUST FOR THE SAKE OF FATEMEH'S SANITY: euc_dist(1, 1) " << sqrt(euclidean_distance(1, 1)) << std::endl; 
                    std::cout << "1: changed ham_ub[" << folan << "] to " << hamerly_upper_bounds[folan];
                    std::cout << " and elkan_lb[" << folan << "][" << labels[folan] << "] to " << elkan_lower_bounds[folan][labels[folan]] << std::endl;
                }*/
                
                for (filan = 0; filan < K; filan++) {
                    // TODO: tentative...
                    if(!is_candidate[folan][filan]) continue;
                    // if(folan == 1){
                    //     std::cout << "am i going crazy? " << filan << " " << labels[folan] << std::endl;
                    // }
                    if (filan == labels[folan]) continue;
                    // new lemma
                    
                    /*if(folan == 4){
                        std::cout << "not pruned by new lemma " << filan << std::endl;
                    }*/
                    if (hamerly_upper_bounds[folan] > elkan_lower_bounds[folan][filan]
                        && hamerly_upper_bounds[folan] >
                            (0.5 * centroid_to_centroid_distances[filan][labels[folan]])) {
                            // HERE
                            // I'm thinking the full dist calculation of the prev lables[folan]
                            // we have to do it all in one go
                            // or the r bool has to be activated for all the levels or sth
                            
                            // if(r_int < int(log2(int(sqrt(D))) + 2)){
                            
                            // TODO: maybe I should have a if(is_candid) around this
                            // we do the TI pruning with the ifs, then the step-wise prunings with the is_candid
                            // or I could just have a if(!is_candid) continue at the start of for K ?
                            // I think this could work
                            calculate_distance_folan_filan_till_level(folan, filan, level);
                            // i'm not sure about this maxing, maybe I should always set it to lb TODO?
                            elkan_lower_bounds[folan][filan] = lower_bounds[folan][filan];
                            /*if(folan == 4){
                                std::cout << "2: changed  elkan_lb[" << folan << "][" << filan << "] to " << elkan_lower_bounds[folan][filan] << " ub was " << upper_bounds[folan][filan] << std::endl;
                            }*/    
                            // elkan_lower_bounds[folan][filan] = std::max(elkan_lower_bounds[folan][filan], lower_bounds[folan][filan]);

                            // then (in non-step-wise) we check if the dist is smaller than the one for labels[folan]
                            // I'm not sure we can do that since we don't exactly have distances, but ub and lb
                            // we also don't have labels[folan], in the same way as before
                            // we could store a smallest_ub_till_now
                            // then use that as the labels[folan]
                            // we could even set the labels[folan] = smallest_ub_till_now and update as we go
                            // then we'd do the ic pruning here to avoid an extra loop
                            // no, that wouldn't really work. since the smallest_ub_till_now is not necessarily the smallest_ub over all
                            // we don't need a separate variable for smallest_ub, it could just be labels[folan]
                            // then the problem is when do we say it is assigned.
                            // TODO
                            // maybe I should update ham_lb separately
                            // like, if they are not from the same point, it would still be correct
                            // I think...
                            // yeah but then ham_ub==labels[folan] would be the smallest_ub
                            // but ham_lb is supposed to be the second smallest lb:-?
                            // or the second largest lb;-? 
                            // I don't think it should be sec largest lb
                            // i'll go with sec smallest lb, fingers crossed --> for this I will have to keep smallest_lb too
                            if(upper_bounds[folan][filan] < hamerly_upper_bounds[folan]){
                                // TODO: how the f am i supposed to update ham_lb now?!
                                // update the ham_lb if necessary
                                // hamerly_lower_bounds[folan] = std::max(hamerly_lower_bounds[folan], hamerly_upper_bounds[folan]);
                                // hamerly_lower_bounds[folan] = hamerly_upper_bounds[folan];
                                labels[folan] = filan;
                                hamerly_upper_bounds[folan] = upper_bounds[folan][filan];
                                /*if(folan == 4){
                                    std::cout << "3: changed  ham_ub[" << folan << "] to " << hamerly_upper_bounds[folan];
                                    std::cout << " and labels[" << folan << "]  to " << labels[folan] << std::endl;
                                }*/
                            }
                            // I have chnaged the def of ham_lb
                            // I will only update the ham_lb in the pruning of stepwise
                            // I only use the ham_lb in the ham_bound if, so I can just update it after the whole for K
                            // will it ever stay unassigned?
                            // after this for K
                            // we go into the step-wise pruning
                            // then, if there are candidates, the smallest lb of the candidates(other than the smallest ub) will be ham_lb
                            // then if there are no candidates, the ham_lb of last level remains
                            // the first level will always have candidates(everybody is a candidate)
                            // so we will always assign it
                            // no im doubting if the def is correct
                            // since the "the smallest lb of the candidates(other than the smallest ub)"
                            // does not really add up to much
                            // idk, let's see what happens:) 
                            // finding sec smallest lb
                            // if(lower_bounds[folan][filan] < fake_smallest_lb){
                            //     // TODO: maybe the max completely ruins the whole thing...
                            //     // i'm not gonna do the max for now
                            //     // hamerly_lower_bounds[folan] = std::max(hamerly_lower_bounds[folan], fake_smallest_lb);
                            //     hamerly_lower_bounds[folan] = fake_smallest_lb;
                            //     fake_smallest_lb = lower_bounds[folan][filan];
                            // } else if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                            //     hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                            // }
                    } 
                    // FOR NOW I WON'T
                    /*else{
                        // TODO: I think I need to change the is_candid
                        if(folan == 4){
                            std::cout << "3.5 skipping this " << filan << std::endl;
                        }
                    }*/

                    // im gonna update smallest_ub here
                    // I think it should be fine
                    // cause after all the hassle, the ub and lb here should be valid
                    //if(upper_bounds[folan][filan] < upper_bounds[folan][smallest_ub[folan]]){ smallest_ub[folan] = filan; }
                } // for filan < K
                // } // if ham_bound
                
                // TODO: now we have to do the pruning based off step-wise... I think...
                // I'll move this to the start of the N loop
                // I added it there, but will keep it here because of the else that can change assigned with ham_boung
                if (assigned[folan] > 0) continue;
                /*if(folan == 4){
                    std::cout << "elkan_lb[4] ";
                    for(filan = 0; filan < K; filan++){
                        std::cout << elkan_lower_bounds[folan][filan] << " ";
                    }
                    std::cout << std::endl;
                }*/

                candidates_exist = false;
                // V9
                if (level == int(log2(int(sqrt(D))) + 1)) {
                    //labels[folan] = smallest_ub[folan];
                    // V9
                    //hamerly_upper_bounds[folan] = upper_bounds[folan][smallest_ub[folan]];
                    assigned[folan] = 1;
                    // v9

                    // but we still need to update the ham_lb
                    for (filan = 0; filan < K; filan++) {
                        // if (filan == smallest_ub[folan]) continue;
                        if (filan == labels[folan]) continue;
                        if (!is_candidate[folan][filan]) continue;
                        if(!candidates_exist){
                            candidates_exist = true;
                            // fake_smallest_lb = DBL_MAX;
                            hamerly_lower_bounds[folan] = DBL_MAX;
                        }

                        // TODO: so im changing the definition of ham_lb
                        // let ham_lb be the smallest lb betweent the last batch of candidates
                        // if there are other candidates at THIS point, then the ham_lb should just be fake_smallest_lb here
                        // but what if there are no other candidates?
                        // then we would have set fake_smallest_lb to DBL_MAX and it stays that way...
                        // I fixed it with the exists_candidate
                        // it's not the label == smallest_ub
                        // so it's the smallest lb to get pruned
                        if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                            hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                        }

                    }

                } else {
                    for (filan = 0; filan < K; filan++) {
                        // if (filan == smallest_ub[folan]) continue;
                        if (filan == labels[folan]) continue;
                        if (!is_candidate[folan][filan]) continue;
                        // if(!candidates_exist){
                        //     candidates_exist = true;
                        //     // fake_smallest_lb = DBL_MAX;
                        //     hamerly_lower_bounds[folan] = DBL_MAX;
                        // }
                        // if (lower_bounds[folan][filan] >= upper_bounds[folan][smallest_ub[folan]]) {
                        if (lower_bounds[folan][filan] >= upper_bounds[folan][labels[folan]]) {
                            is_candidate[folan][filan] = false;
                            if(!candidates_exist){
                                candidates_exist = true;
                                // fake_smallest_lb = DBL_MAX;
                                hamerly_lower_bounds[folan] = DBL_MAX;
                            }
                            /*if(folan == 4){
                                std::cout << "4: pruning " << filan << std::endl;
                            }*/
                            // NOT HERE i THINK
                            // V9
                            // we should also try to find the second closest centroid
                            // TODO: im not sure if i should clear the old ham_lbs from last round
                            // what if it's masked completely, then we would never get here, so the old one should still be valid.
                            // nvm, i think im right, i won't clear it
                            // if (lower_bounds[folan][filan] < hamerly_lower_bounds[folan]) {
                            //     hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                            // }
                            // V9
                            if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                                hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                                /*if(folan == 4){
                                    std::cout << "5: chnaged ham_lb[" << folan << "] to " << hamerly_lower_bounds[folan] << std::endl; 
                                }*/
                            }
                        }
                        
                        // if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                        //     hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                        //     if(folan == 4){
                        //         std::cout << "5: chnaged ham_lb[" << folan << "] to " << hamerly_lower_bounds[folan] << std::endl; 
                        //     }
                        // }
                    }

                    // check if only one is true
                    alaki = 0;
                    for (filan = 0; filan < K; filan++) {
                        if (is_candidate[folan][filan]) alaki++;
                    }
                    // then the only one left is the one with the smallest_ub
                    if (alaki == 1) {
                        //labels[folan] = smallest_ub[folan];
                        // V9
                        //hamerly_upper_bounds[folan] = upper_bounds[folan][smallest_ub[folan]];
                        assigned[folan] = 1;
                        // v9
                        /*if(folan == 4){
                            std::cout << "6: only one candidate chnaged labels[" << folan << "] to " << labels[folan];
                            std::cout << " and ham_ub[" << folan << "] to " << hamerly_upper_bounds[folan] << std::endl; 
                        }*/
                    }
                }
                
            } // for folan < N
            level++;
        } // step-wise loop
        
        // TODO: sanity check, is everybody assigned?
        /*for(folan = 0; folan < N; folan++){
            if(assigned[folan] == 0){
                std::cout << folan << " NOT ASSIGNED" << std::endl;
                exit(2);
            }
        }*/

        //std::cout << "set labels" << std::endl;


        // calc new centroids
        // make copy of old centroids
        // apparently this also does not work, copies pointer somehow I think...
        // memcpy(old_centroids, centroids, sizeof(double) * K * D);
        for (folan = 0; folan < K; folan++) {
            for (filan = 0; filan < D; filan++) {
                old_centroids[folan][filan] = centroids[folan][filan];
            }
        }
        //std::cout << "copied centroids to old centroids" << std::endl;
        // set centroids to 0
        memset(cluster_counts, 0, sizeof(int) * K);
        //std::cout << "set cluster counts to 0" << std::endl;
        for (folan = 0; folan < K; folan++) {
            // just testing
            // cluster_counts[folan] = 0;
            for (filan = 0; filan < D; filan++) {
                centroids[folan][filan] = 0.0;
            }
        }
        // This doesn't work on doubles
        // memset(centroids, 0, sizeof(double) * K * D);
        //std::cout << "after all the memcpys" << std::endl;

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
        //std::cout << "calculated new centroids" << std::endl;
        calculate_centroids_square_sums();
        //td::cout << "ssed them " << std::endl;
        // just to check
        /*int sanity_check = 0;
        std::cout << "cluster counts..." << std::endl;
        for (folan = 0; folan < K; folan++) {
            sanity_check += cluster_counts[folan];
            std::cout << cluster_counts[folan] << " ";
        }
        std::cout << sanity_check << std::endl;*/

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
        //std::cout << "calculated centroid movements" << std::endl;

        /*std::cout << "ham_ub[1] " << hamerly_upper_bounds[1] << " ham_lb[1] " << hamerly_lower_bounds[1] << std::endl;
        for(filan = 0; filan < K; filan++){
            std::cout << elkan_lower_bounds[1][filan] << " " ;
        }
        std::cout << std::endl;*/


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
        //std::cout << "checked convergence" << std::endl;

        // end if converged
        if (has_converged) {
            if (D >= 256*256) {
                int sanity_check = 0;
                for (folan = 0; folan < K; folan++) {
                    sanity_check += cluster_counts[folan];
                    std::cout << cluster_counts[folan] << " ";
                }
                std::cout << sanity_check << std::endl;
            }
            break;
        }

    } // for iter

}


// let's try for the whole thing
// hamerly + elkan + step-wise
// copied from v105
// simpler ham_lb calc
// TODO
void kmeans_v106(){

    int folan = 0;
    int filan = 0; int alaki;
    int ashghal = 0;// for loop usage
    int furthest_moving_centroid, second_furthest_moving_centroid;
    double smallest, second_smallest;
    double tmp, hamerly_bound;
    int hamerly_count = 0;
    bool r; int level;
    int r_int; double fake_smallest_lb; bool candidates_exist;
    int *smallest_ub = (int *) malloc(N * sizeof(int));

    // set initial centroids

    for(folan = 0; folan < K; folan++){
        for(filan = 0; filan < D; filan++){
            centroids[folan][filan] = data_arr[folan][filan];
        }
    }

    // memcpy(centroids, data_arr, sizeof(double) * K * D);
    calculate_centroids_square_sums();
    //std::cout << "copied init centroids" << std::endl;


    bool has_converged = false;
    int cluster_counts[K];


    // let's do the first iteration out here, so the if inside the loop is resolved
    // TODO
    level = 0;
    // calculate closest centroid to each centroid
    for (folan = 0; folan < K; folan++) {
        smallest = DBL_MAX;
        for (filan = folan; filan < K; filan++) {
            tmp = centroids_ss[filan][0] + centroids_ss[folan][0];
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
    //std::cout << "found closest distance to each centroid, and the centr-centr distances" << std::endl;

    // fill out mask
    // set all to 1 first
    for (folan = 0; folan < N; folan++) {
        for (filan = 0; filan < K; filan++) {
            is_candidate[folan][filan] = true;
        }
    }

    //std::cout << "filled out is_candid" << std::endl;

    // set the assigned
    memset(assigned, 0, sizeof(int) * N);
    memset(smallest_ub, 0, sizeof(int) * N);

    //std::cout << "memset out assigned and smallest_ub" << std::endl;

    for(folan = 0; folan < N; folan++){
        for(filan=0; filan <K; filan++){
            last_level_calculated[folan][filan] = -1;
        }
    }
    //std::cout << "set last_level_calc to -1" << std::endl;

    // I think I have to put the step-wise while here

    // TODO:
    // for x in X:
    //  check ham_bound: if true, then make every entroid other than labels[x] not candidate


    while (level < int(log2(int(sqrt(D))) + 2)) {
        //std::cout << "iter 0, in while level =  " << level << std::endl;
        for(folan = 0; folan < N; folan++){
            if(assigned[folan] > 0) continue;
            for(filan = 0; filan < K; filan++){
                if(!is_candidate[folan][filan]) continue;
                calculate_distance_folan_filan_till_level(folan, filan, level);
                // fill elkan lb
                elkan_lower_bounds[folan][filan] = lower_bounds[folan][filan];
                // find smallest_ub
                if(upper_bounds[folan][filan] < upper_bounds[folan][smallest_ub[folan]]) smallest_ub[folan] = filan;
            }
            // std::cout << "doing pruning for point " << folan << " ..." << std::endl;
            // TODO: do step-wise pruning here: DONE
            // also fill out the ham_lb: DONE and elkan_lb: DONE somewhere
            // for ham_lb
            // fake_smallest_lb = lower_bounds[folan][smallest_ub[folan]];
            // fake_smallest_lb = DBL_MAX;
            // everybody is guilty until proven innocent...
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

    // find ham_lb
    for(folan = 0; folan < N; folan++){
        hamerly_lower_bounds[folan] = DBL_MAX;
        for(filan = 0; filan < K; filan++){
            if(labels[folan] != filan && elkan_lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                hamerly_lower_bounds[folan] = elkan_lower_bounds[folan][filan];
            }
        }
    }


    //std::cout << "set labels " << std::endl;


    // sanity check: is everyone assigned?
    /*for(folan = 0; folan < N; folan++){
        if(assigned[folan] == 0){
            std::cout << folan << " NOT ASSIGNED" << std::endl;
            exit(2);
        }
    }*/

    //std::cout << "after all-assigned sanity check" << std::endl;

    // calc new centroids
    // make copy of old centroids
    // apparently this also does not work, copies pointer somehow I think...
    // memcpy(old_centroids, centroids, sizeof(double) * K * D);
    for (folan = 0; folan < K; folan++) {
        for (filan = 0; filan < D; filan++) {
            old_centroids[folan][filan] = centroids[folan][filan];
        }
    }
    //std::cout << "copied centroids to old centroids" << std::endl;
    // set centroids to 0
    memset(cluster_counts, 0, sizeof(int) * K);
    //std::cout << "set cluster counts to 0" << std::endl;
    for (folan = 0; folan < K; folan++) {
        // just testing
        // cluster_counts[folan] = 0;
        for (filan = 0; filan < D; filan++) {
            centroids[folan][filan] = 0.0;
        }
    }
    // This doesn't work on doubles
    // memset(centroids, 0, sizeof(double) * K * D);
    //std::cout << "after all the memcpys" << std::endl;

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
    //std::cout << "calculated new centroids" << std::endl;
    calculate_centroids_square_sums();
    //std::cout << "ssed them " << std::endl;
    // just to check
    /*int sanity_check = 0;
    std::cout << "cluster counts..." << std::endl;
    for (folan = 0; folan < K; folan++) {
        sanity_check += cluster_counts[folan];
        std::cout << cluster_counts[folan] << " ";
    }
    std::cout << sanity_check << std::endl;*/

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
    //std::cout << "calculated centroid movements" << std::endl;


    /*std::cout << "ham_ub[1] " << hamerly_upper_bounds[1] << " ham_lb[1] " << hamerly_lower_bounds[1] <<  std::endl;
    for(filan =0; filan < K; filan++){
        std::cout << elkan_lower_bounds[1][filan] << " ";
    }
    std::cout << std::endl;*/

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

    /*std::cout << "after updating the bounds...\n";

    std::cout << "ham_ub[1] " << hamerly_upper_bounds[1] << " ham_lb[1] " << hamerly_lower_bounds[1] <<  std::endl;
    for(filan =0; filan < K; filan++){
        std::cout << elkan_lower_bounds[1][filan] << " ";
    }
    std::cout << std::endl;*/






    // loop over max_iter
    for (int iter = 1; iter < MAX_ITERATIONS; iter++) {
        level = 0;
        //std::cout << "iteration " << iter << "..." << std::endl;

        // calculate closest centroid to each centroid
        for (folan = 0; folan < K; folan++) {
            smallest = DBL_MAX;
            for (filan = folan; filan < K; filan++) {
                tmp = centroids_ss[filan][0] + centroids_ss[folan][0];
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
        //std::cout << "found closest distance to each centroid, and the centr-centr distances" << std::endl;


        // fill out mask
        // set all to 1 first
        for (folan = 0; folan < N; folan++) {
            for (filan = 0; filan < K; filan++) {
                is_candidate[folan][filan] = true;
            }
        }
        //std::cout << "filled out is_candid" << std::endl;


        // set the assigned
        memset(assigned, 0, sizeof(int) * N);
        //std::cout << "memset assigned" << std::endl;

        for(folan = 0; folan < N; folan++){
            for(filan=0; filan <K; filan++){
                last_level_calculated[folan][filan] = -1;
            }
        }
        //std::cout << "set last_level_calc to -1" << std::endl;

        // I think I have to put the step-wise while here

        // TODO:
        // for x in X:
        //  if ham_bound[x] is true: set is_candid to false for all centroid != labels[x] and set assigned[x] = 1
        for (folan = 0; folan < N; folan++) {
            // hamerly check
            hamerly_bound = ((0.5 * closest_centroid_distance[labels[folan]]) > hamerly_lower_bounds[folan]) ? (
                    0.5 * closest_centroid_distance[labels[folan]]) : hamerly_lower_bounds[folan];
            if (hamerly_bound >= hamerly_upper_bounds[folan]) {
                hamerly_count += K;
                // this one will not move, so none of them are candidates, it is assigned...
                assigned[folan] = 1;
            }
        }

        while (level < int(log2(int(sqrt(D))) + 2)) {
            //std::cout << "in while level " << level << std::endl;
            // set smallest ub
            // let's initiate them to the previous labels so skipping the labels in the for K does not change anything
            // TODO!!!!
            // memset(smallest_ub, 0, sizeof(int) * N);


            for (folan = 0; folan < N; folan++) {
                // TODO: check if sth goes wrong
                smallest_ub[folan] = labels[folan];
                if (assigned[folan] > 0) continue;

                calculate_distance_folan_filan_till_level(folan, labels[folan], level);

                hamerly_upper_bounds[folan] = upper_bounds[folan][labels[folan]];
                elkan_lower_bounds[folan][labels[folan]] = lower_bounds[folan][labels[folan]];


                for (filan = 0; filan < K; filan++) {
                    if(!is_candidate[folan][filan]) continue;

                    if (filan == labels[folan]) continue;

                    if (hamerly_upper_bounds[folan] > elkan_lower_bounds[folan][filan]
                        && hamerly_upper_bounds[folan] >
                            (0.5 * centroid_to_centroid_distances[filan][labels[folan]])) {

                            calculate_distance_folan_filan_till_level(folan, filan, level);
                            // i'm not sure about this maxing, maybe I should always set it to lb TODO?
                            elkan_lower_bounds[folan][filan] = lower_bounds[folan][filan];

                            if(upper_bounds[folan][filan] < hamerly_upper_bounds[folan]){
                                labels[folan] = filan;
                                hamerly_upper_bounds[folan] = upper_bounds[folan][filan];

                            }

                    } 
                    if(upper_bounds[folan][filan] < upper_bounds[folan][smallest_ub[folan]]){ smallest_ub[folan] = filan; }
                } // for filan < K

                // TODO: now we have to do the pruning based off step-wise... I think...
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

            } // for folan < N
            level++;
        } // step-wise loop

        // find ham_lb
        for(folan = 0; folan < N; folan++){
            hamerly_lower_bounds[folan] = DBL_MAX;
            for(filan = 0; filan < K; filan++){
                if(labels[folan] != filan && elkan_lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                    hamerly_lower_bounds[folan] = elkan_lower_bounds[folan][filan];
                }
            }
        }

        // TODO: sanity check, is everybody assigned?
        /*for(folan = 0; folan < N; folan++){
            if(assigned[folan] == 0){
                std::cout << folan << " NOT ASSIGNED" << std::endl;
                exit(2);
            }
        }

        std::cout << "set labels" << std::endl;*/


        // calc new centroids
        // make copy of old centroids
        // apparently this also does not work, copies pointer somehow I think...
        // memcpy(old_centroids, centroids, sizeof(double) * K * D);
        for (folan = 0; folan < K; folan++) {
            for (filan = 0; filan < D; filan++) {
                old_centroids[folan][filan] = centroids[folan][filan];
            }
        }
        //std::cout << "copied centroids to old centroids" << std::endl;
        // set centroids to 0
        memset(cluster_counts, 0, sizeof(int) * K);
        //std::cout << "set cluster counts to 0" << std::endl;
        for (folan = 0; folan < K; folan++) {
            // just testing
            // cluster_counts[folan] = 0;
            for (filan = 0; filan < D; filan++) {
                centroids[folan][filan] = 0.0;
            }
        }
        // This doesn't work on doubles
        // memset(centroids, 0, sizeof(double) * K * D);
        //std::cout << "after all the memcpys" << std::endl;

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
        //std::cout << "calculated new centroids" << std::endl;
        calculate_centroids_square_sums();
        //std::cout << "ssed them " << std::endl;
        // just to check
        /*int sanity_check = 0;
        std::cout << "cluster counts..." << std::endl;
        for (folan = 0; folan < K; folan++) {
            sanity_check += cluster_counts[folan];
            std::cout << cluster_counts[folan] << " ";
        }
        std::cout << sanity_check << std::endl;*/

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
        /*std::cout << "calculated centroid movements" << std::endl;

        std::cout << "ham_ub[1] " << hamerly_upper_bounds[1] << " ham_lb[1] " << hamerly_lower_bounds[1] << std::endl;
        for(filan = 0; filan < K; filan++){
            std::cout << elkan_lower_bounds[1][filan] << " " ;
        }
        std::cout << std::endl;*/


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
        //std::cout << "checked convergence" << std::endl;

        // end if converged
        if (has_converged) {
            if (D >= 256*256) {
                int sanity_check = 0;
                for (folan = 0; folan < K; folan++) {
                    sanity_check += cluster_counts[folan];
                    std::cout << cluster_counts[folan] << " ";
                }
                std::cout << sanity_check << std::endl;
            }
            break;
        }

    } // for iter

}

// Kasper's idea
// hamerly + elkan + step-wise, but hybrid
// the point is that if the level gets high enough(eg more than 1024 dimensions as elkan suggests)
// then we give up on checking the TI, just do step-wise
// copied from v105
// TODO: DONE and tested
void kmeans_v11(){
    //std::cout << "in kmeans_v11 ..." << std::endl;
    int folan = 0;
    int filan = 0; int alaki;
    int ashghal = 0;// for loop usage
    int furthest_moving_centroid, second_furthest_moving_centroid;
    double smallest, second_smallest;
    double tmp, hamerly_bound;
    int hamerly_count = 0;
    bool r; int level;
    int r_int; double fake_smallest_lb; bool candidates_exist;
    int *smallest_ub = (int *) malloc(N * sizeof(int));

    // set initial centroids

    for(folan = 0; folan < K; folan++){
        for(filan = 0; filan < D; filan++){
            centroids[folan][filan] = data_arr[folan][filan];
        }
    }

    // memcpy(centroids, data_arr, sizeof(double) * K * D);
    calculate_centroids_square_sums();
    //std::cout << "copied init centroids" << std::endl;


    bool has_converged = false;
    int cluster_counts[K];


    // let's do the first iteration out here, so the if inside the loop is resolved
    // TODO
    level = 0;
    // calculate closest centroid to each centroid
    for (folan = 0; folan < K; folan++) {
        smallest = DBL_MAX;
        for (filan = folan; filan < K; filan++) {
            tmp = centroids_ss[filan][0] + centroids_ss[folan][0];
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
    //std::cout << "found closest distance to each centroid, and the centr-centr distances" << std::endl;



    // fill out mask
    // set all to 1 first
    for (folan = 0; folan < N; folan++) {
        for (filan = 0; filan < K; filan++) {
            is_candidate[folan][filan] = true;
        }
    }

    //std::cout << "filled out is_candid" << std::endl;

    // set the assigned
    memset(assigned, 0, sizeof(int) * N);
    memset(smallest_ub, 0, sizeof(int) * N);

    //std::cout << "memset out assigned and smallest_ub" << std::endl;

    for(folan = 0; folan < N; folan++){
        for(filan=0; filan <K; filan++){
            last_level_calculated[folan][filan] = -1;
        }
    }
    //std::cout << "set last_level_calc to -1" << std::endl;

    // I think I have to put the step-wise while here
    
    // TODO: I should also do the switch to simple stepwise in iter 0, but for now, I will start with the iterations in the loop
    // DONE we don't do any checks anyway so no need for the if

    while (level < int(log2(int(sqrt(D))) + 2)) {
        //std::cout << "iter 0, in while level =  " << level << std::endl;
        for(folan = 0; folan < N; folan++){
            if(assigned[folan] > 0) continue;
            for(filan = 0; filan < K; filan++){
                if(!is_candidate[folan][filan]) continue;
                calculate_distance_folan_filan_till_level(folan, filan, level);
                // fill elkan lb
                elkan_lower_bounds[folan][filan] = lower_bounds[folan][filan];
                // find smallest_ub
                if(upper_bounds[folan][filan] < upper_bounds[folan][smallest_ub[folan]]) smallest_ub[folan] = filan;
            }
            // std::cout << "doing pruning for point " << folan << " ..." << std::endl;
            // TODO: do step-wise pruning here: DONE
            // also fill out the ham_lb: DONE and elkan_lb: DONE somewhere
            // for ham_lb
            // fake_smallest_lb = lower_bounds[folan][smallest_ub[folan]];
            // fake_smallest_lb = DBL_MAX;
            // everybody is guilty until proven innocent...
            candidates_exist = false;
            if (level == int(log2(int(sqrt(D))) + 1)) {
                labels[folan] = smallest_ub[folan];
                // V9
                hamerly_upper_bounds[folan] = upper_bounds[folan][smallest_ub[folan]];
                assigned[folan] = 1;
                // v9
                // but we still need to update the ham_lb
                for (filan = 0; filan < K; filan++) {
                    if (filan == smallest_ub[folan]) continue;
                    if (!is_candidate[folan][filan]) continue;
                    if(!candidates_exist){
                        candidates_exist = true;
                        // fake_smallest_lb = DBL_MAX;
                        hamerly_lower_bounds[folan] = DBL_MAX;
                    }
                    // if(folan == 0){
                    //     std::cout << "before: HERE: filan: " << filan << " lower_bounds[folan][filan] " << lower_bounds[folan][filan] << " = " << elkan_lower_bounds[folan][filan] << std::endl;
                    //     std::cout << "HERE: fake_smallest_ub: " << fake_smallest_lb << " ham_lb[folan] " << hamerly_lower_bounds[folan] << std::endl;
                    //     std::cout << "----------------" << std::endl;
                    // }
                    
                    // find sec smallest lb
                    // if(lower_bounds[folan][filan] < fake_smallest_lb){
                    //     // keep sec smallest
                    //     hamerly_lower_bounds[folan] = fake_smallest_lb;
                    //     fake_smallest_lb = lower_bounds[folan][filan];
                    // } else if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                    //     hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                    // }

                    // TODO: so im changing the definition of ham_lb
                    // let ham_lb be the smallest lb betweent the last batch of candidates
                    // if there are other candidates at THIS point, then the ham_lb should just be fake_smallest_lb here
                    // but what if there are no other candidates?
                    // then we would have set fake_smallest_lb to DBL_MAX and it stays that way...
                    // I fixed it with the exists_candidate
                    // it's not the label == smallest_ub
                    // so it's the smallest lb to get pruned
                    if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                        hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                    }


                    // if(folan == 0){
                    //     std::cout << "after: HERE: filan: " << filan << " lower_bounds[folan][filan] " << lower_bounds[folan][filan] << " = " << elkan_lower_bounds[folan][filan] << std::endl;
                    //     std::cout << "HERE: fake_smallest_ub: " << fake_smallest_lb << " ham_lb[folan] " << hamerly_lower_bounds[folan] << std::endl;
                    //     std::cout << "----------------" << std::endl;
                    // }

                }
            } else {
                for (filan = 0; filan < K; filan++) {
                    if (filan == smallest_ub[folan]) continue;
                    if (!is_candidate[folan][filan]) continue;
                    if(!candidates_exist){
                        candidates_exist = true;
                        // fake_smallest_lb = DBL_MAX;
                        hamerly_lower_bounds[folan] = DBL_MAX;
                    }
                    if (lower_bounds[folan][filan] >= upper_bounds[folan][smallest_ub[folan]]) {
                        is_candidate[folan][filan] = false;
                        // V9
                        // we should also try to find the second closest centroid
                        // TODO: im not sure if i should clear the old ham_lbs from last round
                        // what if it's masked completely, then we would never get here, so the old one should still be valid.
                        // nvm, i think im right, i won't clear it
                        // if (lower_bounds[folan][filan] < hamerly_lower_bounds[folan]) {
                        //     hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                        // }
                        // V9
                        // I'm fairly sure the above is bullshit, can't remember what I was thinking
                        // so ill make another
                        
                        
                    }
                    /*if(folan == 0){
                        std::cout << "before: HERE: filan: " << filan << " lower_bounds[folan][filan] " << lower_bounds[folan][filan] << " = " << elkan_lower_bounds[folan][filan] << std::endl;
                        std::cout << "HERE: ham_lb[folan] " << hamerly_lower_bounds[folan] << std::endl;
                        std::cout << "----------------" << std::endl;
                    }*/
                    
                    // find sec smallest lb
                    // if(lower_bounds[folan][filan] < fake_smallest_lb){
                    //     // keep sec smallest
                    //     hamerly_lower_bounds[folan] = fake_smallest_lb;
                    //     fake_smallest_lb = lower_bounds[folan][filan];
                    // } else if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                    //     hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                    // }
                    
                    // change of ham_lb definition

                    if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                        hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                    }

                    /*if(folan == 0){
                        std::cout << "after: HERE: filan: " << filan << " lower_bounds[folan][filan] " << lower_bounds[folan][filan] << " = " << elkan_lower_bounds[folan][filan] << std::endl;
                        std::cout << "HERE: ham_lb[folan] " << hamerly_lower_bounds[folan] << std::endl;
                        std::cout << "----------------" << std::endl;
                    }*/
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

    //std::cout << "set labels " << std::endl;
        

    // sanity check: is everyone assigned?
    /*for(folan = 0; folan < N; folan++){
        if(assigned[folan] == 0){
            std::cout << folan << " NOT ASSIGNED" << std::endl;
            exit(2);
        }
    }*/

    //std::cout << "after all-assigned sanity check" << std::endl;

    // calc new centroids
    // make copy of old centroids
    // apparently this also does not work, copies pointer somehow I think...
    // memcpy(old_centroids, centroids, sizeof(double) * K * D);
    for (folan = 0; folan < K; folan++) {
        for (filan = 0; filan < D; filan++) {
            old_centroids[folan][filan] = centroids[folan][filan];
        }
    }
    //std::cout << "copied centroids to old centroids" << std::endl;
    // set centroids to 0
    memset(cluster_counts, 0, sizeof(int) * K);
    //std::cout << "set cluster counts to 0" << std::endl;
    for (folan = 0; folan < K; folan++) {
        // just testing
        // cluster_counts[folan] = 0;
        for (filan = 0; filan < D; filan++) {
            centroids[folan][filan] = 0.0;
        }
    }
    // This doesn't work on doubles
    // memset(centroids, 0, sizeof(double) * K * D);
    //std::cout << "after all the memcpys" << std::endl;

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
    //std::cout << "calculated new centroids" << std::endl;
    calculate_centroids_square_sums();
    //std::cout << "ssed them " << std::endl;
    // just to check
    //int sanity_check = 0;
    /*std::cout << "cluster counts..." << std::endl;
    for (folan = 0; folan < K; folan++) {
        sanity_check += cluster_counts[folan];
        std::cout << cluster_counts[folan] << " ";
    }
    std::cout << sanity_check << std::endl;*/

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
    //std::cout << "calculated centroid movements" << std::endl;

    
    /*std::cout << "ham_ub[1] " << hamerly_upper_bounds[1] << " ham_lb[1] " << hamerly_lower_bounds[1] <<  std::endl;
    for(filan =0; filan < K; filan++){
        std::cout << elkan_lower_bounds[1][filan] << " ";
    }
    std::cout << std::endl;*/

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

    //std::cout << "after updating the bounds...\n";

    /*std::cout << "ham_ub[1] " << hamerly_upper_bounds[1] << " ham_lb[1] " << hamerly_lower_bounds[1] <<  std::endl;
    for(filan =0; filan < K; filan++){
        std::cout << elkan_lower_bounds[1][filan] << " ";
    }
    std::cout << std::endl;*/




    

    // loop over max_iter
    for (int iter = 1; iter < MAX_ITERATIONS; iter++) {
        level = 0;
        //std::cout << "iteration " << iter << "..." << std::endl;

        // calculate closest centroid to each centroid
        for (folan = 0; folan < K; folan++) {
            smallest = DBL_MAX;
            for (filan = folan; filan < K; filan++) {
                tmp = centroids_ss[filan][0] + centroids_ss[folan][0];
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
        //std::cout << "found closest distance to each centroid, and the centr-centr distances" << std::endl;


        // fill out mask
        // set all to 1 first
        for (folan = 0; folan < N; folan++) {
            for (filan = 0; filan < K; filan++) {
                is_candidate[folan][filan] = true;
            }
        }
        //std::cout << "filled out is_candid" << std::endl;


        // set the assigned
        memset(assigned, 0, sizeof(int) * N);
        //std::cout << "memset assigned" << std::endl;

        for(folan = 0; folan < N; folan++){
            for(filan=0; filan <K; filan++){
                last_level_calculated[folan][filan] = -1;
            }
        }
        //std::cout << "set last_level_calc to -1" << std::endl;

        // I think I have to put the step-wise while here

        // TODO:
        // for x in X:
        //  if ham_bound[x] is true: set is_candid to false for all centroid != labels[x] and set assigned[x] = 1
        for (folan = 0; folan < N; folan++) {
            // hamerly check
            hamerly_bound = ((0.5 * closest_centroid_distance[labels[folan]]) > hamerly_lower_bounds[folan]) ? (
                    0.5 * closest_centroid_distance[labels[folan]]) : hamerly_lower_bounds[folan];
            if (hamerly_bound >= hamerly_upper_bounds[folan]) {
                hamerly_count += K;
                // this one will not move, so none of them are candidates, it is assigned...
                assigned[folan] = 1;
            }
        }

        while (level < int(log2(int(sqrt(D))) + 2)) {
            //std::cout << "in while level " << level << std::endl;
            // set smallest ub
            // let's initiate them to the previous labels so skipping the labels in the for K does not change anything
            // TODO!!!!
            // memset(smallest_ub, 0, sizeof(int) * N);

            if(pow(2, 2 * level) < HYBRID_SWITCH_THRESHOLD){
                // do all the TI checks
                for (folan = 0; folan < N; folan++) {
                    // TODO: check if sth goes wrong
                    smallest_ub[folan] = labels[folan];
                    if (assigned[folan] > 0) continue;
                    r = true;
                    r_int = 0; fake_smallest_lb = DBL_MAX;
                    // hamerly_bound = ((0.5 * closest_centroid_distance[labels[folan]]) > hamerly_lower_bounds[folan]) ? (
                    //         0.5 * closest_centroid_distance[labels[folan]]) : hamerly_lower_bounds[folan];
                    // elkan lemma 1 and hamerly
                    // if (hamerly_bound < hamerly_upper_bounds[folan]) {
                    
                    
                    // if(true){
                    /*if(folan == 4){
                        std::cout << "not pruned by ham_bound..." << std::endl;
                    }*/
                    // if((0.5 * closest_centroid_distance[labels[folan]]) < hamerly_upper_bounds[folan]){
                    
                    // update the upper bound
                    // im gonna do it till this level only
                    calculate_distance_folan_filan_till_level(folan, labels[folan], level);
                    // update the ham_ub and elkan_lb if better
                    // TODO: not sure, maybe we have to update it anyway
                    // hamerly_upper_bounds[folan] = std::min(hamerly_upper_bounds[folan], upper_bounds[folan][labels[folan]]);
                    // elkan_lower_bounds[folan][labels[folan]] = std::max(elkan_lower_bounds[folan][labels[folan]], lower_bounds[folan][labels[folan]]);
                    // AS DISCUSSED WITH KASPER, lets's try it this way for now TODO: check
                    hamerly_upper_bounds[folan] = upper_bounds[folan][labels[folan]];
                    elkan_lower_bounds[folan][labels[folan]] = lower_bounds[folan][labels[folan]];
                    /*if(folan == 4){
                        std::cout << "JUST FOR THE SAKE OF FATEMEH'S SANITY: euc_dist(1, 1) " << sqrt(euclidean_distance(1, 1)) << std::endl; 
                        std::cout << "1: changed ham_ub[" << folan << "] to " << hamerly_upper_bounds[folan];
                        std::cout << " and elkan_lb[" << folan << "][" << labels[folan] << "] to " << elkan_lower_bounds[folan][labels[folan]] << std::endl;
                    }*/
                    
                    for (filan = 0; filan < K; filan++) {
                        // TODO: tentative...
                        if(!is_candidate[folan][filan]) continue;
                        // if(folan == 1){
                        //     std::cout << "am i going crazy? " << filan << " " << labels[folan] << std::endl;
                        // }
                        if (filan == labels[folan]) continue;
                        // new lemma
                        
                        /*if(folan == 4){
                            std::cout << "not pruned by new lemma " << filan << std::endl;
                        }*/
                        if (hamerly_upper_bounds[folan] > elkan_lower_bounds[folan][filan]
                            && hamerly_upper_bounds[folan] >
                                (0.5 * centroid_to_centroid_distances[filan][labels[folan]])) {
                                // HERE
                                // I'm thinking the full dist calculation of the prev lables[folan]
                                // we have to do it all in one go
                                // or the r bool has to be activated for all the levels or sth
                                
                                // if(r_int < int(log2(int(sqrt(D))) + 2)){
                                
                                // TODO: maybe I should have a if(is_candid) around this
                                // we do the TI pruning with the ifs, then the step-wise prunings with the is_candid
                                // or I could just have a if(!is_candid) continue at the start of for K ?
                                // I think this could work
                                calculate_distance_folan_filan_till_level(folan, filan, level);
                                // i'm not sure about this maxing, maybe I should always set it to lb TODO?
                                elkan_lower_bounds[folan][filan] = lower_bounds[folan][filan];
                                /*if(folan == 4){
                                    std::cout << "2: changed  elkan_lb[" << folan << "][" << filan << "] to " << elkan_lower_bounds[folan][filan] << " ub was " << upper_bounds[folan][filan] << std::endl;
                                } */   
                                // elkan_lower_bounds[folan][filan] = std::max(elkan_lower_bounds[folan][filan], lower_bounds[folan][filan]);

                                // then (in non-step-wise) we check if the dist is smaller than the one for labels[folan]
                                // I'm not sure we can do that since we don't exactly have distances, but ub and lb
                                // we also don't have labels[folan], in the same way as before
                                // we could store a smallest_ub_till_now
                                // then use that as the labels[folan]
                                // we could even set the labels[folan] = smallest_ub_till_now and update as we go
                                // then we'd do the ic pruning here to avoid an extra loop
                                // no, that wouldn't really work. since the smallest_ub_till_now is not necessarily the smallest_ub over all
                                // we don't need a separate variable for smallest_ub, it could just be labels[folan]
                                // then the problem is when do we say it is assigned.
                                // TODO
                                // maybe I should update ham_lb separately
                                // like, if they are not from the same point, it would still be correct
                                // I think...
                                // yeah but then ham_ub==labels[folan] would be the smallest_ub
                                // but ham_lb is supposed to be the second smallest lb:-?
                                // or the second largest lb;-? 
                                // I don't think it should be sec largest lb
                                // i'll go with sec smallest lb, fingers crossed --> for this I will have to keep smallest_lb too
                                if(upper_bounds[folan][filan] < hamerly_upper_bounds[folan]){
                                    // TODO: how the f am i supposed to update ham_lb now?!
                                    // update the ham_lb if necessary
                                    // hamerly_lower_bounds[folan] = std::max(hamerly_lower_bounds[folan], hamerly_upper_bounds[folan]);
                                    // hamerly_lower_bounds[folan] = hamerly_upper_bounds[folan];
                                    labels[folan] = filan;
                                    hamerly_upper_bounds[folan] = upper_bounds[folan][filan];
                                    /*if(folan == 4){
                                        std::cout << "3: changed  ham_ub[" << folan << "] to " << hamerly_upper_bounds[folan];
                                        std::cout << " and labels[" << folan << "]  to " << labels[folan] << std::endl;
                                    }*/
                                }
                                // I have chnaged the def of ham_lb
                                // I will only update the ham_lb in the pruning of stepwise
                                // I only use the ham_lb in the ham_bound if, so I can just update it after the whole for K
                                // will it ever stay unassigned?
                                // after this for K
                                // we go into the step-wise pruning
                                // then, if there are candidates, the smallest lb of the candidates(other than the smallest ub) will be ham_lb
                                // then if there are no candidates, the ham_lb of last level remains
                                // the first level will always have candidates(everybody is a candidate)
                                // so we will always assign it
                                // no im doubting if the def is correct
                                // since the "the smallest lb of the candidates(other than the smallest ub)"
                                // does not really add up to much
                                // idk, let's see what happens:) 
                                // finding sec smallest lb
                                // if(lower_bounds[folan][filan] < fake_smallest_lb){
                                //     // TODO: maybe the max completely ruins the whole thing...
                                //     // i'm not gonna do the max for now
                                //     // hamerly_lower_bounds[folan] = std::max(hamerly_lower_bounds[folan], fake_smallest_lb);
                                //     hamerly_lower_bounds[folan] = fake_smallest_lb;
                                //     fake_smallest_lb = lower_bounds[folan][filan];
                                // } else if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                                //     hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                                // }
                        } 
                        // FOR NOW I WON'T
                        /*else{
                            // TODO: I think I need to change the is_candid
                            if(folan == 4){
                                std::cout << "3.5 skipping this " << filan << std::endl;
                            }
                        }*/

                        // im gonna update smallest_ub here
                        // I think it should be fine
                        // cause after all the hassle, the ub and lb here should be valid
                        if(upper_bounds[folan][filan] < upper_bounds[folan][smallest_ub[folan]]){ smallest_ub[folan] = filan; }
                    } // for filan < K
                    // } // if ham_bound
                    
                    // TODO: now we have to do the pruning based off step-wise... I think...
                    // I'll move this to the start of the N loop
                    // I added it there, but will keep it here because of the else that can change assigned with ham_boung
                    if (assigned[folan] > 0) continue;
                    /*if(folan == 4){
                        std::cout << "elkan_lb[4] ";
                        for(filan = 0; filan < K; filan++){
                            std::cout << elkan_lower_bounds[folan][filan] << " ";
                        }
                        std::cout << std::endl;
                    }*/

                    candidates_exist = false;
                    // V9
                    if (level == int(log2(int(sqrt(D))) + 1)) {
                        labels[folan] = smallest_ub[folan];
                        // V9
                        hamerly_upper_bounds[folan] = upper_bounds[folan][smallest_ub[folan]];
                        assigned[folan] = 1;
                        // v9

                        // but we still need to update the ham_lb
                        for (filan = 0; filan < K; filan++) {
                            if (filan == smallest_ub[folan]) continue;
                            if (!is_candidate[folan][filan]) continue;
                            if(!candidates_exist){
                                candidates_exist = true;
                                // fake_smallest_lb = DBL_MAX;
                                hamerly_lower_bounds[folan] = DBL_MAX;
                            }   

                            // TODO: so im changing the definition of ham_lb
                            // let ham_lb be the smallest lb betweent the last batch of candidates
                            // if there are other candidates at THIS point, then the ham_lb should just be fake_smallest_lb here
                            // but what if there are no other candidates?
                            // then we would have set fake_smallest_lb to DBL_MAX and it stays that way...
                            // I fixed it with the exists_candidate
                            // it's not the label == smallest_ub
                            // so it's the smallest lb to get pruned
                            if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                                hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                            }

                        }

                    } else {
                        for (filan = 0; filan < K; filan++) {
                            if (filan == smallest_ub[folan]) continue;
                            if (!is_candidate[folan][filan]) continue;
                            
                            if (lower_bounds[folan][filan] >= upper_bounds[folan][smallest_ub[folan]]) {
                                is_candidate[folan][filan] = false;
                                if(!candidates_exist){
                                    candidates_exist = true;
                                    // fake_smallest_lb = DBL_MAX;
                                    hamerly_lower_bounds[folan] = DBL_MAX;
                                }
                                /*if(folan == 4){
                                    std::cout << "4: pruning " << filan << std::endl;
                                }*/
                                // NOT HERE i THINK
                                // V9
                                // we should also try to find the second closest centroid
                                // TODO: im not sure if i should clear the old ham_lbs from last round
                                // what if it's masked completely, then we would never get here, so the old one should still be valid.
                                // nvm, i think im right, i won't clear it
                                // if (lower_bounds[folan][filan] < hamerly_lower_bounds[folan]) {
                                //     hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                                // }
                                // V9
                                if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                                    hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                                    /*if(folan == 4){
                                        std::cout << "5: chnaged ham_lb[" << folan << "] to " << hamerly_lower_bounds[folan] << std::endl; 
                                    }*/
                                }
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
                            /*if(folan == 4){
                                std::cout << "6: only one candidate chnaged labels[" << folan << "] to " << labels[folan];
                                std::cout << " and ham_ub[" << folan << "] to " << hamerly_upper_bounds[folan] << std::endl; 
                            }*/
                        }
                    }
                    
                } // for folan < N
            }
            else{
                // no extra TI checks, only fill out the Ham bounds so we can use them to set assigneds
                calculate_sqrt_distances_till_level_with_assigned_and_ll(level);
                memset(smallest_ub, 0, sizeof(int) * N);
                for (folan = 0; folan < N; folan++) {
                    for (filan = 0; filan < K; filan++) {
                        if(!is_candidate[folan][filan]) continue;
                        if (upper_bounds[folan][filan] < upper_bounds[folan][smallest_ub[folan]]) smallest_ub[folan] = filan;
                    }
                }
                for(folan = 0; folan < N; folan++){
                    if (assigned[folan] > 0) continue;
                    
                    candidates_exist = false;
                    // V9
                    if (level == int(log2(int(sqrt(D))) + 1)) {
                        labels[folan] = smallest_ub[folan];
                        // V9
                        hamerly_upper_bounds[folan] = upper_bounds[folan][smallest_ub[folan]];
                        assigned[folan] = 1;
                        // v9

                        // but we still need to update the ham_lb
                        for (filan = 0; filan < K; filan++) {
                            if (filan == smallest_ub[folan]) continue;
                            if (!is_candidate[folan][filan]) continue;
                            if(!candidates_exist){
                                candidates_exist = true;
                                // fake_smallest_lb = DBL_MAX;
                                hamerly_lower_bounds[folan] = DBL_MAX;
                            }

                            // TODO: so im changing the definition of ham_lb
                            // let ham_lb be the smallest lb betweent the last batch of candidates
                            // if there are other candidates at THIS point, then the ham_lb should just be fake_smallest_lb here
                            // but what if there are no other candidates?
                            // then we would have set fake_smallest_lb to DBL_MAX and it stays that way...
                            // I fixed it with the exists_candidate
                            // it's not the label == smallest_ub
                            // so it's the smallest lb to get pruned
                            if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                                hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                            }

                        }

                    } else {
                        for (filan = 0; filan < K; filan++) {
                            if (filan == smallest_ub[folan]) continue;
                            if (!is_candidate[folan][filan]) continue;
                            // if(!candidates_exist){
                            //     candidates_exist = true;
                            //     // fake_smallest_lb = DBL_MAX;
                            //     hamerly_lower_bounds[folan] = DBL_MAX;
                            // }
                            if (lower_bounds[folan][filan] >= upper_bounds[folan][smallest_ub[folan]]) {
                                is_candidate[folan][filan] = false;
                                /*if(folan == 4){
                                    std::cout << "4: pruning " << filan << std::endl;
                                }*/
                                if(!candidates_exist){
                                    candidates_exist = true;
                                    // fake_smallest_lb = DBL_MAX;
                                    hamerly_lower_bounds[folan] = DBL_MAX;
                                }
                                // NOT HERE i THINK
                                // V9
                                // we should also try to find the second closest centroid
                                // TODO: im not sure if i should clear the old ham_lbs from last round
                                // what if it's masked completely, then we would never get here, so the old one should still be valid.
                                // nvm, i think im right, i won't clear it
                                // if (lower_bounds[folan][filan] < hamerly_lower_bounds[folan]) {
                                //     hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                                // }
                                // V9
                                if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                                    hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                                    /*if(folan == 4){
                                        std::cout << "5: chnaged ham_lb[" << folan << "] to " << hamerly_lower_bounds[folan] << std::endl; 
                                    }*/
                                }
                            }
                            
                            // if(lower_bounds[folan][filan] < hamerly_lower_bounds[folan]){
                            //     hamerly_lower_bounds[folan] = lower_bounds[folan][filan];
                            //     if(folan == 4){
                            //         std::cout << "5: chnaged ham_lb[" << folan << "] to " << hamerly_lower_bounds[folan] << std::endl; 
                            //     }
                            // }
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
                            /*if(folan == 4){
                                std::cout << "6: only one candidate chnaged labels[" << folan << "] to " << labels[folan];
                                std::cout << " and ham_ub[" << folan << "] to " << hamerly_upper_bounds[folan] << std::endl; 
                            }*/
                        }
                    }
                    
                } // for folan < N

            }
            level++;
        } // step-wise loop
        
        // TODO: sanity check, is everybody assigned?
        /*for(folan = 0; folan < N; folan++){
            if(assigned[folan] == 0){
                std::cout << folan << " NOT ASSIGNED" << std::endl;
                exit(2);
            }
        }

        std::cout << "set labels" << std::endl;*/


        // calc new centroids
        // make copy of old centroids
        // apparently this also does not work, copies pointer somehow I think...
        // memcpy(old_centroids, centroids, sizeof(double) * K * D);
        for (folan = 0; folan < K; folan++) {
            for (filan = 0; filan < D; filan++) {
                old_centroids[folan][filan] = centroids[folan][filan];
            }
        }
        //std::cout << "copied centroids to old centroids" << std::endl;
        // set centroids to 0
        memset(cluster_counts, 0, sizeof(int) * K);
        //std::cout << "set cluster counts to 0" << std::endl;
        for (folan = 0; folan < K; folan++) {
            // just testing
            // cluster_counts[folan] = 0;
            for (filan = 0; filan < D; filan++) {
                centroids[folan][filan] = 0.0;
            }
        }
        // This doesn't work on doubles
        // memset(centroids, 0, sizeof(double) * K * D);
        //std::cout << "after all the memcpys" << std::endl;

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
        //std::cout << "calculated new centroids" << std::endl;
        calculate_centroids_square_sums();
        //std::cout << "ssed them " << std::endl;
        // just to check
        /*int sanity_check = 0;
        std::cout << "cluster counts..." << std::endl;
        for (folan = 0; folan < K; folan++) {
            sanity_check += cluster_counts[folan];
            std::cout << cluster_counts[folan] << " ";
        }
        std::cout << sanity_check << std::endl;*/

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
        //std::cout << "calculated centroid movements" << std::endl;

        /*std::cout << "ham_ub[1] " << hamerly_upper_bounds[1] << " ham_lb[1] " << hamerly_lower_bounds[1] << std::endl;
        for(filan = 0; filan < K; filan++){
            std::cout << elkan_lower_bounds[1][filan] << " " ;
        }
        std::cout << std::endl;
        */

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
        //std::cout << "checked convergence" << std::endl;

        // end if converged
        if (has_converged) {
            if (D >= 256*256) {
                int sanity_check = 0;
                for (folan = 0; folan < K; folan++) {
                    sanity_check += cluster_counts[folan];
                    std::cout << cluster_counts[folan] << " ";
                }
                std::cout << sanity_check << std::endl;
            }
            break;
        }

    } // for iter

}


// hamerly + step-wise + elkan Mohammad-style, without ham_lb
// copied from v9
// TODO: DONE
void kmeans_v12() {
    //std::cout << "in kmeans_v12 ..." << std::endl;
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
    //std::cout << "copied init centroids" << std::endl;


    bool has_converged = false;
    int cluster_counts[K];
    // V2 DIFF start
    double centroid_squares[K];
    // V2 DIFF end

    // loop over max_iter
    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
        //std::cout << "iteration " << iter << "..." << std::endl;

        // calculate the square sum of centroids
        // for (folan = 0; folan < K; folan++) {
        //     centroid_squares[folan] = 0;
        //     for (filan = 0; filan < D; filan++) {
        //         centroid_squares[folan] += (centroids[folan][filan] * centroids[folan][filan]);
        //     }
        // }

        // calculate closest centroid to each centroid
        for (folan = 0; folan < K; folan++) {
            smallest = DBL_MAX;
            for (filan = folan; filan < K; filan++) {
                tmp = centroids_ss[filan][0] + centroids_ss[folan][0];
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
        //std::cout << "found closest distance to each centroid, and the centr-centr distances" << std::endl;


        if (iter == 0) {
            //std::cout << "in if iter == 0" << std::endl;
            
            memset(assigned, 0, sizeof(int) * N);
            // calculate_labels_with_sqrt_hamerly_elkan_integrated();
            calculate_labels_with_sqrt();
            //std::cout << "calculated distances" << std::endl;
            //std::cout << "filled out hamerly ub and lb for first time" << std::endl;

            // let's fill the ham_ub and elkan_lb
            for(folan = 0; folan < N; folan++){
                hamerly_upper_bounds[folan] = upper_bounds[folan][labels[folan]];
                for(filan = 0; filan < K; filan++){
                    elkan_lower_bounds[folan][filan] = lower_bounds[folan][filan];
                }
            }

        } else {
            //std::cout << "iter was not 0" << std::endl;
            // V9
            // fill out mask
            // set all to 1 first
            for (folan = 0; folan < N; folan++) {
                for (filan = 0; filan < K; filan++) {
                    is_candidate[folan][filan] = true;
                }
            }

            memset(assigned, 0, sizeof(int) * N);
            // now let's do pruning based on TI
            for (folan = 0; folan < N; folan++) {
                // hamerly check
                hamerly_bound = ((0.5 * closest_centroid_distance[labels[folan]]) > hamerly_lower_bounds[folan]) ? (
                        0.5 * closest_centroid_distance[labels[folan]]) : hamerly_lower_bounds[folan];
                // if (hamerly_bound >= hamerly_upper_bounds[folan]) {
                if ((0.5 * closest_centroid_distance[labels[folan]]) >= hamerly_upper_bounds[folan]) {
                    hamerly_count += K;
                    // this one will not move, so none of them are candidates
                    assigned[folan] = 1;
                    continue;
                } 
                else  {
                    for(filan = 0; filan < K; filan++){
                        if (filan == labels[folan]) continue;
                        if (hamerly_upper_bounds[folan] <= elkan_lower_bounds[folan][filan]
                            || hamerly_upper_bounds[folan] <=
                                (0.5 * centroid_to_centroid_distances[filan][labels[folan]])){
                                    is_candidate[folan][filan] = false;
                                }
                    }
                }
            }
            //std::cout << "hamerly pruned " << hamerly_count << std::endl;
            
            // calculate_labels_with_sqrt_hamerly_elkan_integrated();
            calculate_labels_with_sqrt();

            // let's fill the ham_ub and elkan_lb
            for(folan = 0; folan < N; folan++){
                hamerly_upper_bounds[folan] = upper_bounds[folan][labels[folan]];
                for(filan = 0; filan < K; filan++){
                    elkan_lower_bounds[folan][filan] = lower_bounds[folan][filan];
                }
            }

        }

        //std::cout << "set labels" << std::endl;




        // calc new centroids
        // make copy of old centroids
        // apparently this also does not work, copies pointer somehow I think...
        // memcpy(old_centroids, centroids, sizeof(double) * K * D);
        for (folan = 0; folan < K; folan++) {
            for (filan = 0; filan < D; filan++) {
                old_centroids[folan][filan] = centroids[folan][filan];
            }
        }
        //std::cout << "copied centroids to old centroids" << std::endl;
        // set centroids to 0
        memset(cluster_counts, 0, sizeof(int) * K);
        //std::cout << "set cluster counts to 0" << std::endl;
        for (folan = 0; folan < K; folan++) {
            // just testing
            // cluster_counts[folan] = 0;
            for (filan = 0; filan < D; filan++) {
                centroids[folan][filan] = 0.0;
            }
        }
        // This doesn't work on doubles
        // memset(centroids, 0, sizeof(double) * K * D);
        //std::cout << "after all the memcpys" << std::endl;

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
        //std::cout << "calculated new centroids" << std::endl;
        calculate_centroids_square_sums();

        // just to check
        /*int sanity_check = 0;
        std::cout << "cluster counts..." << std::endl;
        for (folan = 0; folan < K; folan++) {
            sanity_check += cluster_counts[folan];
            std::cout << cluster_counts[folan] << " ";
        }
        std::cout << sanity_check << std::endl;*/

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
        //std::cout << "calculated centroid movements" << std::endl;

        // update upper and lower hamerly bounds based on centroid movements
        // for (folan = 0; folan < N; folan++) {
        //     hamerly_upper_bounds[folan] += centroid_movement[labels[folan]];
        //     if (labels[folan] == furthest_moving_centroid) {
        //         hamerly_lower_bounds[folan] -= centroid_movement[second_furthest_moving_centroid];
        //     } else {
        //         hamerly_lower_bounds[folan] -= centroid_movement[furthest_moving_centroid];
        //     }

        //     for (filan = 0; filan < K; filan++) {
        //         elkan_lower_bounds[folan][filan] -= centroid_movement[filan];
        //     }
        // }
        // FROM MOHAMMAD'S PSEUDO
        for (folan = 0; folan < N; folan++) {
            hamerly_upper_bounds[folan] += centroid_movement[furthest_moving_centroid];
            hamerly_lower_bounds[folan] -= centroid_movement[furthest_moving_centroid];
            for (filan = 0; filan < K; filan++) {
                elkan_lower_bounds[folan][filan] -= centroid_movement[furthest_moving_centroid];
            }
        }


        // FATEMEH DEBUG
        //std::cout << "after updating the bounds, hamerly_ub[0] " << hamerly_upper_bounds[0] << " hamerly_lb[0] "
        //     << hamerly_lower_bounds[0] << std::endl;
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
        //std::cout << "checked convergence" << std::endl;

        // end if converged
        if (has_converged) {
            if (D >= 256*256) {
                int sanity_check = 0;
                for (folan = 0; folan < K; folan++) {
                    sanity_check += cluster_counts[folan];
                    std::cout << cluster_counts[folan] << " ";
                }
                std::cout << sanity_check << std::endl;
            }
            break;
        }
    }

    //std::cout << "hamerly pruned " << hamerly_count << std::endl;
}



// hamerly + step-wise + elkan Mohammad-style, with ham_lb
// copied from v12
// TODO
void kmeans_v121() {
    //std::cout << "in kmeans_v121 ..." << std::endl;
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
    //std::cout << "copied init centroids" << std::endl;


    bool has_converged = false;
    int cluster_counts[K];
    // V2 DIFF start
    double centroid_squares[K];
    // V2 DIFF end

    // loop over max_iter
    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
        //std::cout << "iteration " << iter << "..." << std::endl;

        // calculate the square sum of centroids
        // for (folan = 0; folan < K; folan++) {
        //     centroid_squares[folan] = 0;
        //     for (filan = 0; filan < D; filan++) {
        //         centroid_squares[folan] += (centroids[folan][filan] * centroids[folan][filan]);
        //     }
        // }

        // calculate closest centroid to each centroid
        for (folan = 0; folan < K; folan++) {
            smallest = DBL_MAX;
            for (filan = folan; filan < K; filan++) {
                tmp = centroids_ss[filan][0] + centroids_ss[folan][0];
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
        //std::cout << "found closest distance to each centroid, and the centr-centr distances" << std::endl;


        if (iter == 0) {
            //std::cout << "in if iter == 0" << std::endl;
            
            memset(assigned, 0, sizeof(int) * N);
            // calculate_labels_with_sqrt_hamerly_elkan_integrated();
            calculate_labels_with_sqrt();
            //std::cout << "calculated distances" << std::endl;
            //std::cout << "filled out hamerly ub and lb for first time" << std::endl;

            // let's fill the ham_ub and elkan_lb
            for(folan = 0; folan < N; folan++){
                hamerly_upper_bounds[folan] = upper_bounds[folan][labels[folan]];
                second_smallest = DBL_MAX;
                for(filan = 0; filan < K; filan++){
                    elkan_lower_bounds[folan][filan] = lower_bounds[folan][filan];
                    // this will be the smallest of lbs besides the lb of the label, so it should be the second smallest lb overall
                    if(elkan_lower_bounds[folan][filan] < second_smallest && filan != labels[folan]){
                        second_smallest = elkan_lower_bounds[folan][filan];
                    }
                }
                hamerly_lower_bounds[folan] = second_smallest;
            }

        } else {
            //std::cout << "iter was not 0" << std::endl;
            // V9
            // fill out mask
            // set all to 1 first
            for (folan = 0; folan < N; folan++) {
                for (filan = 0; filan < K; filan++) {
                    is_candidate[folan][filan] = true;
                }
            }

            memset(assigned, 0, sizeof(int) * N);
            // now let's do pruning based on TI
            for (folan = 0; folan < N; folan++) {
                // hamerly check
                hamerly_bound = ((0.5 * closest_centroid_distance[labels[folan]]) > hamerly_lower_bounds[folan]) ? (
                        0.5 * closest_centroid_distance[labels[folan]]) : hamerly_lower_bounds[folan];
                if (hamerly_bound >= hamerly_upper_bounds[folan]) {
                // if ((0.5 * closest_centroid_distance[labels[folan]]) >= hamerly_upper_bounds[folan]) {
                    hamerly_count += K;
                    // this one will not move, so none of them are candidates
                    assigned[folan] = 1;
                    continue;
                } 
                else  {
                    for(filan = 0; filan < K; filan++){
                        if (filan == labels[folan]) continue;
                        if (hamerly_upper_bounds[folan] <= elkan_lower_bounds[folan][filan]
                            || hamerly_upper_bounds[folan] <=
                                (0.5 * centroid_to_centroid_distances[filan][labels[folan]])){
                                    is_candidate[folan][filan] = false;
                                }
                    }
                }
            }
            //std::cout << "hamerly pruned " << hamerly_count << std::endl;
            
            // calculate_labels_with_sqrt_hamerly_elkan_integrated();
            calculate_labels_with_sqrt();

            // let's fill the ham_ub and elkan_lb
            // for(folan = 0; folan < N; folan++){
            //     hamerly_upper_bounds[folan] = upper_bounds[folan][labels[folan]];
            //     for(filan = 0; filan < K; filan++){
            //         elkan_lower_bounds[folan][filan] = lower_bounds[folan][filan];
            //     }
            // }

            for(folan = 0; folan < N; folan++){
                hamerly_upper_bounds[folan] = upper_bounds[folan][labels[folan]];
                second_smallest = DBL_MAX;
                for(filan = 0; filan < K; filan++){
                    elkan_lower_bounds[folan][filan] = lower_bounds[folan][filan];
                    // this will be the smallest of lbs besides the lb of the label, so it should be the second smallest lb overall
                    if(elkan_lower_bounds[folan][filan] < second_smallest && filan != labels[folan]){
                        second_smallest = elkan_lower_bounds[folan][filan];
                    }
                }
                hamerly_lower_bounds[folan] = second_smallest;
            }

        }

        //std::cout << "set labels" << std::endl;




        // calc new centroids
        // make copy of old centroids
        // apparently this also does not work, copies pointer somehow I think...
        // memcpy(old_centroids, centroids, sizeof(double) * K * D);
        for (folan = 0; folan < K; folan++) {
            for (filan = 0; filan < D; filan++) {
                old_centroids[folan][filan] = centroids[folan][filan];
            }
        }
        //std::cout << "copied centroids to old centroids" << std::endl;
        // set centroids to 0
        memset(cluster_counts, 0, sizeof(int) * K);
        //std::cout << "set cluster counts to 0" << std::endl;
        for (folan = 0; folan < K; folan++) {
            // just testing
            // cluster_counts[folan] = 0;
            for (filan = 0; filan < D; filan++) {
                centroids[folan][filan] = 0.0;
            }
        }
        // This doesn't work on doubles
        // memset(centroids, 0, sizeof(double) * K * D);
        //std::cout << "after all the memcpys" << std::endl;

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
        //std::cout << "calculated new centroids" << std::endl;
        calculate_centroids_square_sums();

        // just to check
        /*int sanity_check = 0;
        std::cout << "cluster counts..." << std::endl;
        for (folan = 0; folan < K; folan++) {
            sanity_check += cluster_counts[folan];
            std::cout << cluster_counts[folan] << " ";
        }
        std::cout << sanity_check << std::endl;*/

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
        //std::cout << "calculated centroid movements" << std::endl;

        // update upper and lower hamerly bounds based on centroid movements
        // for (folan = 0; folan < N; folan++) {
        //     hamerly_upper_bounds[folan] += centroid_movement[labels[folan]];
        //     if (labels[folan] == furthest_moving_centroid) {
        //         hamerly_lower_bounds[folan] -= centroid_movement[second_furthest_moving_centroid];
        //     } else {
        //         hamerly_lower_bounds[folan] -= centroid_movement[furthest_moving_centroid];
        //     }

        //     for (filan = 0; filan < K; filan++) {
        //         elkan_lower_bounds[folan][filan] -= centroid_movement[filan];
        //     }
        // }
        // FROM MOHAMMAD'S PSEUDO
        for (folan = 0; folan < N; folan++) {
            hamerly_upper_bounds[folan] += centroid_movement[furthest_moving_centroid];
            hamerly_lower_bounds[folan] -= centroid_movement[furthest_moving_centroid];
            for (filan = 0; filan < K; filan++) {
                elkan_lower_bounds[folan][filan] -= centroid_movement[furthest_moving_centroid];
            }
        }


        // FATEMEH DEBUG
        //std::cout << "after updating the bounds, hamerly_ub[0] " << hamerly_upper_bounds[0] << " hamerly_lb[0] "
        //     << hamerly_lower_bounds[0] << std::endl;
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
        //std::cout << "checked convergence" << std::endl;

        // end if converged
        if (has_converged) {
            if (D >= 256*256) {
                int sanity_check = 0;
                for (folan = 0; folan < K; folan++) {
                    sanity_check += cluster_counts[folan];
                    std::cout << cluster_counts[folan] << " ";
                }
                std::cout << sanity_check << std::endl;
            }
            break;
        }
    }

    //std::cout << "hamerly pruned " << hamerly_count << std::endl;
}


// STUFF FOR COMPLETE MARIGOLD END

};






static void BM_Kmeans(benchmark::State& state) {
    //std::cout << "reading state" << std::endl;
        //std::cout << state.range(0) << " " << state.range(1) << " " << state.range(2) << " " << std::endl;
        //std::cout << "done reading state" << std::endl;

    for (auto _ : state) { 
        state.PauseTiming();
        
        Kmeans_bench lloyd(state.range(0), state.range(1), state.range(2));
        //std::cout << "woop" << std::endl;
        lloyd.init(argc_, argv_); 
        //std::cout << "woop" << std::endl;
        state.ResumeTiming();
        //std::cout << "woop" << std::endl;
        lloyd.kmeans();
        //std::cout << "woop" << std::endl;
        state.PauseTiming();
        
        std::ofstream label_file;
        std::string file = "/mnt/c/Users/kaspe/OneDrive/Skrivebord/Reps/scalable-kmeans/cpp/result_Lloyd"+std::to_string(state.range(0))+"_"+std::to_string(state.range(1))+"_"+std::to_string(state.range(2))+"_"+".txt";
        label_file.open(file.c_str());

        for (int i = 0; i < state.range(1); i++) {
            label_file << lloyd.labels[i] << "\n";
        }
        label_file.close();

        lloyd.clean();

        state.ResumeTiming();
    }
}

static void BM_Kmeans_stepwise(benchmark::State& state) {
    

    for (auto _ : state) { 
        state.PauseTiming();
        Kmeans_bench v4(state.range(0), state.range(1), state.range(2));
        v4.init(argc_, argv_); 
        state.ResumeTiming();

        v4.kmeans_v4();

        state.PauseTiming();
        
        std::ofstream label_file;
        std::string file = "/mnt/c/Users/kaspe/OneDrive/Skrivebord/Reps/scalable-kmeans/cpp/result_Stepwise"+std::to_string(state.range(0))+"_"+std::to_string(state.range(1))+"_"+std::to_string(state.range(2))+"_"+".txt";
        label_file.open(file.c_str());

        for (int i = 0; i < state.range(1); i++) {
            label_file << v4.labels[i] << "\n";
        }
        label_file.close();

        v4.clean();

        state.ResumeTiming();

    }
}

static void BM_Kmeans_elkan(benchmark::State& state) { 

    for (auto _ : state) { 
        state.PauseTiming();
        Kmeans_bench v6(state.range(0), state.range(1), state.range(2));
        v6.init(argc_, argv_); 
        state.ResumeTiming();
        v6.kmeans_v6();

        state.PauseTiming();
        std::ofstream label_file;
        std::string file = "/mnt/c/Users/kaspe/OneDrive/Skrivebord/Reps/scalable-kmeans/cpp/result_Elkan"+std::to_string(state.range(0))+"_"+std::to_string(state.range(1))+"_"+std::to_string(state.range(2))+"_"+".txt";
        label_file.open(file.c_str());

        for (int i = 0; i < state.range(1); i++) {
            label_file << v6.labels[i] << "\n";
        }
        label_file.close();

        v6.clean();
        state.ResumeTiming();
    }
}

static void BM_Kmeans_hamerly(benchmark::State& state) { 

    for (auto _ : state) { 
        state.PauseTiming();
        Kmeans_bench v5(state.range(0), state.range(1), state.range(2));
        v5.init(argc_, argv_); 
        state.ResumeTiming();

        v5.kmeans_v5();

        state.PauseTiming();
        v5.clean();
        state.ResumeTiming();
    }
}

static void BM_Kmeans_ElkanNewTE(benchmark::State& state) { 

    for (auto _ : state) { 
        state.PauseTiming();
        Kmeans_bench v7(state.range(0), state.range(1), state.range(2));
        v7.init(argc_, argv_); 
        state.ResumeTiming();

        v7.kmeans_v7();

        state.PauseTiming();
        v7.clean();
        state.ResumeTiming();
    }
}

static void BM_Kmeans_ElkanHamerly(benchmark::State& state) { 

    for (auto _ : state) { 
        state.PauseTiming();
        Kmeans_bench v85(state.range(0), state.range(1), state.range(2));
        v85.init(argc_, argv_); 
        state.ResumeTiming();

        v85.kmeans_v85();

        state.PauseTiming();
        std::ofstream label_file;
        std::string file = "/mnt/c/Users/kaspe/OneDrive/Skrivebord/Reps/scalable-kmeans/cpp/result_ElkHam"+std::to_string(state.range(0))+"_"+std::to_string(state.range(1))+"_"+std::to_string(state.range(2))+"_"+".txt";
        label_file.open(file.c_str());

        for (int i = 0; i < state.range(1); i++) {
            label_file << v85.labels[i] << "\n";
        }
        label_file.close();

        v85.clean();
        state.ResumeTiming();
    }
}

static void BM_Kmeans_HamerlyStepwise(benchmark::State& state) { 

    for (auto _ : state) { 
        state.PauseTiming();
        Kmeans_bench v9(state.range(0), state.range(1), state.range(2));
        v9.init(argc_, argv_); 
        state.ResumeTiming();

        v9.kmeans_v9();

        state.PauseTiming();
        std::ofstream label_file;
        std::string file = "/mnt/c/Users/kaspe/OneDrive/Skrivebord/Reps/scalable-kmeans/cpp/result_HamStep"+std::to_string(state.range(0))+"_"+std::to_string(state.range(1))+"_"+std::to_string(state.range(2))+"_"+".txt";
        label_file.open(file.c_str());

        for (int i = 0; i < state.range(1); i++) {
            label_file << v9.labels[i] << "\n";
        }
        label_file.close();

        v9.clean();
        state.ResumeTiming();
    }
}

static void BM_Kmeans_HamerlyNewTE(benchmark::State& state) {

    for (auto _ : state) {
        state.PauseTiming();
        Kmeans_bench v7_5(state.range(0), state.range(1), state.range(2));
        v7_5.init(argc_, argv_);
        state.ResumeTiming();

        v7_5.kmeans_v7_5();

        state.PauseTiming();
        

        v7_5.clean();
        state.ResumeTiming();
    }
}

static void BM_Kmeans_MARIGOLD(benchmark::State& state) { 

    for (auto _ : state) { 
        state.PauseTiming();
        Kmeans_bench v105(state.range(0), state.range(1), state.range(2));
        v105.init(argc_, argv_); 
        state.ResumeTiming();

        v105.kmeans_v105();
        

        state.PauseTiming();
        
        std::ofstream label_file;
        std::string file = "/mnt/c/Users/kaspe/OneDrive/Skrivebord/Reps/scalable-kmeans/cpp/result_MARIGOLD"+std::to_string(state.range(0))+"_"+std::to_string(state.range(1))+"_"+std::to_string(state.range(2))+"_"+".txt";
        label_file.open(file.c_str());

        for (int i = 0; i < state.range(1); i++) {
            label_file << v105.labels[i] << "\n";
        }
        label_file.close();

        v105.clean();
        state.ResumeTiming();
    }
}

static void BM_Kmeans_GOLDSWICH(benchmark::State& state) { 

    for (auto _ : state) { 
        state.PauseTiming();
        Kmeans_bench v11(state.range(0), state.range(1), state.range(2));
        v11.init(argc_, argv_); 
        state.ResumeTiming();

        v11.kmeans_v11();
        

        state.PauseTiming();
        
        std::ofstream label_file;
        std::string file = "/mnt/c/Users/kaspe/OneDrive/Skrivebord/Reps/scalable-kmeans/cpp/result_GOLDSWICH"+std::to_string(state.range(0))+"_"+std::to_string(state.range(1))+"_"+std::to_string(state.range(2))+"_"+".txt";
        label_file.open(file.c_str());

        for (int i = 0; i < state.range(1); i++) {
            label_file << v11.labels[i] << "\n";
        }
        label_file.close();

        v11.clean();
        state.ResumeTiming();
    }
}

static void BM_Kmeans_v106(benchmark::State& state) { 
    for (auto _ : state) { 
        state.PauseTiming();
        Kmeans_bench v106(state.range(0), state.range(1), state.range(2));
        v106.init(argc_, argv_); 
        state.ResumeTiming();

        v106.kmeans_v106();
        

        state.PauseTiming();
        
        std::ofstream label_file;
        std::string file = "/mnt/c/Users/kaspe/OneDrive/Skrivebord/Reps/scalable-kmeans/cpp/result_MARIGOLD106"+std::to_string(state.range(0))+"_"+std::to_string(state.range(1))+"_"+std::to_string(state.range(2))+"_"+".txt";
        label_file.open(file.c_str());

        for (int i = 0; i < state.range(1); i++) {
            label_file << v106.labels[i] << "\n";
        }
        label_file.close();

        v106.clean();
        state.ResumeTiming();
    }
}

static void BM_Kmeans_v12(benchmark::State& state) { 
    for (auto _ : state) { 
        state.PauseTiming();
        Kmeans_bench v12(state.range(0), state.range(1), state.range(2));
        v12.init(argc_, argv_); 
        state.ResumeTiming();

        v12.kmeans_v12();
        

        state.PauseTiming();
        
        std::ofstream label_file;
        std::string file = "/mnt/c/Users/kaspe/OneDrive/Skrivebord/Reps/scalable-kmeans/cpp/result_v12"+std::to_string(state.range(0))+"_"+std::to_string(state.range(1))+"_"+std::to_string(state.range(2))+"_"+".txt";
        label_file.open(file.c_str());

        for (int i = 0; i < state.range(1); i++) {
            label_file << v12.labels[i] << "\n";
        }
        label_file.close();

        v12.clean();
        state.ResumeTiming();
    }
}

static void BM_Kmeans_v121(benchmark::State& state) { 
    for (auto _ : state) { 
        state.PauseTiming();
        Kmeans_bench v121(state.range(0), state.range(1), state.range(2));
        v121.init(argc_, argv_); 
        state.ResumeTiming();

        v121.kmeans_v121();
        

        state.PauseTiming();
        
        std::ofstream label_file;
        std::string file = "/mnt/c/Users/kaspe/OneDrive/Skrivebord/Reps/scalable-kmeans/cpp/result_v121"+std::to_string(state.range(0))+"_"+std::to_string(state.range(1))+"_"+std::to_string(state.range(2))+"_"+".txt";
        label_file.open(file.c_str());

        for (int i = 0; i < state.range(1); i++) {
            label_file << v121.labels[i] << "\n";
        }
        label_file.close();

        v121.clean();
        state.ResumeTiming();
    }
}




void data_load(std::string data_file_name_, double D, double N, int K) {
    data_arr = (double **) malloc(N * sizeof(double *));
        for (int i = 0; i < N; i++) {
            data_arr[i] = (double *) malloc(D * sizeof(double));
        }

    // read data from file
        // put into "data_arr" array
        // TODO: DONE
        std::string data_file_name = data_file_name_;
        //std::string label_file_name = argv[2];

        //std::cout << "data file name: " << data_file_name << std::endl;
        //std::cout << "label file name: " << label_file_name << std::endl;


        std::ifstream data_file(data_file_name.c_str());
        int data_count = std::count(std::istreambuf_iterator<char>(data_file),
                                    std::istreambuf_iterator<char>(), '\n');
        if (data_count < N) {
            DEBUGPRINT("%d      %f", data_count, N);

            std::cout << "NOT ENOUGH DATA!!\n";
            exit(3);
        }

        //std::cout << "count of data in file: " << data_count << " N: " << N << std::endl;

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
}


void free_data() {
    free(data_arr);
}

static void BM_Load_Data_VB(benchmark::State& state) {
    for (auto _ : state) { 
        std::string name = "../Data/steinn_14_jun/processed/misfit_VB_"+ std::to_string(state.range(0)) + "_h5_dct.txt";
        free_data();
        data_load(name, state.range(1), state.range(2), state.range(3));
        std::cout << "#############VB############" << std::endl;
    }
}

static void BM_Load_Data_Bi5d(benchmark::State& state) {
    for (auto _ : state) { 
        std::string name = "../Data/steinn_14_jun/processed/misfit_Bi5d_"+ std::to_string(state.range(0)) + "_h5_dct.txt";
        free_data();
        data_load(name, state.range(1), state.range(2), state.range(3));
        std::cout << "#############Bi5d############" << std::endl;
    }
}

static void BM_Load_Data_Se3d(benchmark::State& state) {
    for (auto _ : state) { 
        std::string name = "../Data/steinn_14_jun/processed/misfit_Se3d_"+ std::to_string(state.range(0)) + "_h5_dct.txt";
        free_data();
        data_load(name, state.range(1), state.range(2), state.range(3));
        std::cout << "#############Se3d###########" << std::endl;
    }
}

static void BM_Load_Data_flake(benchmark::State& state) {
    for (auto _ : state) { 
        std::string name = "../Data/steinn_14_jun/processed/gr_flake_"+ std::to_string(state.range(0)) + "_h5_dct.txt";
        free_data();
        data_load(name, state.range(1), state.range(2), state.range(3));
        std::cout << "#############Flake############" << std::endl;
    }
}

void add_all_VB(void (*load)(benchmark::State&),int dct_len_max, int N, int K) {
     //const std::string VB_8 = "../Data/steinn_14_jun/processed/misfit_VB_8_h5_dct.txt";
    int D = 8*8;
    //data_load("../Data/steinn_14_jun/processed/misfit_VB_8_h5_dct.txt", D, N, K);
    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);

    D=16*16;

    BENCHMARK(*load)->Args({16, D, N, K})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);

    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);

    D=32*32;

    BENCHMARK(*load)->Args({32, D, N, K})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);

    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);

    D=64*64;
    BENCHMARK(*load)->Args({64, D, N, K})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);


    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);

    D = 128*128;
    BENCHMARK(*load)->Args({128, D, N, K})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);


    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);

    D = 256*256;
    BENCHMARK(*load)->Args({256, D, N, K})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);


    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);

    D = 512*512;
    BENCHMARK(*load)->Args({512, D, N, K})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);


    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    
    D = 1024*1024;
    BENCHMARK(*load)->Args({1024, D, N, K})->Iterations(1)->Repetitions(1)->Unit(benchmark::kMillisecond);

    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
}  

void add_all_Bi5d(void (*load)(benchmark::State&),int dct_len_max, int N, int K) {
    int D = 8*8;
    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);

    D=16*16;

    BENCHMARK(*load)->Args({16, D, N, K})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);

    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);

    D=32*32;

    BENCHMARK(*load)->Args({32, D, N, K})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);

    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);

    D=64*64;
    BENCHMARK(*load)->Args({64, D, N, K})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);


    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);

    D = 128*128;
    BENCHMARK(*load)->Args({128, D, N, K})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);


    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);

    D = 256*256;
    BENCHMARK(*load)->Args({256, D, N, K})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);


    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);

    D = 512*512;
    BENCHMARK(*load)->Args({512, D, N, K})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);


    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    
    D = 1024*1024;
    BENCHMARK(*load)->Args({1024, D, N, K})->Iterations(1)->Repetitions(1)->Unit(benchmark::kMillisecond);

    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
}

void add_all_Se3d(void (*load)(benchmark::State&),int dct_len_max, int N, int K) {
    int D = 8*8;
    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);

    D=16*16;

    BENCHMARK(*load)->Args({16, D, N, K})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);

    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);

    D=32*32;

    BENCHMARK(*load)->Args({32, D, N, K})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);

    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);

    D=64*64;
    BENCHMARK(*load)->Args({64, D, N, K})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);


    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);

    D = 128*128;
    BENCHMARK(*load)->Args({128, D, N, K})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);


    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);

    D = 256*256;
    BENCHMARK(*load)->Args({256, D, N, K})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);


    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);

    D = 512*512;
    BENCHMARK(*load)->Args({512, D, N, K})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);


    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    
    D = 1024*1024;
    BENCHMARK(*load)->Args({1024, D, N, K})->Iterations(1)->Repetitions(1)->Unit(benchmark::kMillisecond);

    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
}

void add_all_flake(void (*load)(benchmark::State&),int dct_len_max, int N, int K) {
     //const std::string VB_8 = "../Data/steinn_14_jun/processed/misfit_VB_8_h5_dct.txt";
    int D = 8*8;
    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);

    D=16*16;

    BENCHMARK(*load)->Args({16, D, N, K})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);

    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);

    D=32*32;

    BENCHMARK(*load)->Args({32, D, N, K})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);

    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);

    D=64*64;
    BENCHMARK(*load)->Args({64, D, N, K})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);


    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);

    D = 128*128;
    BENCHMARK(*load)->Args({128, D, N, K})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);


    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);

    D = 256*256;
    BENCHMARK(*load)->Args({256, D, N, K})->Repetitions(1)->Iterations(1)->Unit(benchmark::kMillisecond);


    BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(5);

}



int main(int argc, char **argv) {

    argv_ = argv;
    argc_ = argc;
    int D = 128*128;
    int N = 168;
    int K = 5; 

    const std::string VB = "../Data/steinn_14_jun/processed/misfit_VB_128_h5_dct.txt";
    const std::string FLAKE = "../Data/steinn_14_jun/processed/gr_flake_128_h5_dct.txt";
    const std::string Se3d = "../Data/steinn_14_jun/processed/misfit_Se3d_128_h5_dct.txt";
    const std::string Bi5d = "../Data/steinn_14_jun/processed/misfit_Bi5d_128_h5_dct.txt";

    data_load(VB, D, N, K);
    //add_all_Bi5d(&BM_Load_Data_Bi5d, 1024, N, K);
    //add_all_Se3d(&BM_Load_Data_Se3d, 1024, N, K);
    //
    //add_all_VB(&BM_Load_Data_VB, 1024, N, K);
    
    //N = 1911;
    //data_load(FLAKE, D, N, K);
    //add_all_flake(&BM_Load_Data_flake, 256, N, K);
    BENCHMARK(BM_Kmeans)->Args({D, N, 5})->Args({D, N, 10})->Args({D, N, 15})->Args({D, N, 20})->Args({D, N, 25})->Args({D, N, 30})->Args({D, N, 35})->Args({D, N, 40})->Unit(benchmark::kMillisecond)->Iterations(1);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, 5})->Args({D, N, 10})->Args({D, N, 15})->Args({D, N, 20})->Args({D, N, 25})->Args({D, N, 30})->Args({D, N, 35})->Args({D, N, 40})->Unit(benchmark::kMillisecond)->Iterations(1);
    BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, 5})->Args({D, N, 10})->Args({D, N, 15})->Args({D, N, 20})->Args({D, N, 25})->Args({D, N, 30})->Args({D, N, 35})->Args({D, N, 40})->Unit(benchmark::kMillisecond)->Iterations(1);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, 5})->Args({D, N, 10})->Args({D, N, 15})->Args({D, N, 20})->Args({D, N, 25})->Args({D, N, 30})->Args({D, N, 35})->Args({D, N, 40})->Unit(benchmark::kMillisecond)->Iterations(1);
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, 5})->Args({D, N, 10})->Args({D, N, 15})->Args({D, N, 20})->Args({D, N, 25})->Args({D, N, 30})->Args({D, N, 35})->Args({D, N, 40})->Unit(benchmark::kMillisecond)->Iterations(1);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, 5})->Args({D, N, 10})->Args({D, N, 15})->Args({D, N, 20})->Args({D, N, 25})->Args({D, N, 30})->Args({D, N, 35})->Args({D, N, 40})->Unit(benchmark::kMillisecond)->Iterations(1);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, 5})->Args({D, N, 10})->Args({D, N, 15})->Args({D, N, 20})->Args({D, N, 25})->Args({D, N, 30})->Args({D, N, 35})->Args({D, N, 40})->Unit(benchmark::kMillisecond)->Iterations(1);
    //BENCHMARK(BM_Kmeans_GOLDSWICH)->Args({D, N, 5})->Args({D, N, 10})->Args({D, N, 15})->Args({D, N, 20})->Args({D, N, 25})->Args({D, N, 30})->Args({D, N, 35})->Args({D, N, 40})->Unit(benchmark::kMillisecond)->Iterations(1);
    //BENCHMARK(BM_Kmeans_v106)->Args({D, N, 5})->Args({D, N, 10})->Args({D, N, 15})->Args({D, N, 20})->Args({D, N, 25})->Args({D, N, 30})->Args({D, N, 35})->Args({D, N, 40})->Unit(benchmark::kMillisecond)->Iterations(1);
    //BENCHMARK(BM_Kmeans_v12)->Args({D, N, 5})->Args({D, N, 10})->Args({D, N, 15})->Args({D, N, 20})->Args({D, N, 25})->Args({D, N, 30})->Args({D, N, 35})->Args({D, N, 40})->Unit(benchmark::kMillisecond)->Iterations(1);
    //BENCHMARK(BM_Kmeans_v121)->Args({D, N, 5})->Args({D, N, 10})->Args({D, N, 15})->Args({D, N, 20})->Args({D, N, 25})->Args({D, N, 30})->Args({D, N, 35})->Args({D, N, 40})->Unit(benchmark::kMillisecond)->Iterations(1);
    


    //BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(1);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(1);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(1);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(1);
    //BENCHMARK(BM_Kmeans_GOLDSWICH)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(1);

    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, 5})->Unit(benchmark::kMillisecond);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(1);
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond);
    //BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond);
    
    
    //BENCHMARK(BM_Kmeans_ElkanHamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond);
    //BENCHMARK(BM_Kmeans_HamerlyStepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond);
    
    //BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond)->Iterations(1);
    //BENCHMARK(BM_Load_Data_Bi5d)->Args({256, D, N, K})->Iterations(1)->Unit(benchmark::kMillisecond);
    
    




    //BENCHMARK(BM_Kmeans)->Iterations(20)->Args({D, N, K})->Unit(benchmark::kMillisecond);
    //BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Iterations(20)->Unit(benchmark::kMillisecond);
    //BENCHMARK(BM_Kmeans_hamerly)->Iterations(10)->Unit(benchmark::kMillisecond);
    //BENCHMARK(BM_Kmeans_ElkanNewTE)->Iterations(2000)->Unit(benchmark::kMillisecond);
    //BENCHMARK(BM_Kmeans_HamerlyNewTE)->Iterations(1)->Unit(benchmark::kMillisecond);
    
    //BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Iterations(20)->Unit(benchmark::kMillisecond);
    //BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Iterations(1)->Unit(benchmark::kMillisecond);


    
    //D = 8*8;
    //data_load(VB_8, D, N, K);
    /*for (int dct_len = 8; dct_len >= 1024; dct_len *= 2) {
        D = dct_len*dct_len;
        BENCHMARK(BM_Load_Data_VB)->Args({dct_len, D, N, K})->Iterations(1)->Unit(benchmark::kMillisecond);
        
        BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond);
        BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond);
        BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond);
        BENCHMARK(BM_Kmeans_ElkanNewTE)->Args({D, N, K})->Unit(benchmark::kMillisecond);
        BENCHMARK(BM_Kmeans_HamerlyNewTE)->Args({D, N, K})->Unit(benchmark::kMillisecond);
        BENCHMARK(BM_Kmeans_ElkanHamerlyNewTE)->Args({D, N, K})->Unit(benchmark::kMillisecond);
        BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond);
        BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond);
    }*/

    /*BENCHMARK(BM_Load_Data_VB)->Args({8, D, N, K})->Iterations(1)->Unit(benchmark::kMillisecond);
        
        BENCHMARK(BM_Kmeans)->Args({D, N, K})->Unit(benchmark::kMillisecond);
        BENCHMARK(BM_Kmeans_elkan)->Args({D, N, K})->Unit(benchmark::kMillisecond);
        BENCHMARK(BM_Kmeans_hamerly)->Args({D, N, K})->Unit(benchmark::kMillisecond);
        BENCHMARK(BM_Kmeans_ElkanNewTE)->Args({D, N, K})->Unit(benchmark::kMillisecond);
        BENCHMARK(BM_Kmeans_HamerlyNewTE)->Args({D, N, K})->Unit(benchmark::kMillisecond);
        BENCHMARK(BM_Kmeans_ElkanHamerlyNewTE)->Args({D, N, K})->Unit(benchmark::kMillisecond);
        BENCHMARK(BM_Kmeans_stepwise)->Args({D, N, K})->Unit(benchmark::kMillisecond);
        BENCHMARK(BM_Kmeans_MARIGOLD)->Args({D, N, K})->Unit(benchmark::kMillisecond);
    */

    benchmark::Initialize(&argc, argv);
    //std::cout << "test" << std::endl;
    benchmark::RunSpecifiedBenchmarks();
    //std::cout << "test" << std::endl;

    free_data();

}




