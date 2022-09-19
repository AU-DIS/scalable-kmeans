#Importing required modules
import math
import numpy as np

import kmeans_common_func 
#In case of need: Freshly reload the import modules
import importlib
importlib.reload(kmeans_common_func)


def dist_masked(x, centroids, x_squared, centroid_squared, is_candidate, labels, level, incremental_dots):
    ## x_squared.shape === N, log(sqrt(d))
    ## centroid_squared.shape == k, log(sqrt(d))
    ## dist returned .shape == N, k
    ## is_candidate.shape == N, k
    ## labels.shape == 1, N ---> just to check whether a cluster has already been assigned or not
    ## only the values for not-assigned and candidates is correct, the rest is
    
    lb_dist = [[0 for _ in range(len(centroids))] for _ in range(len(x))]
    ub_dist = [[0 for _ in range(len(centroids))] for _ in range(len(x))]
    
    d_sqrt = int(math.sqrt(len(x[0])))
    two_p_lv_m1 = int(2**(level - 1))
    two_p_lv = int(2**level)
    
    
    for i in range(len(x_squared)):
        if labels[i] >= 0: continue
        for j in range(len(centroid_squared)):
            if is_candidate[i][j] == 0: 
                lb_dist[i][j] = float('inf')
                ub_dist[i][j] = float('inf')
                continue
            ## the known parts
            
            
            ## sum of squareds
            lb_dist[i][j] += x_squared[i][level]
            ub_dist[i][j] += x_squared[i][level]
            
            lb_dist[i][j] += centroid_squared[j][level]
            ub_dist[i][j] += centroid_squared[j][level]
            
            ## dot product
            ## Known part
            if level == 0:
                incremental_dots[i][j] = (x[i][0] * centroids[j][0])
            else:
                for l in range(two_p_lv_m1):
                    incremental_dots[i][j] += np.dot(x[i][d_sqrt*(two_p_lv_m1 + l): d_sqrt * (two_p_lv_m1 + l) + two_p_lv],  centroids[j][d_sqrt*(two_p_lv_m1 + l): d_sqrt * (two_p_lv_m1 + l) + two_p_lv])
                    incremental_dots[i][j] += np.dot(x[i][l*d_sqrt + two_p_lv_m1:l*d_sqrt + two_p_lv],  centroids[j][l*d_sqrt + two_p_lv_m1:l*d_sqrt + two_p_lv])
                                                      
            lb_dist[i][j] -= 2*incremental_dots[i][j]
            ub_dist[i][j] -= 2*incremental_dots[i][j]
            
            
#           DEBUGGING
#           if i == 0:
#               print('known')
#               print('x_sq: ' + str(x_squared[i][level]) + ' c_sq:  ' + str(centroid_squared[j][level]) + ' this_dot: ' + str(this_dot))
            
            
            ## unknown parts
            ## sum of squareds
            lb_dist[i][j] += (x_squared[i][-1] - x_squared[i][level])
            ub_dist[i][j] += (x_squared[i][-1] - x_squared[i][level])
            
            lb_dist[i][j] += (centroid_squared[j][-1] - centroid_squared[j][level])
            ub_dist[i][j] += (centroid_squared[j][-1] - centroid_squared[j][level])
            
            ## dot product approx
            this_dot = math.sqrt((x_squared[i][-1] - x_squared[i][level]) * (centroid_squared[j][-1] - centroid_squared[j][level]))
            lb_dist[i][j] -= 2 * this_dot
            ub_dist[i][j] += 2 * this_dot
            
            if lb_dist[i][j] < 0: lb_dist[i][j] = 0
            if ub_dist[i][j] < 0: ub_dist[i][j] = 0
            lb_dist[i][j] = math.sqrt(lb_dist[i][j])
            ub_dist[i][j] = math.sqrt(ub_dist[i][j])

#           DEBUGGING
#           if i == 0:
#               print('unknown')
#               print('x_sq: ' + str((x_squared[i][-1] - x_squared[i][level])) + ' c_sq:  ' + str((centroid_squared[j][-1] - centroid_squared[j][level])) + ' this_dot: ' + str(this_dot))
    
    return lb_dist, ub_dist, incremental_dots


def get_labels(x, k, centroids, x_squared, centroid_squared, mask, is_print = 'no'):
    ## init labels as -1, then assign to right labels in the subsequent codes 
    labels = np.array([-1 for _ in range(len(x))])
    level = 0

    ## shows if that centroid is a candidate for being the closest to that point; 1 if true, 0 is false
    #is_candidate = np.ones((len(x), k)) #When there is no mask for triangled version
    is_candidate = mask
    incremental_dots = [[0 for _ in range(k)] for _ in range(len(x))]

    lb_return_dists, ub_return_dists = [], []

    cal_count_itr = 0

    while(level < len(x_squared[0])):
        cal_count_each_level = 0
        
        these_dists_lb, these_dists_ub, incremental_dots = dist_masked(x, centroids, x_squared, centroid_squared, is_candidate, labels, level, incremental_dots)
        # have to do this before pruning, aka updating is_candidate
        if level == 0:
            lb_return_dists, ub_return_dists = these_dists_lb.copy(), these_dists_ub.copy()
            cal_count_each_level += pow(2, (2 * level)) * sum(sum(is_candidate)) #Number of features times number of candidate distances
            #Fatemeh's way of number of features: 2^(2*level) - 2^(2*(level - 1))

        else:
            cal_count_each_level += (pow(2, (2 * level)) - pow(2, (2* (level-1)))) * sum(sum(is_candidate)) #Number of features times number of candidate distances
            for i in range(len(x)):
                for j in range(k):
                    if is_candidate[i][j] == 1:
                        lb_return_dists[i][j] = these_dists_lb[i][j]
                        ub_return_dists[i][j] = these_dists_ub[i][j]
                        #cal_count_each_level += 1

                        
        smallest_ub = [0 for _ in range(len(x))]
        for i in range(len(x)):
            for j in range(k):
                if these_dists_ub[i][j] < these_dists_ub[i][smallest_ub[i]]: smallest_ub[i] = j
         
        for i in range(len(x)):
            if labels[i] >= 0: continue

            if level == (len(x_squared[0]) - 1): labels[i]= smallest_ub[i]

            else:
                for j in range(k):
                    if j == smallest_ub[i]: continue
                    if is_candidate[i][j] == 0: continue

                    if these_dists_lb[i][j] >= these_dists_ub[i][smallest_ub[i]]:
                        ## def not a candidate
                        is_candidate[i][j] = 0
                if sum(is_candidate[i]) == 1:
                    ## then the only one left is the one with the smallest_ub
                    #print('assigned label ' + str(smallest_ub[i]))
                    labels[i] = smallest_ub[i]
        
        level +=1
        if is_print == 'yes': print('Stepwise number of features for cal_count_each_level for level ' + str(level) + ' is: ' + str(cal_count_each_level))
        cal_count_itr += cal_count_each_level
    
    if print == 'yes': print('Stepwise cal_count_itr for this iteration: ' + str(cal_count_itr))
    if print == 'yes': print('Traditional total stepwise count for this iteration: ' + str(len(x)*len(x_squared[0])*k))
    return labels, lb_return_dists, ub_return_dists, cal_count_itr

def get_square_squared_sums(data):
    d_sqrt = int(math.sqrt(len(data[0])))
    log_d_sqrt = int(math.log(d_sqrt, 2))
    
    res = [[0 for _ in range(log_d_sqrt + 1)] for _ in range(len(data))]
    
    for i in range(len(data)):
        res[i][0] = (data[i][0] ** 2)
        for level in range(1, log_d_sqrt+1):
            if level > 0: res[i][level] += res[i][level-1]
            two_p_level_m1 = int(2**(level-1))
            two_p_level = int(2**level)
            for l in range(two_p_level_m1):
                res[i][level] += sum([folan**2 for folan in data[i][(d_sqrt*(two_p_level_m1 + l)): (d_sqrt*(two_p_level_m1 + l)) + two_p_level]])
                res[i][level] += sum([folan**2 for folan in data[i][l*d_sqrt + two_p_level_m1: l*d_sqrt + two_p_level]])
    return res





class kmeans_class:
    def __init__(self, data_points, k, max_iter=100, initialCentroids=None, kmeans_threshold=0, is_print='no'):
        # k: no. of clusters
        # data_points: data to cluster
        # max_iter: iteration before termination
        # kmeans_threshold: termination condition to detect convergence if the centroids do not move more than the 'kmeans_threshold' value
        
        self.k = k
        self.data_points = data_points
        self.max_iter = max_iter
        self.kmeans_threshold = kmeans_threshold
        self.method = 'MARIGOLD'

        if initialCentroids != None:
            self.centroids = dict(zip(list(range(self.k)),initialCentroids))
        else:
            self.centroids = dict(zip(list(range(self.k)), self.selectSeeds(self.k)))

        self.execute_clustering()
    

    def selectSeeds(self, k):
        #Not in use for random function below
        #random.seed(1)
        #seeds = random.sample(self.pointList, k)
        
        #For same initialization with other methods
        idx = [i for i in range(k)]
        seeds = [[] for _ in range(k)]
        for i in range(k):
            seeds[i] = self.data_points[idx[i]]
        return seeds


    #main function for clustering
    def execute_clustering(self):
        
        data_x = np.array(self.data_points)

        self.labels = [0 for point in data_x]
        iter_count = 1
        if print == 'yes': print('Current iteration: ', iter_count)

        ## The proposed MARIGOLD
        if self.method == 'MARIGOLD':
            
            self.classifications = {} ## Points coordinates by centroid
            self.pointsClassif = {} ## Points indices by centroid
            
            lowerBounds = np.zeros((data_x.shape[0],self.k)) ## Lower bounds matrix. Dimensions : (nb_points, nb_centroids)
            upperBounds = np.zeros((data_x.shape[0])) ## Upper bounds vector : Dimension : (nb_points)
            lowerBounds_Hamerly = np.zeros((data_x.shape[0])) ## Hamerly's lower bounds vector : Dimension : (nb_points)


            for i in range(self.k):
                self.classifications[i] = []
                self.pointsClassif[i] = []
            
            mask = np.ones((data_x.shape[0],self.k)) #Setting triangle inequality filter for distance candidates
            x_squared = get_square_squared_sums(data_x)
            centroid_squared = get_square_squared_sums(self.centroids)
            labels, lb_dists, ub_dists, cal_count_itr = get_labels(data_x, self.k, self.centroids, x_squared, centroid_squared, mask)

                
            i=0
            for point in data_x:
                #distances = [np.linalg.norm(point-self.centroids[centroid]) for centroid in self.centroids]
                classification = labels[i]
                    
                ## Centroid assigned to each point
                self.classifications[classification].append(point)
                self.pointsClassif[classification].append(i)
                i+=1
                
            prevCentroids = dict(self.centroids.copy())
            prevClassifications = dict(self.classifications.copy())
            prevPointsClassif = dict(self.pointsClassif.copy())

            for classification in self.classifications:
                if len(self.classifications[classification]) == 0:
                    pass

                else:
                    self.centroids[classification] = np.average(self.classifications[classification],axis=0)

            optimized = [True for centroid in self.centroids]
            
            centroidDistanceChange = {}

            for centroid in self.centroids:
                original_centroid = prevCentroids[centroid]
                current_centroid = self.centroids[centroid]
                centroidDistanceChange[centroid] = np.linalg.norm(original_centroid-current_centroid)
                
                #if abs(np.sum((current_centroid-original_centroid)/original_centroid*100.0)) > self.kmeans_threshold :
                if abs(np.sum(np.divide((current_centroid-original_centroid), original_centroid, out=np.zeros_like(current_centroid), where=original_centroid!=0)*100.0)) > self.kmeans_threshold :
                    optimized[centroid] = False

            if False not in optimized:
                
                for centroid in self.pointsClassif:
                    for point in self.pointsClassif[centroid]:
                        self.labels[point] = centroid
                return
            

            ## Lower bound distance between the point and each center
            ## Initialized as distance between the point and each initial centroid
            lowerBounds = lb_dists
            i=0
            for point in data_x:
                # For Hamerly's second lower bound
                ## Second Lower bound distance between the point and the second closest center
                ## Initialized as distance between the point and second closest initial centroid
                if len(set(lowerBounds[i])) == 1:
                    lowerBounds_Hamerly[i] = lowerBounds[i][0]
                else:
                    lowerBounds_Hamerly[i] = sorted(set(lowerBounds[i]))[1]
                #lowerBounds_Hamerly[i] = sorted(set([item for item in lowerBounds[i] if item >= lowerBounds[i][labels[i]]]))[1]

                ## Upper bound distance between the point and assigned centroid
                upperBounds[i] = ub_dists[i][labels[i]]                    
                i+=1


            maxCentroidDistanceChange = max(centroidDistanceChange.values())
            ## Update lower and upper bound distances
            for centroid in self.pointsClassif:
                for i in list(range(data_x.shape[0])):
                    lowerBounds[i][centroid] -= centroidDistanceChange[centroid]
                    lowerBounds[i][centroid] -= maxCentroidDistanceChange
                    
                for i in self.pointsClassif[centroid]:
                    upperBounds[i] += maxCentroidDistanceChange
                    #upperBounds[i] += centroidDistanceChange[centroid]
                    furthestMovingCentroid = max(centroidDistanceChange, key=centroidDistanceChange.get)
                    if furthestMovingCentroid == centroid:
                        #lowerBounds_Hamerly[i] -= sorted(set(centroidDistanceChange.values()), reverse=True)[1]
                        lowerBounds_Hamerly[i] -= maxCentroidDistanceChange
                    else:
                        lowerBounds_Hamerly[i] -= maxCentroidDistanceChange

                    
            #For ecokmeans efficiency measure
            dist_count = data_x.shape[0]*self.k #Already done above
            itr_discard_count_total = 0 #See how many dist calculation are discarded
            cal_count_total = 0 #Total number of features passed through distance calculation
            cal_count_total += cal_count_itr #cal_count_itr stands for number of data features passed through distance calc module in each iteration 

            ## Repeat until convergence
            for it in range(self.max_iter):
                iter_count += 1
                if print == 'yes': print('Current iteration: ', iter_count)

                listCentroids = list(range(self.k))
                self.classifications = {} ## Points coordinates by centroid
                self.pointsClassif = {} ## Points indices by centroid
                centroidDistances = {} ## Distances to other centroids for each centroid
                closestCentroidDistances = {} ## Distance to closest centroid for each centroid
                
                for i in range(self.k):
                    self.classifications[i] = []
                    self.pointsClassif[i] = []
                    centroidDistances[i] = [np.linalg.norm(self.centroids[i]-self.centroids[c_prime]) for c_prime in self.centroids]
                    closestCentroidDistances[i] = min(centroidDistances[i][:i]+centroidDistances[i][i+1:])
                
                #Preparing mask to discard calcualations using triangle inequality
                mask = np.ones((data_x.shape[0],self.k)) #Setting triangle inequality filter for distance candidates

                dist_count += data_x.shape[0]*self.k
                this_itr_dist_discard = 0 #for counting number of discarded dist calc in each inner iteration of kmeans

                for centroid in prevPointsClassif:
                    for i in prevPointsClassif[centroid]:
                    #for i in range(data_x.shape[0]):
                        
                        #Elkan Lemma1 
                        ## Check if upper bound lower than 1/2 of distance with closest centroid
                        hamerly_bound = max(0.5*closestCentroidDistances[centroid], lowerBounds_Hamerly[i])
                        if upperBounds[i] <= hamerly_bound: #For Elkan lemma1 and Hamerly
                        #if upperBounds[i] <= 0.5*closestCentroidDistances[centroid]: #for Elkan lemma1
                            ## If condition is met : said point keeps its centroid with no further computation needed
                            #mask[i] = 0
                            #mask[i][centroid] = 1
                            #this_itr_dist_discard += 1
                            #itr_discard_count_total += 1
                            pass
                        
                        else:
                            assigned_centroid = centroid
                            for c_prime in (listCentroids[:centroid]+listCentroids[centroid+1:]):
                                #Proposed New_Lemma (tighter version of Elkan lemma1 for individual distances)
                                if lowerBounds[i][c_prime] > (upperBounds[i] + centroidDistanceChange[assigned_centroid] + centroidDistanceChange[c_prime]):
                                    ## If condition is met : said point keeps its centroid with no further computation needed
                                    this_itr_dist_discard += 1
                                    itr_discard_count_total += 1
                                    mask[i][c_prime] = 0

                                #Elkan lemma2 covering Hammerly lemma partially
                                ## Check if lower bound between point and c_prime < upper bound between point and its current centroid
                                ## AND if (0.5*distance between current centroid and c_prime) < upper bound between point and its current centroid
                                elif ((upperBounds[i] < lowerBounds[i][c_prime]) or (upperBounds[i] < 0.5*centroidDistances[assigned_centroid][c_prime])): 
                                    this_itr_dist_discard += 1
                                    itr_discard_count_total += 1
                                    mask[i][c_prime] = 0
                                    pass
                                else:
                                    #Distance calculation necessary
                                    pass
                        
                            #self.classifications[assigned_centroid].append(data_x[i])
                            #self.pointsClassif[assigned_centroid].append(i)

                centroid_squared = get_square_squared_sums(self.centroids)
                labels, lb_dists, ub_dists, cal_count_itr = get_labels(data_x, self.k, self.centroids, x_squared, centroid_squared, mask)

                cal_count_total += cal_count_itr #For performance measurement

                
                i=0
                for point in data_x:
                    #distances = [np.linalg.norm(point-self.centroids[centroid]) for centroid in self.centroids]
                    classification = labels[i]
                    
                    ## Centroid assigned to each point
                    self.classifications[classification].append(point)
                    self.pointsClassif[classification].append(i)
                    i += 1

                prevCentroids = dict(self.centroids.copy())
                prevClassifications = dict(self.classifications.copy())
                prevPointsClassif = dict(self.pointsClassif.copy())
                
                for classification in self.classifications:
                    if len(self.classifications[classification]) == 0:
                        pass

                    else:
                        self.centroids[classification] = np.average(self.classifications[classification],axis=0)

                optimized = [True for centroid in self.centroids]
                
                centroidDistanceChange = {}

                for centroid in self.centroids:
                    original_centroid = prevCentroids[centroid]
                    current_centroid = self.centroids[centroid]
                    centroidDistanceChange[centroid] = np.linalg.norm(original_centroid-current_centroid)

                    #if abs(np.sum((current_centroid-original_centroid)/original_centroid*100.0)) > self.kmeans_threshold :
                    if abs(np.sum(np.divide((current_centroid-original_centroid), original_centroid, out=np.zeros_like(current_centroid), where=original_centroid!=0)*100.0)) > self.kmeans_threshold :
                        optimized[centroid] = False

                if False not in optimized:
                    break
                    

                i=0
                for point in data_x:
                    ## Lower bound distance between the point and each center
                    ## Initialized as distance between the point and each initial centroid
                    for centroid in range(self.k):
                        #if (mask[i][centroid] == 1) and (sum(mask[i]!=1)): 
                        if mask[i][centroid] == 1: 
                            lowerBounds[i][centroid] = lb_dists[i][centroid]
                    i+=1

                i=0
                for point in data_x:
                    # For Hamerly's second lower bound
                    ## Second Lower bound distance between the point and the second closest center
                    ## Initialized as distance between the point and second closest initial centroid
                    lowerBounds_Hamerly[i] = sorted(set(lowerBounds[i]))[1]
                    #lowerBounds_Hamerly[i] = sorted(set([item for item in lowerBounds[i] if item >= lowerBounds[i][labels[i]]]))[1]

                    ## Upper bound distance between the point and assigned centroid
                    upperBounds[i] = ub_dists[i][labels[i]]                    
                    i+=1


                maxCentroidDistanceChange = max(centroidDistanceChange.values())
                ## Update lower and upper bound distances
                for centroid in self.pointsClassif:
                    for i in list(range(data_x.shape[0])):
                        #lowerBounds[i][centroid] -= centroidDistanceChange[centroid]
                        lowerBounds[i][centroid] -= maxCentroidDistanceChange
                    
                    for i in self.pointsClassif[centroid]:
                        upperBounds[i] += maxCentroidDistanceChange
                        #upperBounds[i] += centroidDistanceChange[centroid]
                        furthestMovingCentroid = max(centroidDistanceChange, key=centroidDistanceChange.get)
                        if furthestMovingCentroid == centroid:
                            #lowerBounds_Hamerly[i] -= sorted(set(centroidDistanceChange.values()), reverse=True)[1]
                            lowerBounds_Hamerly[i] -= maxCentroidDistanceChange
                        else:
                            lowerBounds_Hamerly[i] -= maxCentroidDistanceChange


        
        ## Update labels (cluster) for each point
        for centroid in self.pointsClassif:
            for point in self.pointsClassif[centroid]:
                self.labels[point] = centroid
        
        labels = self.labels
        centroids = self.centroids
        print('converged in ' + str(iter_count) + " iterations")
        print('Distance calc access completely discarded through trangle inequality: ', (dist_count-itr_discard_count_total))

        datasize = data_x.shape[1]
        total_lloyd_feature_count = (data_x.shape[0] * self.k) * iter_count * datasize
        total_stepwise_feature_count = cal_count_total
        print('Traditional Lloyd total features count for all iterations: ' + str(total_lloyd_feature_count))
        print('Stepwise and Triangle Inequality based feature count in distance calculation: ' + str(total_stepwise_feature_count))
        print('Percentage of discarded features by the stepwise method for all iterations: ' + str(100 - ((total_stepwise_feature_count/total_lloyd_feature_count)*100))+ '%')
                

        return labels, centroids
                        

def kmeans_eco_manual(x, k, max_iter=100, initialCentroids=None, kmeans_threshold=0, is_print='no'):
    results = kmeans_class(x, k, max_iter=100, initialCentroids=None, kmeans_threshold=0)
    labels = results.labels
    centroids = results.centroids
    return labels, centroids


