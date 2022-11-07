#Importing required modules
import math
import numpy as np

import kmeans_common_func 
#In case of need: Freshly reload the import modules
import importlib
importlib.reload(kmeans_common_func)



###########################################################
#Done! def getSquareSquaredSums(data)
#Done! labels[point], lower_bounds[point], upper_bounds[point] = setLabelMG(point, centroids, x_squared[point], centroid_squared, step_level, labels[point], lower_bounds[point], upper_bounds[point])
#Done! converged, centroid_divergence = kmeansRecalculate(data_x, centroids, labels)
#Done! lower_bounds, upper_bounds, lower_bounds_hamerly, near_centroid_distances = updateBounds(data_x, centroids, lower_bounds, upper_bounds, lower_bounds_hamerly, centroid_divergence, near_centroid_distances)   
#Done! step_lb, step_ub = distToLevel(point, centroids[c_current], x_squared[c_current], centroid_squared[c_current], level, step_level)


#def getSquareSquaredSums(data)
def getSquareSquaredSums(data):
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



#step_lb, step_ub, incremental_dots[c_current] = distToLevel(point, centroids[c_current], x_squared[c_current], centroid_squared[c_current], level, step_level, idots)
def distToLevel(point, centroid, x_squared, centroid_squared, level, step_level, idots, cal_count_each_level):
    d_sqrt = int(math.sqrt(point.shape[0]))
    two_p_lv_m1 = int(2**(level - 1))
    two_p_lv = int(2**level)
    
    
    ## sum of squareds
    lb_dist = x_squared[-1] + centroid_squared[-1]
    ub_dist = x_squared[-1] + centroid_squared[-1]


    ########## the known parts
    ## sum of squareds
    #lb_dist += x_squared[level] + centroid_squared[level]
    #ub_dist += x_squared[level] + centroid_squared[level]
            
    ## dot product
    ## Known part
    if level == 0:
        idots = (point[0] * centroid[0])
        cal_count_each_level += 1 #pow(2, (2 * level))
    else:
        cal_count_each_level += (pow(2, (2 * level)) - pow(2, (2* (level-1))))
        for l in range(two_p_lv_m1):
            idots += np.dot(point[d_sqrt*(two_p_lv_m1 + l): d_sqrt * (two_p_lv_m1 + l) + two_p_lv],  centroid[d_sqrt*(two_p_lv_m1 + l): d_sqrt * (two_p_lv_m1 + l) + two_p_lv])
            idots += np.dot(point[l*d_sqrt + two_p_lv_m1:l*d_sqrt + two_p_lv],  centroid[l*d_sqrt + two_p_lv_m1:l*d_sqrt + two_p_lv])
                                                      
    lb_dist -= 2*idots
    ub_dist -= 2*idots
            
            
    ## unknown parts
    ## sum of squareds
    #lb_dist += (x_squared[-1] - x_squared[level])
    #ub_dist += (x_squared[-1] - x_squared[level])
            
    #lb_dist += (centroid_squared[-1] - centroid_squared[level])
    #ub_dist += (centroid_squared[-1] - centroid_squared[level])
            
    ## dot product approx
    this_dot = math.sqrt((x_squared[-1] - x_squared[level]) * (centroid_squared[-1] - centroid_squared[level]))
    lb_dist -= 2 * this_dot
    ub_dist += 2 * this_dot
            
    if lb_dist < 0: lb_dist = 0
    if ub_dist < 0: ub_dist = 0
    lb_dist = math.sqrt(lb_dist)
    ub_dist = math.sqrt(ub_dist)

    return lb_dist, ub_dist, idots, cal_count_each_level
    


#converged, centroid_divergence = kmeansRecalculate(data_x, centroids, labels)
def kmeansRecalculate(data_x, centroids, labels, kmeans_threshold=0):
    centroids_old = centroids.copy()
    converged = True
    centroid_divergence = np.zeros((centroids.shape[0]))

    classifications = {}
    for i in range(centroids.shape[0]):
        classifications[i] = []
    for point in range(data_x.shape[0]):
        classifications[labels[point]].append(data_x[point])
    for c_current in range(centroids.shape[0]):
        centroids[c_current] = np.average(classifications[c_current],axis=0)
        centroid_divergence[c_current] = np.linalg.norm(centroids[c_current]-centroids_old[c_current])
        if centroid_divergence[c_current] > kmeans_threshold: converged = False

    return converged, centroids



#labels[point], cal_count_each_point = setLabelMG(data_x[point], centroids, x_squared[point], centroid_squared, step_level)
def setLabelMG(point, centroids, x_squared, centroid_squared, step_level):
    level = 0
    lower_bounds = np.zeros((centroids.shape[0])) ## Lower bounds matrix. Dimensions : (nb_points, nb_centroids)
    upper_bound = np.inf ## Upper bounds vector : Dimension : (nb_points)
    mask = np.ones((centroids.shape[0]))

    incremental_dots = [0 for _ in range(centroids.shape[0])]
    cal_count_each_point = 0

    while (level <= step_level) and (sum(mask) > 1):
        for c_current in range(centroids.shape[0]):
            if (mask[c_current] != 1): continue
        
            if upper_bound < lower_bounds[c_current]:
                mask[c_current] = 0
            else:
                step_lb, step_ub, incremental_dots[c_current], cal_count_each_point = distToLevel(point, centroids[c_current], x_squared, centroid_squared[c_current], level, step_level, incremental_dots[c_current], cal_count_each_point)
                #if step_lb > lower_bounds[c_current]: 
                lower_bounds[c_current] = step_lb #Keep highest lb per c
                if step_ub < upper_bound: 
                    label = c_current
                    upper_bound = step_ub

        level += 1
    
    return label, cal_count_each_point 





class kmeans_class:
    def __init__(self, data_points, k, max_iter=100, initialCentroids=None, kmeans_threshold=0, is_print='no', mode='fast'):
        # k: no. of clusters
        # data_points: data to cluster
        # max_iter: iteration before termination
        # kmeans_threshold: termination condition to detect convergence if the centroids do not move more than the 'kmeans_threshold' value
        
        self.k = k
        self.data_points = data_points
        self.max_iter = max_iter
        self.kmeans_threshold = kmeans_threshold
        self.method = 'stepwise'
        self.mode = mode #Based on distance mode

        if initialCentroids != None:
            #self.centroids = dict(zip(list(range(self.k)),initialCentroids))
            self.centroids = initialCentroids
        else:
            #self.centroids = dict(zip(list(range(self.k)), self.selectSeeds(self.k)))
            self.centroids = self.selectSeeds(self.k)

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
        centroids = np.array(self.centroids)

        self.labels = [0 for point in data_x]
        self.classifications = {} ## Points coordinates by centroid
        self.pointsClassif = {} ## Points indices by centroid
            
        x_squared = getSquareSquaredSums(data_x)
        step_level = int(math.log(data_x.shape[1], 4)) ##'L'
        #step_level = int(len(x_squared[0])) ##'L'

        labels = [0 for point in data_x]

        #For stepwise efficiency measure
        cal_count_total = 0
        iter_count = 0
        if print == 'yes': print('Current iteration: ', iter_count)

        ## Repeat until convergence
        for _ in range(self.max_iter):
            iter_count += 1
            if print == 'yes': print('Current iteration: ', iter_count)

            centroid_squared = getSquareSquaredSums(centroids)
            
            for point in range(data_x.shape[0]): 
                labels[point], cal_count_each_point = setLabelMG(data_x[point], centroids, x_squared[point], centroid_squared, step_level)
                cal_count_total += cal_count_each_point

                
            converged, centroids = kmeansRecalculate(data_x, centroids, labels) 

            if converged: 
                break    

        print('Total iterations taken to converge: ', iter_count)
        print('Total stepwise feature count: ', cal_count_total)
        self.labels = labels
        self.centroids = centroids
        return labels, centroids        
                        

def kmeans_stepwise_manual(x, k, max_iter=100, initialCentroids=None, kmeans_threshold=0, is_print='no', mode='fast'):
    results = kmeans_class(x, k, max_iter=100, initialCentroids=None, kmeans_threshold=0, mode=mode)
    labels = results.labels
    centroids = results.centroids
    return labels, centroids

