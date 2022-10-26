import numpy as np
#Importing required modules
import math

import kmeans_common_func 
#In case of need: Freshly reload the import modules
import importlib
importlib.reload(kmeans_common_func)

class kmeans_class:
    def __init__(self, data_points, k, max_iter=100, initialCentroids=None, kmeans_threshold=0, print='no', mode='fast'):
        # k: no. of clusters
        # data_points: data to cluster
        # max_iter: iteration before termination
        # kmeans_threshold: termination condition to detect convergence if the centroids do not move more than the 'kmeans_threshold' value
        
        self.k = k
        self.data_points = data_points
        self.max_iter = max_iter
        self.kmeans_threshold = kmeans_threshold
        self.method = 'Lloyd'
        self.mode = mode #Based on distance mode

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

        ## Lloyd's Kmeans
        if self.method == 'Lloyd':
            
            self.classifications = {} ## Points coordinates by centroid
            self.pointsClassif = {} ## Points indices by centroid
            
            for i in range(self.k):
                self.classifications[i] = []
                self.pointsClassif[i] = []
                
            x_sq_sum = np.zeros((data_x.shape[0]))
            c_sq_sum = np.zeros((self.k))
            if self.mode == 'semi_optimized': x_sq_sum = kmeans_common_func.squared_sum(data_x) #Based on distance mode
            if self.mode == 'semi_optimized': c_sq_sum = kmeans_common_func.squared_sum(self.centroids)

            i=0
            for point in data_x:
                #distances = [np.linalg.norm(point-self.centroids[centroid]) for centroid in self.centroids]
                #distances = [kmeans_common_func.euclidean_distance_raw(point, self.centroids[centroid]) for centroid in self.centroids]
                distances = [kmeans_common_func.euclidean_distance(point, self.centroids[centroid], x_sq_sum[i], c_sq_sum[centroid], mode=self.mode) for centroid in self.centroids]
                classification = distances.index(min(distances))
                    
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
            if self.mode == 'semi_optimized': c_sq_sum = kmeans_common_func.squared_sum(self.centroids)

            optimized = [True for centroid in self.centroids]
            
            centroidDistanceChange = {}

            for centroid in self.centroids:
                original_centroid = prevCentroids[centroid]
                current_centroid = self.centroids[centroid]
                centroidDistanceChange[centroid] = kmeans_common_func.euclidean_distance(original_centroid, current_centroid)
                #centroidDistanceChange[centroid] = np.linalg.norm(original_centroid-current_centroid)
                
                #if abs(np.sum((current_centroid-original_centroid)/original_centroid*100.0)) > self.kmeans_threshold :
                if abs(np.sum(np.divide((current_centroid-original_centroid), original_centroid, out=np.zeros_like(current_centroid), where=original_centroid!=0)*100.0)) > self.kmeans_threshold :
                    optimized[centroid] = False

            if False not in optimized:
                
                for centroid in self.pointsClassif:
                    for point in self.pointsClassif[centroid]:
                        self.labels[point] = centroid
                return
            
                    
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
                    #centroidDistances[i] = [np.linalg.norm(self.centroids[i]-self.centroids[c_prime]) for c_prime in self.centroids]
                    centroidDistances[i] = [kmeans_common_func.euclidean_distance(self.centroids[i], self.centroids[c_prime]) for c_prime in self.centroids]
                    closestCentroidDistances[i] = min(centroidDistances[i][:i]+centroidDistances[i][i+1:])
                
                for centroid in prevPointsClassif:
                    for i in prevPointsClassif[centroid]:
                        
                            r = True
                            assigned_centroid = centroid
                            for c_prime in (listCentroids[:centroid]+listCentroids[centroid+1:]):
                                    if r:
                                        #distToCurrentCentroid = np.linalg.norm(data_x[i] - self.centroids[assigned_centroid])
                                        #distToCurrentCentroid = kmeans_common_func.euclidean_distance(data_x[i], self.centroids[assigned_centroid])
                                        distToCurrentCentroid = kmeans_common_func.euclidean_distance(data_x[i], self.centroids[assigned_centroid], x_sq_sum[i], c_sq_sum[assigned_centroid], mode=self.mode)
                                        r = False
                                        
                                    #distToCPrime = np.linalg.norm(data_x[i]-self.centroids[c_prime])
                                    #distToCPrime = kmeans_common_func.euclidean_distance(data_x[i], self.centroids[c_prime])
                                    distToCPrime = kmeans_common_func.euclidean_distance(data_x[i], self.centroids[c_prime], x_sq_sum[i], c_sq_sum[c_prime], mode=self.mode)
                                        
                                    if distToCurrentCentroid > distToCPrime:
                                        assigned_centroid = c_prime 
                                        distToCurrentCentroid = distToCPrime
                        
                            self.classifications[assigned_centroid].append(data_x[i])
                            self.pointsClassif[assigned_centroid].append(i)
                                            
                prevCentroids = dict(self.centroids.copy())
                prevClassifications = dict(self.classifications.copy())
                prevPointsClassif = dict(self.pointsClassif.copy())
                
                for classification in self.classifications:
                    if len(self.classifications[classification]) == 0:
                        pass

                    else:
                        self.centroids[classification] = np.average(self.classifications[classification],axis=0)
                if self.mode == 'semi_optimized': c_sq_sum = kmeans_common_func.squared_sum(self.centroids)

                optimized = [True for centroid in self.centroids]
                
                centroidDistanceChange = {}

                for centroid in self.centroids:
                    original_centroid = prevCentroids[centroid]
                    current_centroid = self.centroids[centroid]
                    #centroidDistanceChange[centroid] = np.linalg.norm(original_centroid-current_centroid)
                    centroidDistanceChange[centroid] = kmeans_common_func.euclidean_distance(original_centroid, current_centroid)

                    #if abs(np.sum((current_centroid-original_centroid)/original_centroid*100.0)) > self.kmeans_threshold :
                    if abs(np.sum(np.divide((current_centroid-original_centroid), original_centroid, out=np.zeros_like(current_centroid), where=original_centroid!=0)*100.0)) > self.kmeans_threshold :
                        optimized[centroid] = False

                if False not in optimized:
                    break
                    
        
        ## Update labels (cluster) for each point
        for centroid in self.pointsClassif:
            for point in self.pointsClassif[centroid]:
                self.labels[point] = centroid
        
        labels = self.labels
        centroids = self.centroids
        print('converged in ' + str(iter_count) + " iterations")

        datasize = data_x.shape[1]
        total_lloyd_dist_count = (data_x.shape[0] * self.k) * iter_count * datasize
        print('Total Lloyd features in dist calculation: ', total_lloyd_dist_count)
                

        return labels, centroids
                        

def kmeans_lloyd_manual(x, k, max_iter=100, initialCentroids=None, kmeans_threshold=0, mode='fast'):
    results = kmeans_class(x, k, max_iter=100, initialCentroids=None, kmeans_threshold=0, mode=mode)
    labels = results.labels
    centroids = results.centroids
    return labels, centroids


