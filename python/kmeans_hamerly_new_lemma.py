import numpy as np

class kmeans_class:
    def __init__(self, data_points, k, max_iter=100, initialCentroids=None, kmeans_threshold=0, print='no'):
        # k: no. of clusters
        # data_points: data to cluster
        # max_iter: iteration before termination
        # kmeans_threshold: termination condition to detect convergence if the centroids do not move more than the 'kmeans_threshold' value
        
        self.k = k
        self.data_points = data_points
        self.max_iter = max_iter
        self.kmeans_threshold = kmeans_threshold
        self.method = 'Hamerly_New_Lemma'

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
        iter_count = 0   

        ## Hamerly's Kmeans
        if self.method == 'Hamerly_New_Lemma':
            
            self.classifications = {} ## Points coordinates by centroid
            self.pointsClassif = {} ## Points indices by centroid
            
            lowerBounds = np.zeros((data_x.shape[0],self.k)) ## Lower bounds matrix. Dimensions : (nb_points, nb_centroids)
            upperBounds = np.zeros((data_x.shape[0])) ## Upper bounds vector : Dimension : (nb_points)
            lowerBounds_Hamerly = np.zeros((data_x.shape[0])) ## Hamerly's lower bounds vector : Dimension : (nb_points)


            for i in range(self.k):
                self.classifications[i] = []
                self.pointsClassif[i] = []
                
            i=0
            for point in data_x:
                distances = [np.linalg.norm(point-self.centroids[centroid]) for centroid in self.centroids]
                classification = distances.index(min(distances))
                    
                ## Centroid assigned to each point
                self.classifications[classification].append(point)
                self.pointsClassif[classification].append(i)
                    
                ## Lower bound distance between the point and each center
                ## Initialized as distance between the point and each initial centroid
                lowerBounds[i] = distances
                    
                # For Hamerly's second lower bound
                ## Second Lower bound distance between the point and the second closest center
                ## Initialized as distance between the point and second closest initial centroid
                lowerBounds_Hamerly[i] = sorted(set(distances))[1]
                #second_nearest_cluster = distances.index(lowerBounds_Hamerly[i])

                ## Upper bound distance between the point and assigned centroid
                upperBounds[i] = min(distances)
                    
                i+=1
                
            prevCentroids = dict(self.centroids.copy())
            prevClassifications = dict(self.classifications.copy())
            prevPointsClassif = dict(self.pointsClassif.copy())
            

#            for classification in self.classifications:
#                if len(self.classifications[classification]) == 0:
#                    pass

#                else:
#                    self.centroids[classification] = np.average(self.classifications[classification],axis=0)

#            optimized = [True for centroid in self.centroids]
            
            centroidDistanceChange = {}

            for centroid in self.centroids:
                original_centroid = prevCentroids[centroid]
                current_centroid = self.centroids[centroid]
                centroidDistanceChange[centroid] = np.linalg.norm(original_centroid-current_centroid)
                
#                if abs(np.sum((current_centroid-original_centroid)/original_centroid*100.0)) > self.kmeans_threshold :
#                    optimized[centroid] = False

#            if False not in optimized:
                
#                for centroid in self.pointsClassif:
#                    for point in self.pointsClassif[centroid]:
#                        self.labels[point] = centroid
#                return
            
            ## Update lower and upper bound distances
            #for centroid in self.pointsClassif:
            #    for i in list(range(data_x.shape[0])):
            #        lowerBounds[i][centroid] -= centroidDistanceChange[centroid]
                    
            #    for i in self.pointsClassif[centroid]:
            #        upperBounds[i] += centroidDistanceChange[centroid]
            #        #lowerBounds_Hamerly[i] += centroidDistanceChange[second_nearest_cluster]
            #        lowerBounds_Hamerly[i] += max(centroidDistanceChange)

                    
            ## Repeat until convergence
            new_lemma_count = 0  #Number of discarded distance calculations
            hamerly_count = 0  #Number of discarded distance calculations
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
                
                for centroid in prevPointsClassif:
                    for i in prevPointsClassif[centroid]:
                    #for i in range(data_x.shape[0]):
                        
                        r = True
                        distToCurrentCentroid = upperBounds[i]
                        
                        #Elkan Lemma1 
                        ## Check if upper bound lower than 1/2 of distance with closest centroid
                        #if upperBounds[i] <= 0.5*closestCentroidDistances[centroid]:
                        #    ## If condition is met : said point keeps its centroid with no further computation needed
                        #    self.classifications[centroid].append(data_x[i])
                        #    self.pointsClassif[centroid].append(i)
                        #    new_lemma_count += 1

                        hamerly_bound = max(0.5*closestCentroidDistances[centroid], lowerBounds_Hamerly[i])
                        if upperBounds[i] <= hamerly_bound: #For Hamerly first check
                        #elif upperBounds[i] <= lowerBounds_Hamerly[i]: #For Hamerly
                            ## If condition is met : said point keeps its centroid with no further computation needed
                            self.classifications[centroid].append(data_x[i])
                            self.pointsClassif[centroid].append(i)
                            hamerly_count += 1
                        
                        else:
                            if r:
                                distToCurrentCentroid = np.linalg.norm(data_x[i] - self.centroids[centroid])
                                #lowerBounds[i][centroid] = distToCurrentCentroid
                                r = False
                            upperBounds[i] = distToCurrentCentroid #Tighten the upper bound
                            if upperBounds[i] <= hamerly_bound: #For Hamerly second check
                            #elif upperBounds[i] <= lowerBounds_Hamerly[i]: #For Hamerly
                                ## If condition is met : said point keeps its centroid with no further computation needed
                                self.classifications[centroid].append(data_x[i])
                                self.pointsClassif[centroid].append(i)
                                hamerly_count += 1
                                
                            else:
                                assigned_centroid = centroid
                                for c_prime in (listCentroids[:centroid]+listCentroids[centroid+1:]):
                                    #Proposed New Lemma
                                    if lowerBounds[i][c_prime] > (upperBounds[i] + centroidDistanceChange[assigned_centroid] + centroidDistanceChange[c_prime]):
                                        ## If condition is met : said point keeps its centroid with no further computation needed
                                        new_lemma_count += 1

                                    #Elkan lemma2
                                    ## Check if lower bound between point and c_prime < upper bound between point and its current centroid
                                    ## AND if (0.5*distance between current centroid and c_prime) < upper bound between point and its current centroid
                                    #if ((distToCurrentCentroid > lowerBounds[i][c_prime]) and (distToCurrentCentroid > 0.5*centroidDistances[centroid][c_prime])): 
                                        #if r:
                                        #    distToCurrentCentroid = np.linalg.norm(data_x[i] - self.centroids[centroid])
                                        #    lowerBounds[i][centroid] = distToCurrentCentroid
                                        #    r = False
                                    else:    
                                        distToCPrime = np.linalg.norm(data_x[i]-self.centroids[c_prime])
                                        #lowerBounds[i][c_prime] = distToCPrime
                                        
                                        if distToCurrentCentroid > distToCPrime:
                                            #second_nearest_cluster = assigned_centroid
                                            #lowerBounds_Hamerly[i] = distToCurrentCentroid
                                            assigned_centroid = c_prime 
                                            upperBounds[i] = distToCPrime
                                            lowerBounds_Hamerly[i] = distToCurrentCentroid
                                            distToCurrentCentroid = distToCPrime
                                        else:
                                            if lowerBounds_Hamerly[i] > distToCPrime:
                                                lowerBounds_Hamerly[i] = distToCPrime
                                            
                                        
                                    #else:
                                    #    new_lemma_count += 1 #Meaning one distance caluculation has been discarded
                        
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
                    
                ## Update of lower and upper bound distances
                maxCentroidDistanceChange = max(centroidDistanceChange.values())
                for centroid in self.pointsClassif:
                    for i in list(range(data_x.shape[0])):
                        lowerBounds[i][centroid] -= centroidDistanceChange[centroid]
                    
                    for i in self.pointsClassif[centroid]:
                        upperBounds[i] += centroidDistanceChange[centroid]
                        #lowerBounds_Hamerly[i] += centroidDistanceChange[second_nearest_cluster]
                        furthestMovingCentroid = max(centroidDistanceChange, key=centroidDistanceChange.get)
                        if furthestMovingCentroid == centroid:
                            lowerBounds_Hamerly[i] -= sorted(set(centroidDistanceChange.values()), reverse=True)[1]
                        else:
                            lowerBounds_Hamerly[i] -= maxCentroidDistanceChange
        
        ## Update labels (cluster) for each point
        for centroid in self.pointsClassif:
            for point in self.pointsClassif[centroid]:
                self.labels[point] = centroid
        
        labels = self.labels
        centroids = self.centroids
        print('converged in ' + str(iter_count) + " iterations")
        print('Hamerly dist calculation discarded: ', hamerly_count)
        print('New_Lemma dist calculation discarded: ', new_lemma_count)

        datasize = data_x.shape[1]
        total_lloyd_dist_count = (data_x.shape[0] * self.k) * iter_count * datasize
        total_hamerly_dist_count = total_lloyd_dist_count - (hamerly_count * datasize)
        total_new_lemma_dist_count = total_lloyd_dist_count - ((hamerly_count + new_lemma_count) * datasize)
        print('Total Lloyd dist features calculation: ', total_lloyd_dist_count)
        print('Total Hamerly dist features in calculation: ', total_hamerly_dist_count)
        print('Total New_Lemma dist features in calculation: ', total_new_lemma_dist_count)
                

        return labels, centroids
                        

def kmeans_hamerly_new_lemma_manual(x, k, max_iter=100, initialCentroids=None, kmeans_threshold=0):
    results = kmeans_class(x, k, max_iter=100, initialCentroids=None, kmeans_threshold=0)
    labels = results.labels
    centroids = results.centroids
    return labels, centroids


