#In case of need: Freshly reload the import modules
#import importlib
#importlib.reload(kmeans_eco)
import math
import numpy as np


def euclidean_distance_raw(p1, p2):
    dist = 0
    for i in range(len(p1)):
        dist += ((p1[i] - p2[i]) ** 2)
    return math.sqrt(dist)

def euclidean_distance(p1, p2, sq_p1=0, sq_p2=0, mode='fast'):
    #mode can be 'fast', 'manual', 'semi_optimized'
    #print('Distmode: ', mode)
    dist = 0
    if mode == 'manual': #Everything needs to be done manually
        for i in range(len(p1)):
            dist += ((p1[i] - p2[i]) ** 2)
        dist = math.sqrt(dist)
    
    elif mode == 'fast': #Using optimized library
        dist += np.linalg.norm(p1 - p2)
    
    elif mode == 'semi_optimized': #Using optimized library
        dist += sq_p1 + sq_p2
        for i in range(len(p1)):
            dist -= 2*p1[i]*p2[i]
        if dist < 0: dist = 0 #Due to floating point precision error
        dist = math.sqrt(dist)

    return dist

def squared_sum(data):
    sum_squared = np.zeros(len(data))
    for i in range(len(data)):
        for j in range(len(data[i])):
            sum_squared[i] += (data[i][j] ** 2)
    return sum_squared

def has_converged(old_medoids, medoids):
    return set([tuple(x) for x in old_medoids]) == set([tuple(x) for x in medoids])
  

