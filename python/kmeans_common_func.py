#In case of need: Freshly reload the import modules
#import importlib
#importlib.reload(kmeans_eco)


def euclidean_distance(p1, p2):
    dist = 0
    for i in range(len(p1)):
        dist += ((p1[i] - p2[i]) ** 2)
    return dist

def has_converged(old_medoids, medoids):
    return set([tuple(x) for x in old_medoids]) == set([tuple(x) for x in medoids])
  

