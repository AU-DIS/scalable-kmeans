

double euclidian_distance(int x, int c, int d, int k, double data[], double centroids[]) {
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