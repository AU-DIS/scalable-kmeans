#ifndef KMEANSSTRATEGY_H
#define KMEANSSTRATEGY_H

class KmeansStrategy {
    public:
      // pure virtual function
        virtual ~KmeansStrategy() = default;
        virtual void run();
        virtual int get_max_iter();
        virtual int get_data_dims();
        virtual int get_data_points();
        virtual int get_num_clusters();
      
    private:

        int max_iter;       // Maximum number of iterations
        int data_dims;      // Dimensions of each datapoint 
        int data_points;    // Number of datapoints   
        int num_clusters;   // Number of clusters
                          
};

#endif