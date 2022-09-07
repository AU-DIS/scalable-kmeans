# scalable-kmeans
KMeans but really fast


# Datasets
Datafiles will be ignored so please take care of them on your own, or add a download script for others to use.

Some can be found on odin at ```/storage/ARPESDATA/Data```.


# Converting .h5 to .csv
If csv is prefered a conversion script can be found in /Data.

Run it as
```
./h5tocsv.sh [/path/to/.h5]
```
Will output a csv with the same name at the same location.

Depends on pip packages h5py, numpy, pandas, itertools

> **Warning**
> Keep in mind that .csv will be a lot bigger, so use with care.
