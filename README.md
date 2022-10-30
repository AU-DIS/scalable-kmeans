# scalable-kmeans
KMeans. Includes different python and cpp implementations 


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

Output example:
```
,x,y,z,k,v
0,0,0,0,0,14.288461685180664
1,0,0,0,1,22.503585815429688
2,0,0,0,2,21.361770629882812
3,0,0,0,3,19.15595054626465
4,0,0,0,4,17.359394073486328
5,0,0,0,5,18.26104736328125
6,0,0,0,6,15.660025596618652
7,0,0,0,7,13.00551700592041
8,0,0,0,8,16.275680541992188
9,0,0,0,9,16.462200164794922
...
11010047,11,13,255,255,7.41770601272583
```
Where x,y,z,k is the coordinate and v is the value.

> **Warning**
> Keep in mind that .csv will be a lot bigger, so use with care.

