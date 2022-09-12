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

# Notes about C++
## Dependencies
C packages can be confusing to install if you do it directly form the source, but most can be found with ```apt``` instead.  

Install Blas with ```sudo apt-get install libblas-dev liblapack-dev``` and include as ```#include <cblas.h>```. Then you are done.

## Print statements
We all know prints are constly, but there are a few tricks we can utilize in C++.

**First**, if a line does not need to be printet in real-time, use ```'\n'``` instead of ```std::endl```.
```std::cout << std::endl;``` is the same as ```std::cout << '\n' << std::flush;```.

```c++
for (int n : {0, 1, 2, 3, 4, 5}) 
        std::cout << n << std::endl;
```
vs 
```c++
for (int n : {0, 1, 2, 3, 4, 5}) 
        std::cout << n << '\n';
std::cout << std::flush;
```

**Secondly**, We can utilize macros to only include conditional prints when debugging instead of having a flag check on runtime.
```c++
#ifdef DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 0
#endif
#define DEBUGPRINT(fmt, ...) \
            do { if (DEBUG_TEST) fprintf(stderr, fmt, __VA_ARGS__); } while (0)
```

Replacing ```std::cout << "Hello" << std::endl;``` with ```DEBUGPRINT("Hello")``` will only add the prints if compiled with flag ```DEBUG=1```.

Example with clang: ```clang++ kmeans.cpp -D DEBUG=1```
>**Note** the ```do{} while(0)``` and the condition check on ```DEBUG_TEST``` might seem redundant, but they are failsafes. 

