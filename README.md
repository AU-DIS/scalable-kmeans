# scalable-kmeans
KMeans. Includes different python and cpp implementations.

A cleaner fast python wrapper of MARIGOLD can be found here: 
https://github.com/AU-DIS/pyMARIGOLD

# Description of python code
Python code includes six algorithms (mentioned in [1]): Lloyed, Elkan, Hamerly, Elkan with Hamerly, Stepwise and MARIGOLD. The codes can be executed from the notebook. 


# Data
Raw .h5 files for misfit data sets have been made available through Git LFS: https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-git-large-file-storage  
Installation guide: https://docs.github.com/en/repositories/working-with-files/managing-large-files/installing-git-large-file-storage

They can be found in Data/ folder and converted to another format using pandas. *NOTE:* The provided files have not been DCT transformed. You must do so first to obtain the best performance from MARIGOLD.  


[DeepGlobe](https://openaccess.thecvf.com/content_cvpr_2018_workshops/papers/w4/Demir_DeepGlobe_2018_A_CVPR_2018_paper.pdf) Land Cover can be found online. Here is a current link: https://www.kaggle.com/datasets/balraj98/deepglobe-land-cover-classification-dataset/data
