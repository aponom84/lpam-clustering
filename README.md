# Link Partitioning Around Medoids (lpam) – Tool for Overlaping Communities Detection in Networks

## What it is? ##

This and java implementation of the method for finding overllaping communities  in network.
Below you can find an example of method output for "School friendship" dataset.
### School friendship ###
![school friendship](https://github.com/aponom84/lpam-clustering/blob/master/final_pictures/school-2_ACM_pmp_7_out.png)

### American Football Club ###
![fooball club](https://github.com/aponom84/lpam-clustering/blob/master/final_pictures/footballTSEinput_CM_kmd_12_out.png)

./Scripts – directory for all Jupyter notebook that were used for computation experimenst
./final_pictures – contains all pictures and gephi files.
./lpam – source code of lpam method
./datasets – contains realworld datasets and sysntetic datasets 
./related_methids – 
./literature – papers about overllaping communities detection methods in networks

## Compilation ##
In order to repoduce all computation experiments you should compile all method in the related_methods directory 
and tool for measuring ONMI in the directory Overlapping-NMI. 
Just go to corresponding folder and type "cmake".

```
cd Overlapping-NMI
make
```
...TODO

## Contributors ##

* Leonid Peculias: idea of overlaping communities detection based on the link partitoing with the help of non-overlaping communities detection methods based on the partition around medoids
* Marat Shamshetdinov (m.shamshetdinov@gmail.com): Implementation of exact model for lp_solver and cplex solvers.
* Nikit Putehin: implementation of amplified commute distance and heuristics method for findg k-medoids (Clarance and k-meanns)
* Alexander Ponomarenko (aponom84@gmail.com): Basic implementation of link partitionig algorithm on java, calculating commute distance, jupyter notebooks and collecting everything together
