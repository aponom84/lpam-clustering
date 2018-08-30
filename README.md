# Link Partitioning Around Medoids (lpam) – Tool for Overlaping Communities Detection in Networks

## What it is? ##

This and java implementation of the method for finding overllaping communities  in networks

## Compilation ##
In order to repoduce all computation experiments you should compile all method in the related_methods directory 
and tool for measuring ONMI in the directory Overlapping-NMI. Just go to corresponding folder and type "cmake".

## Contributors ##

* Leonid Peculias: idea of overlaping communities detection based on the link partitoing with the help of non-overlaping communities detection methods based on the partition around medoids
* Marat Shamshetdinov (m.shamshetdinov@gmail.com): Implementation of exact model for lp_solver and cplex solvers.
* Nikit Putehin: implementation of amplified commute distance and heuristics method for findg k-medoids (Clarance and k-meanns)
* Alexander Ponomarenko (aponom84@gmail.com): Basic implementation of link partitionig algorithm on java, calculating commute distance, jupyter notebooks and collecting everything together
