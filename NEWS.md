# Rphenograph 0.99.1.9003


## New features

* S. Thomas Kelly added pruning and graphs clustering methods. This ahs been turned into calling igraph functions.

* Etienne K. Becht added approximate HNSW nearest neighbors for speed.

* S. Granjeaud improves the the Jaccard_coefficient function by pre-sorting nearest indices. Now the computation takes only a few seconds for a dataset of 300 k datapoints and 30 nearest neighbors.

* Louvain is the default graph clustering method. Any functions of the (r)igraph package can be specified.


## Bug fixes and minor improvements

* EK Becht reported that some points are missing in the grap. Self loop are added in the Jaccard C code.


# Rphenograph 0.99.1

Release from H. Chen who created and published this package.
