#' RphenoGraph clustering
#'
#' R implementation of the PhenoGraph algorithm
#'
#' A simple R implementation of the
#' [PhenoGraph](http://www.cell.com/cell/abstract/S0092-8674(15)00637-6) algorithm, which is a
#' clustering method designed for high-dimensional single-cell data analysis. It works by creating a
#' graph ("network") representing phenotypic similarities between cells by calclating the Jaccard
#' coefficient between nearest-neighbor sets, and then identifying communities using the well known
#' [Louvain method](https://sites.google.com/site/findcommunities/) in this graph.
#'
#' @param data matrix; input data matrix.
#' @param k integer; number of nearest neighbours (default:30).
#' @param directed logical; whether to use a symmetric (FALSE, default) or asymmetric ("directed",
#'   TRUE) graph.
#' @param edge_comb character map multiple edges to single edges; it uses igraph
#'   \code{igraph-attribute-combination}. If graph is directed, there should be no duplicated edge.
#'   If graph is undirected, there are numerous edges that exist in both directions. Those edges
#'   will be combined using the \code{simplify} function, that defines which function to apply to
#'   the weights. Available values are "sum", "prod", "mean", "min", "max"... "first" is used by
#'   default as duplicated links are getting the same weight (same intersection) and is equivalent
#'   to "mean". "sum" and "prod" are reinforcing the weight of links that exist in both directions.
#' @param inc_self logical; whether to include the 1st nearest neighbour, i.e. the data point
#'   itself, when computing the shared neighbours (intersection) with its nearest neighbours. By
#'   default, it is FALSE, as the original Rphenograph implementation. No self link is created,
#'   unless the data point is linked to none of its neighbours because they share no common
#'   neighbours.
#' @param clust_fun function; community detection algorithm. Defaults to cluster_louvain. Louvain is
#'   not implemented for directed graphs. Other options: cluster_walktrap, cluster_spinglass,
#'   cluster_leading_eigen, cluster_edge_betweenness, cluster_fast_greedy, cluster_label_prop.
#' @param verbose logical; verbosity (default=FALSE)
#' @param knn_fun function or character; the function used to search the nearest neighbors. NULL
#'   points to RANN::nn2 as original code, which is configured for exact nearest neighbours (eps=0).
#'   "hnsw" points to RcppHNSW::hnsw_knn as proposed by E. Becht. Any function that returns a matrix
#'   of index to nearest neighbors could be used. NB: the first column contains the point itself and
#'   is removed before building the graph unless inc_self is TRUE.
#' @param knn_report vector of integers; this vector defines the indices within the k nearest
#'   neighbours to report for each data points. This should be a vector of integer from 1 to k. A
#'   single value of 0 indicates that all neighbours (indices from 1 to k) will be reported, as by
#'   default.
#'
#' @return a list contains an igraph graph object for \code{graph_from_data_frame} and a communities
#'   object, for which the most interesting operations are the following:
#'
#'   \item{membership}{returns a numeric vector, one number for each vertex in the graph that was
#'   the input of the community detection.}
#'
#'   \item{length}{returns an integer scalar.}
#'
#'   \item{sizes}{returns a numeric vector.}
#'
#'   \item{modularity}{returns a numeric scalar.}
#'
#'   There are also other operations: \item{algorithm}{returns a character scalar.},
#'   \item{crossing}{returns a logical vector.} \item{is_hierarchical}{returns a logical scalar.},
#'   \item{merges}{returns a two-column numeric matrix.}, \item{cut_at}{returns a numeric vector,
#'   the membership vector of the vertices.}, \item{as.dendrogram}{returns a dendrogram object.},
#'   \item{show_trace}{returns a character vector.}, \item{code_len}{returns a numeric scalar for
#'   communities found with the InfoMAP method and NULL for other methods.}, \item{print}{returns
#'   the communities object itself, invisibly.}, \item{plot}{for communities objects returns NULL,
#'   invisibly.}
#'
#' @details
#'
#' The Python Phenograph algorithm defines a "prune" parameter. prune is a logical that determines
#' whether the graph is simplified (i.e. removing reciprocal edges) by addition (FALSE, default) or
#' multiplication (TRUE). Phenograph is using matrix operations, adjacency matrix of the graph added
#' or multiplied to its transpose. Originally, Rphenograph does not simplify the graph, keeping
#' reciprocal edges. Now Rphenograph offers igraph simplifications, which might lead to different
#' results from Phenograph. In Phenograph, it's not clear what's happening when a link exists in one
#' direction: does product lead to 0 when there is only one link?
#'
#' There is no rational whether to include or not to the data point itself in the list of nearest
#' neighbours. Phenograph does not, but the scran package does. The choice is driven here by
#' \code{inc_self}. Whatever the choice, the graph could be assymmetric when a point A lists B as
#' its k NN, but B does not, because B has k neighbours nearer than A. So A has a link to B, but B
#' has not link to A. When A and B are in the list of each other, there are two reverse links; A and
#' B have the same intersection and thus the same Jaccard coefficient.
#'
#' @references Jacob H. Levine and et.al. Data-Driven Phenotypic Dissection of AML Reveals
#'   Progenitor-like Cells that Correlate with Prognosis. Cell, 2015.
#'
#' @examples
#' library(igraph)
#' iris_unique <- unique(iris) # Remove duplicates
#' data <- as.matrix(iris_unique[,1:4])
#' Rphenograph_out <- Rphenograph(data, k = 45)
#' modularity(Rphenograph_out[[2]])
#' membership(Rphenograph_out[[2]])
#' iris_unique$phenograph_cluster <- factor(membership(Rphenograph_out[[2]]))
#' if(require(ggplot2)) {
#'     ggplot(iris_unique, aes(x=Sepal.Length, y=Sepal.Width, col=Species, shape=phenograph_cluster)) +
#'         geom_point(size = 3) + theme_bw()
#' } else {
#'     with(iris_unique, plot(x=Sepal.Length, y=Sepal.Width,
#'     col=Species, pch=(15:17)[phenograph_cluster],
#'     cex=c(2.5, 2, 1.5)[phenograph_cluster]))
#' }
#'
#' @importFrom igraph graph.data.frame cluster_louvain modularity membership simplify
#' @useDynLib Rphenograph
#'
#' @export
Rphenograph <- function(data, k=30, 
                        directed=FALSE,
                        edge_comb="none",
                        inc_self=FALSE,
                        clust_fun=NULL, 
                        verbose=FALSE, 
                        knn_fun=NULL,
                        knn_report=0)
{
    if(is.data.frame(data))
        data <- as.matrix(data)
    
    if(!is.matrix(data))
        stop("Wrong input data, should be a data frame or matrix!")
    
    if(k<1){
        stop("k must be a positive integer!")
    }else if (k > nrow(data)-2){
        stop("k must be smaller than the total number of points!")
    }

    if(verbose){
        message("\nRun Rphenograph starts:\n", 
                "  -Input data of ", nrow(data)," rows and ", ncol(data), " columns","\n",
                "  -k is set to ", k)
        cat("  Finding nearest neighbors...")
    }
    if(is.null(knn_fun)) knn_fun <- find_neighbors
    if(is.character(knn_fun) && knn_fun == "hnsw")
      if (requireNamespace("RcppHNSW",  quietly = TRUE)) {
        knn_fun <- function(data, k) RcppHNSW::hnsw_knn(data, k=k)$idx
      } else
        stop("RcppHNSW must be installed for using hnsw nearest neighbors search.")
    if(!is.function(knn_fun))
        stop("knn_fun must be a function or a recognized keyword.")
    ## t1 <- system.time(neighborMatrix <- find_neighbors(data, k=k+1)[,-1])
    t1 <- system.time({
      if (inc_self) {
        neighborMatrix <- knn_fun(data, k)
      } else {
        neighborMatrix <- knn_fun(data, k+1)[,-1]
      }
    })

    if(verbose){
      cat("DONE ~",t1[3],"s\n")
      cat("  Compute jaccard coefficient between nearest-neighbor sets...")
    }
    t2 <- system.time({
      links <- jaccard_coeff_4(neighborMatrix, knn_report)
      links <- links[links[,1]>0, ]  # remove unused allocation (unreported zero Jaccard coeff)
      ## Fix if data point goes missing (due to all of its associated jaccard
      ## coefficients being 0 and if it ever appears as another points' nearest
      ## neighbor, the corresponding jaccard coefficient also being 0.
      ## Now done in the C code
    })

    if(verbose){
        cat("DONE ~",t2[3],"s\n", " Build a graph from the weighted links...")
    }

    t3 <- system.time({
      relations <- as.data.frame(links)
      colnames(relations) <- c("from","to","weight")
      g <- graph.data.frame(relations, directed=directed)
      if (directed && edge_comb != "none") {
        cat("  There are typically no duplicated edge when the graph is directed.")
      }
      if (edge_comb == "none") {
        if(verbose)
          cat("\n", " Links are not combined...")
      } else {
        if(verbose)
          cat("\n", " Combining links using", edge_comb, "method...")
        g <- simplify(g, edge.attr.comb = edge_comb, remove.loops = FALSE)
        # Loops are kept to preserve unlinked data points
      }
    })

    if(verbose) {
        cat("DONE ~",t3[3],"s\n")
        cat("  Run ", "clustering on the graph ...")  # TODO: add function name
    }

    # Other community detection algorithms: 
    #    cluster_walktrap, cluster_spinglass, 
    #    cluster_leading_eigen, cluster_edge_betweenness, 
    #    cluster_fast_greedy, cluster_label_prop  
    if(is.null(clust_fun)) clust_fun <- cluster_louvain
    
    t4 <- system.time(community <- clust_fun(g))

    if(verbose){
        cat("DONE ~",t4[3],"s\n")
        
        message("Run Rphenograph DONE, totally takes ", 
                format(sum(c(t1[3],t2[3],t3[3],t4[3])), digits = 5), " s.")
        cat("  Return a community class\n")
        cat("  -Modularity value:", modularity(community), "\n")
        cat("  -Number of clusters:", length(unique(membership(community))), "\n")
    }
    
    return(list(g, community))
}

#' K Nearest Neighbour Search
#'
#' Uses a kd-tree to find the p number of near neighbours for each point in an input/output dataset.
#' 
#' Use the nn2 function from the RANN package, utilizes the Approximate Near Neighbor (ANN) C++ library, 
#' which can give the exact near neighbours or (as the name suggests) approximate near neighbours 
#' to within a specified error bound. For more information on the ANN library please 
#' visit http://www.cs.umd.edu/~mount/ANN/.
#' 
#' @describeIn Rphenograph Uses a kd-tree to find the p number of near neighbours for each point in an input/output dataset.
#' 
#' @return a n-by-k matrix of neighbor indices
#' 
#' @examples
#' iris_unique <- unique(iris) # Remove duplicates
#' data <- as.matrix(iris_unique[,1:4])
#' neighbors <- find_neighbors(data, k=10)
#' 
#' @importFrom RANN nn2
#' @export
find_neighbors <- function(data, k){
    nearest <- nn2(data, data, k, searchtype = "standard")
    return(nearest[[1]])
}
