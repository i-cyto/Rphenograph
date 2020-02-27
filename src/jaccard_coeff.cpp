#include <Rcpp.h>
using namespace Rcpp;

// Compute jaccard coefficient between nearest-neighbor sets
//
// Weights of both i->j and j->i are recorded if they have intersection. In this case
// w(i->j) should be equal to w(j->i). In some case i->j has weights while j<-i has no
// intersections, only w(i->j) is recorded. This is determinded in code `if(u>0)`. 
// In this way, the undirected graph is symmetrized by halfing the weight 
// in code `weights(r, 2) = u/(2.0*ncol - u)/2`.
//
// Author: Chen Hao, Date: 25/09/2015


// [[Rcpp::export]]
NumericMatrix jaccard_coeff(NumericMatrix idx) {
    int nrow = idx.nrow(), ncol = idx.ncol();
    NumericMatrix weights(nrow*ncol, 3);
    int r = 0;
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            int k = idx(i,j)-1;
            NumericVector nodei = idx(i,_);
            NumericVector nodej = idx(k,_);
            int u = intersect(nodei, nodej).size();  // count intersection number
            if(u>0){ 
                weights(r, 0) = i+1;
                weights(r, 1) = k+1;
                weights(r, 2) = u/(2.0*ncol - u)/2;  // symmetrize the graph
                r++;
            }
        }
    }
    
    return weights;
}

// Author: Samuel Granjeaud, Date: 17/02/2020
// Slightly improve ressources usage

// [[Rcpp::export]]
NumericMatrix jaccard_coeff_1(NumericMatrix idx) {
    int nrow = idx.nrow(), ncol = idx.ncol();
    NumericMatrix weights(nrow*ncol, 3);
    int r = 0;
    NumericVector nodei (ncol);
    NumericVector nodej (ncol);
    for (int i = 0; i < nrow; i++) {
        nodei = idx(i,_);
        for (int j = 0; j < ncol; j++) {
            int k = idx(i,j)-1;
            nodej = idx(k,_);
            int u = intersect(nodei, nodej).size();  // count intersection number
            if (u > 0) {
                weights(r, 0) = i+1;
                weights(r, 1) = k+1;
                weights(r, 2) = u/(4.0*ncol - 2*u);  // symmetrize the graph
                r++;
            }
        }
    }
    
    return weights;
}


// Author: Samuel Granjeaud, Date: 17/02/2020
// Improve intersection computation

// [[Rcpp::export]]
NumericMatrix jaccard_coeff_3(NumericMatrix idx) {
    int nrow = idx.nrow(), ncol = idx.ncol();
    NumericMatrix weights(nrow*ncol, 3);
    int r = 0;
    // allocate memory
    int** i_idx = new int*[nrow];
    if (nrow) {
        i_idx[0] = new int[nrow * ncol];
        for (int i = 1; i < nrow; ++i)
            i_idx[i] = i_idx[0] + i * ncol;
    }
    // copy to integer matrix and sort by row
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            i_idx[i][j] = (int)idx(i,j);
        }
        std::sort(i_idx[i], i_idx[i]+ncol);
    }
    // compute Jaccard
    int *p_i, *p_k, k;
    for (int i = 0; i < nrow; i++) {
        p_i = i_idx[i];
        for (int j = 0; j < ncol; j++) {
            k = i_idx[i][j]-1;
            p_k = i_idx[k];
            // compute intersection of sorted vectors
            int u = 0;
            int i_i = 0, i_k = 0;
            for (;;) {
                if (p_i[i_i] == p_k[i_k]) {
                    u++;
                    i_i++; i_k++;
                    if (i_i >= ncol || i_k >= ncol) break;
                }
                else if (p_i[i_i] < p_k[i_k]) {
                    i_i++;
                    if (i_i >= ncol) break;
                }
                else {
                    i_k++;
                    if (i_k >= ncol) break;
                }
            }
            // compute Jaccard
            if (u < 1) continue;
            weights(r, 0) = i+1;
            weights(r, 1) = k+1;
            weights(r, 2) = u/(4.0*ncol - 2*u);  // symmetrize the graph
            r++;
        }
    }
    // free matrix
    if (nrow) delete [] i_idx[0];
    delete [] i_idx;
    // done
    return weights;
}

// iris_data = iris[,-5]; nn = Rphenograph::find_neighbors(iris_data, 15); links = Rphenograph:::jaccard_coeff(nn); links1 = Rphenograph:::jaccard_coeff_1(nn); links3 = Rphenograph:::jaccard_coeff_3(nn); cbind(links, links1, links3)
// cbind(links[order(links[,1], links[,2]),], links1[order(links1[,1], links1[,2]),], links3[order(links3[,1], links3[,2]),])
