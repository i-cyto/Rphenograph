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


// Slightly improve ressources usage
// Author: S. Granjeaud, Date: 17/02/2020

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


// Improve intersection computation by pre-sorting the indices, what requires more memory
// Author: S. Granjeaud, Date: 17/02/2020

// [[Rcpp::export]]
NumericMatrix jaccard_coeff_3(NumericMatrix idx) {
    int nrow = idx.nrow(), ncol = idx.ncol();
    NumericMatrix weights(nrow*ncol, 3);
    int r = 0;
    // allocate memory
    int** i_idx = new int*[nrow];
    if (nrow) {
        i_idx[0] = new int[nrow * ncol];  // mem alloc
        for (int i = 1; i < nrow; ++i)
            i_idx[i] = i_idx[0] + i * ncol;  // fill pointer on 1st element of each row
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


// Improve intersection computation by pre-sorting the indices, what requires more memory
// Report a lower number of nearest neighbors, what reduces the graph size and lowers errors on iris
// Author: S. Granjeaud, Date: 24/05/2020

// [[Rcpp::export]]
NumericMatrix jaccard_coeff_4(NumericMatrix idx, IntegerVector ks) {
    int nrow = idx.nrow(), ncol = idx.ncol();
    int nks = ks.length();
    // ks can't be empty
    if (nks == 0)
        stop("ks must be a vector integer with one element at least.");
    int* kx = new int[ncol];  // ncol is the maximal length
    int nkx;
    if (nks == 1 && ks[0] == 0) {
        // if ks is c(0), back to full list of column indices, aka auto mode
        for (int i = 0; i < ncol ; i++) {
            kx[i] = i;
        }
        nkx = ncol;
    } else {
        // subtract 1 from ks indices for C computation, and check boundaries
        bool err_lower = 0, err_greater = 0;
        for (int i = 0; i < nks ; i++) {
            kx[i] = ks[i] - 1;
            if (kx[i] < 0 && !err_lower) err_lower = 1;
            if (kx[i] >= ncol && !err_greater) err_greater = 1;
        }
        nkx = nks;
        if (err_lower || err_greater)
            stop("ks must a vector of integers from 1 to ncol(idx).");
    }
    // Result and its index
    NumericMatrix weights(nrow*nkx, 3);
    int r = 0;
    // allocate memory
    int** i_idx = new int*[nrow];
    i_idx[0] = new int[nrow * ncol];  // mem alloc
    for (int i = 1; i < nrow; ++i)
        i_idx[i] = i_idx[0] + i * ncol;  // fill pointer on 1st element of each row
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
        for (int j = 0; j < nkx; j++) {
            //k = i_idx[i][ks[j]]-1;
            k = idx(i,kx[j])-1;
            if (i == k) continue;  // avoid loop
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
            if (u < 1) continue;  // skip this link if no intersection
            weights(r, 0) = i+1;
            weights(r, 1) = k+1;
            weights(r, 2) = u/(4.0*ncol - 2*u);  // symmetrize the graph
            // symmetrization introduces an extra division by 2, but not clear why this is needed
            r++;
        }
    }
    // free matrix
    delete [] i_idx[0];
    delete [] i_idx;
    delete [] kx;
    // report all nodes (E.K. Becht's issue)
    unsigned long int *reported = new unsigned long int[nrow+1];  // using R indices, 0 unused
    for (int i = 1; i <= nrow; i++)
        reported[i] = 0;
    // loop to check all nodes
    for (int i = 0; i < r; i++) {
        reported[(unsigned long int)weights(i, 0)]++;
        reported[(unsigned long int)weights(i, 1)]++;
    }
    // add all unreported nodes
    for (int i = 1; i <= nrow; i++) {
        if (reported[i] > 0) continue;
        weights(r, 0) = i;
        weights(r, 1) = i;
        weights(r, 2) = 1/2.0;
        // extra division by 2, but not clear why this is needed
        r++;
    }
    delete [] reported;
    // End issue
    // done
    return weights;
}

// DEBUGGING R code, see vignette
