//------------------------------------------------------------------------------
// LAGraph_dfs_traverse_bin_tree:  binary tree traverse in order which is providen by dfs
//------------------------------------------------------------------------------

/*
    LAGraph:  graph algorithms based on GraphBLAS

    Copyright 2019 LAGraph Contributors.

    (see Contributors.txt for a full list of Contributors; see
    ContributionInstructions.txt for information on how you can Contribute to
    this project).

    All Rights Reserved.

    NO WARRANTY. THIS MATERIAL IS FURNISHED ON AN "AS-IS" BASIS. THE LAGRAPH
    CONTRIBUTORS MAKE NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
    AS TO ANY MATTER INCLUDING, BUT NOT LIMITED TO, WARRANTY OF FITNESS FOR
    PURPOSE OR MERCHANTABILITY, EXCLUSIVITY, OR RESULTS OBTAINED FROM USE OF
    THE MATERIAL. THE CONTRIBUTORS DO NOT MAKE ANY WARRANTY OF ANY KIND WITH
    RESPECT TO FREEDOM FROM PATENT, TRADEMARK, OR COPYRIGHT INFRINGEMENT.

    Released under a BSD license, please see the LICENSE file distributed with
    this Software or contact permission@sei.cmu.edu for full terms.

    Created, in part, with funding and support from the United States
    Government.  (see Acknowledgments.txt file).

    This program includes and/or can make use of certain third party source
    code, object code, documentation and other files ("Third Party Software").
    See LICENSE file for more details.

*/

//------------------------------------------------------------------------------

// LAGraph_dfs_traverse_bin_tree:  TODO: WRITE IT
//  based on the breadth-first search in the GraphBLAS C
// API Specification by Scott McMillan, CMU.  Modified by Tim Davis, Texas A&M.

// Perform a single-source BFS, starting at a source node.  Returns a dense
// vector v such that v(i) > 0 if node is reachable from the source node.
// v(source)=1 and v(i)=k if the path with the fewest edges from the source
// node to i has k-1 edges.  If i is not reachable from the source node, then
// v(i) is zero.

// This method is a simple one for illustration only, and works well in
// practice, except in the following cases:

// (1) It takes Omega(n) time.  If nvals(v) << n is expected, use
// LAGraph_bfs_pushpull instead, which is much faster if v is expected to be
// very sparse.

// (2) It assumes that vxm(q,A) is fast, and implemented in a 'push' fashion,
// using saxpy operations instead of dot products.  This requires that the rows
// A(i,:) are efficient to access, which is the case if A is in CSR format
// internally in GraphBLAS (or perhaps in both CSR and CSC formats).  This
// method will be *exceedlingly* slow if A is a MATLAB matrix in CSC format.

// See LAGraph_bfs_pushpull, which handles the above two cases.

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_ALL    \
{                           \
    GrB_free (&v) ;         \
}


GrB_Info build_permutation_matrix(GrB_Matrix *perm_matrix, int *children_amount, const GrB_Vector *children_pos, GrB_Index tay, GrB_Index beta) {
    // as we don't have ability to compute children_pos * children_pos.transpose()
    // we will use reduce with plus to calculate amount of not zero values
    GrB_Vector v = NULL;

    LAGr_reduce(&children_amount, GrB_NULL, GxB_PLUS_INT64_MONOID, *children_pos, GrB_NULL);
}

GrB_Info LAGraph_dfs_traverse_bin_tree     // traversal of the binary tree in pre-order DFS
        (
                GrB_Vector *v_output,   // v(i) is the number of node in traverse order
                GrB_Matrix A,           // input graph, treated as if boolean in semiring
                GrB_Index source        // starting node of the DFS
        ) {

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    GrB_Vector v = NULL;

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    // create an empty vector v, and make it dense
    GrB_Index n;
    GrB_BinaryOp empty_oper = NULL;
    LAGr_Matrix_nrows(&n, A);
    LAGr_Vector_new (&v, GrB_INT64, n);
    for (GrB_Index i = 0; i < n; i++) {
        LAGr_Vector_setElement(v, i, i);
    }

    int children_amount;
    GrB_Index tay = 0, beta = 0;
    GrB_Matrix permutation_matrix = NULL;
    // create matrix with shape like A.shape

    build_permutation_matrix(&permutation_matrix, &children_amount, &v, tay, beta);
    printf("CGHHHH: %d\n", children_amount);


    // create a boolean vector q, and set q(source) to true
//    LAGr_Vector_new (&q, GrB_BOOL, n) ;
//    LAGr_Vector_setElement (q, true, source) ;

    //--------------------------------------------------------------------------
    // BFS traversal and label the nodes
    //--------------------------------------------------------------------------

//    for (int64_t level = 1 ; level <= n ; level++)
//    {
//        // v<q> = level
//        LAGr_assign (v, q, NULL, level, GrB_ALL, n, desc_s) ;
//
//        // break if q is empty
//        LAGr_Vector_nvals (&nvals, q) ;
//        if (nvals == 0) break ;
//
//        // q'<!v> = q'*A
//        LAGr_vxm (q, v, NULL, semiring, q, A, desc_rc) ;
//    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*v_output) = v;       // return result
    v = NULL;              // set to NULL so LAGRAPH_FREE_ALL doesn't free it
    LAGRAPH_FREE_ALL;      // free all workspace (except for result v)
    return (GrB_SUCCESS);
}

