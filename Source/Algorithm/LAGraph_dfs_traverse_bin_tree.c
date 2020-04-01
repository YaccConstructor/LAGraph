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


// it is empty because all functions can call it but only main has real ones
#define LAGRAPH_FREE_ALL


GrB_Info build_identity_matrix(GrB_Matrix *identity, GrB_Index size) {
    LAGr_Matrix_new(identity, GrB_INT32, size, size);
    LAGr_assign (*identity, NULL, NULL, 0, GrB_ALL, size, GrB_ALL, size, NULL);
    for (int i = 0; i < size; i++) {
        LAGr_Matrix_setElement(*identity, 1, i, i);
    }
}

GrB_Info shift_cols_in_matrix(GrB_Matrix matrix, GrB_Index start, GrB_Index end, GrB_Index shift_step, GrB_Index size) {
    // create tmp matrix
    GrB_Matrix tmp_matrix = NULL;
    GrB_Index tmp_matrix_ncols = end - start;
    LAGr_Matrix_new(&tmp_matrix, GrB_INT32, size, tmp_matrix_ncols);

    // create indexes of columns which will be shifted
    GrB_Index * col_indexes = malloc(sizeof(GrB_Index) * tmp_matrix_ncols);
    for (int j = 0; j < tmp_matrix_ncols; j++) {
        col_indexes[j] = start + j;
    }

    // extract all columns from tmp_beta column to ith column for shifting columns
    LAGr_extract (tmp_matrix, NULL, NULL, matrix, GrB_ALL, size, col_indexes, tmp_matrix_ncols, GrB_NULL);

    // shift columns in matrix
    for (int j = 0; j < tmp_matrix_ncols; j++) {
        col_indexes[j] += shift_step;
    }
    LAGr_assign(matrix, GrB_NULL, GrB_NULL, tmp_matrix, GrB_ALL, size, col_indexes, tmp_matrix_ncols, GrB_NULL);

    // free
    free(col_indexes);
    GrB_free(&tmp_matrix);
}


GrB_Info build_permutation_matrix(
        GrB_Matrix *perm_matrix,
        int32_t *children_amount,
        const GrB_Vector children_pos,
        GrB_Index tay,
        GrB_Index beta,
        GrB_Index size
) {
    // as we don't have ability to compute children_pos * children_pos.transpose()
    // we will use reduce with plus to calculate amount of not zero values
    LAGr_reduce(children_amount, GrB_NULL, GxB_PLUS_INT32_MONOID, children_pos, GrB_NULL);
    GrB_Matrix step_1_matrix = NULL;
    build_identity_matrix(&step_1_matrix, size);
    GrB_Index tmp_beta = beta;

    // step 1
    for (GrB_Index i = 0; i < size; i++) {
        int32_t elem;
        LAGr_Vector_extractElement(&elem, children_pos, i)
        if (elem == 1) {
            // create tmp vector with all zeroes
            GrB_Vector tmp = NULL;
            LAGr_Vector_new (&tmp, GrB_INT32, size);

            // extract the column corresponding to ith position
            LAGr_extract(tmp, GrB_NULL, GrB_NULL, step_1_matrix, GrB_ALL, size, i, GrB_NULL)

            // shift columns in matrix
            shift_cols_in_matrix(step_1_matrix, tmp_beta, i, 1, size);

            // put the column corresponding to ith position
            LAGr_assign(step_1_matrix, GrB_NULL, GrB_NULL, tmp, GrB_ALL, size, tmp_beta, GrB_NULL);

            // increment stack counter as we now work with s_u without this node
            tmp_beta++;

            // free after iter
            GrB_free(&tmp);
        }
    }


    // step 2
    GrB_Matrix step_2_matrix = NULL;
    build_identity_matrix(&step_2_matrix, size);

    // create tmp_matrix_1 for copying columns which have to be moved to top of s_p
    GrB_Matrix tmp_matrix_1 = NULL;
    GrB_Index tmp_matrix_1_ncols = tmp_beta - beta;
    LAGr_Matrix_new(&tmp_matrix_1, GrB_INT32, size, tmp_matrix_1_ncols);

    // create indexes of columns which have to be moved to top of s_p
    GrB_Index * col_moved_indexes = malloc(sizeof(GrB_Index) * tmp_matrix_1_ncols);
    for (int j = 0; j < tmp_matrix_1_ncols; j++) {
        col_moved_indexes[j] = beta + j;
    }

    // extract all columns which have to be moved to top of s_p
    LAGr_extract (tmp_matrix_1, NULL, NULL, step_2_matrix, GrB_ALL, size, col_moved_indexes, tmp_matrix_1_ncols, GrB_NULL);

    // shift columns in matrix
    shift_cols_in_matrix(step_2_matrix, tay, beta, tmp_matrix_1_ncols, size);

    // put the columns which have to be moved to top of s_p to the top of s_p
    for (int j = 0; j < tmp_matrix_1_ncols; j++) {
        col_moved_indexes[j] = tay + j;
    }
    LAGr_assign(step_2_matrix, GrB_NULL, GrB_NULL, tmp_matrix_1, GrB_ALL, size, col_moved_indexes, tmp_matrix_1_ncols, GrB_NULL);

    // summarize permutations by matrix multiplication
    LAGr_Matrix_new(perm_matrix, GrB_INT32, size, size);
    LAGr_mxm(*perm_matrix, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_INT32, step_1_matrix, step_2_matrix, GrB_NULL);

    // free
    free(col_moved_indexes);
    LAGr_free(&tmp_matrix_1);
    LAGr_free(&step_1_matrix);
    LAGr_free(&step_2_matrix);
}


GrB_Info get_vector_e(GrB_Vector *res, GrB_Index size, GrB_Index position) {
    LAGr_Vector_new (res, GrB_INT32, size);
    LAGr_assign (*res, NULL, NULL, 0, GrB_ALL, size, NULL);
    LAGr_Vector_setElement (*res, 1, position);
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
    LAGr_Matrix_nrows(&n, A);
    LAGr_Vector_new (&v, GrB_INT32, n);
    for (GrB_Index i = 0; i < n; i++) {
        LAGr_Vector_setElement(v, i, i);
    }

    int32_t children_amount;
    GrB_Index tay = 0, beta = 0;
    GrB_Matrix permutation_matrix = NULL;

    // descriptor for transposing left matrix in multiplication
    GrB_Descriptor left_transpose_descr = NULL;
    LAGr_Descriptor_new(&left_transpose_descr);
    LAGr_Descriptor_set(left_transpose_descr, GrB_INP0, GrB_TRAN);

    // push and pop root
    GrB_Vector children = NULL;
    get_vector_e(&children, n, source);
    build_permutation_matrix(&permutation_matrix, &children_amount, children, tay, beta, n);

    LAGr_mxv(v, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_INT32, permutation_matrix, v, left_transpose_descr);
    LAGr_mxm(A, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_INT32, permutation_matrix, A, left_transpose_descr);
    LAGr_mxm(A, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_INT32, A, permutation_matrix, GrB_NULL);

    beta += children_amount;

    // free before loop
    GrB_free(&children);
    GrB_free(&permutation_matrix);

    while (tay < n) {
        // get children
        get_vector_e(&children, n, tay);
        LAGr_mxv(children, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_INT32, A, children, left_transpose_descr);

        // pop
        tay++;
        if (beta < tay) {
            beta++;
        }

        // push
        build_permutation_matrix(&permutation_matrix, &children_amount, children, tay, beta, n);
        LAGr_mxv(v, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_INT32, permutation_matrix, v, left_transpose_descr);
        LAGr_mxm(A, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_INT32, permutation_matrix, A, left_transpose_descr);
        LAGr_mxm(A, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_INT32, A, permutation_matrix, GrB_NULL);
        beta += children_amount;

        // free before next iteration
        GrB_free(&children);
        GrB_free(&permutation_matrix);
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*v_output) = v;       // return result
    v = NULL;              // set to NULL so LAGRAPH_FREE_ALL doesn't free it
    GrB_free (&left_transpose_descr);
    return (GrB_SUCCESS);
}

