//------------------------------------------------------------------------------
// LAGraph_dfs_topsort: Topological sorting
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

/* Put comments here */

//------------------------------------------------------------------------------

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_ALL

void shift_cols_in_matrix(GrB_Matrix M, GrB_Index start, GrB_Index end, GrB_Index shift_step, GrB_Index dim)
{
    // Create matrix TMP[dim x (end - start)]
    GrB_Matrix TMP = GrB_NULL;
    GrB_Index column_amount = end - start;
    LAGr_Matrix_new(&TMP, GrB_UINT64, dim, column_amount);

    // create indexes of columns which will be shifted
    GrB_Index* col_indexes = malloc(sizeof(GrB_Index) * column_amount);
    for (int j = 0; j < column_amount; j++)
        col_indexes[j] = start + j;

    // extract all columns from tmp_beta column to ith column for shifting columns
    LAGr_extract(TMP, NULL, NULL, M, GrB_ALL, dim, col_indexes, column_amount, GrB_NULL);

    // shift columns in matrix
    for (int j = 0; j < column_amount; j++)
        col_indexes[j] += shift_step;

    LAGr_assign(M, GrB_NULL, GrB_NULL, TMP, GrB_ALL, dim, col_indexes, column_amount, GrB_NULL);

    // Cleanup
    // *death* <--- col_indexes, TMP
    free(col_indexes);
    GrB_free(&TMP);
}

/*
 * Create identity matrix P[dim x dim]
 */
void build_id_permute(GrB_Matrix* P, GrB_Index dim)
{
    LAGr_Matrix_new(P, GrB_UINT64, dim, dim);
    LAGr_assign(*P, GrB_NULL, GrB_NULL, 0, GrB_ALL, dim, GrB_ALL, dim, GrB_NULL);
    for (int i = 0; i < dim; i++)
        LAGr_Matrix_setElement(*P, 1, i, i);
}

/*
 * Create permutation matrix P[dim x dim]
 * n: boolean vector where 1 indicate child of current node
 * v: amount of nodes moved from s_u to s_p
 */
void build_permute(GrB_Matrix *P, GrB_Index *v, GrB_Vector n, GrB_Index dim, GrB_Index tau, GrB_Index beta)
{
    //------------------------------------------------------------------------------
    // Matrix P1 creation
    // Moving from s_p to top of s_p
    //------------------------------------------------------------------------------

    GrB_Matrix P1 = GrB_NULL;
    build_id_permute(&P1, dim);

    for (GrB_Index i = 0; i < dim; i++)
    {
        uint64_t x;
        LAGr_Vector_extractElement(&x, n, i);

        if (x == 1 && i >= tau && i < beta)
        {
            // tmp := (0,..,0)
            GrB_Vector tmp = GrB_NULL;
            LAGr_Vector_new(&tmp, GrB_UINT64, dim);

            // tmp <--- i-th column of P1
            LAGr_extract(tmp, GrB_NULL, GrB_NULL, P1, GrB_ALL, dim, i, GrB_NULL);

            // shift columns in P1
            shift_cols_in_matrix(P1, tau, i, 1, dim);

            // tmp_beta-th column of P1 <--- tmp
            LAGr_assign(P1, GrB_NULL, GrB_NULL, tmp, GrB_ALL, dim, tau, GrB_NULL);

            // *death* <--- tmp
            GrB_free(&tmp);
        }
    }

    //------------------------------------------------------------------------------
    // Matrix P2 creation
    // Moving from s_u to bottom of s_p
    //------------------------------------------------------------------------------

    GrB_Matrix P2 = GrB_NULL;
    build_id_permute(&P2, dim);
    GrB_Index tmp_beta = beta;
    *v = 0;

    for (GrB_Index i = 0; i < dim; i++)
    {
        uint64_t x;
        LAGr_Vector_extractElement(&x, n, i);

        if (x == 1 && i >= beta)
        {
            // tmp := (0,..,0)
            GrB_Vector tmp = GrB_NULL;
            LAGr_Vector_new(&tmp, GrB_UINT64, dim);

            // tmp <--- i-th column of P2
            LAGr_extract(tmp, GrB_NULL, GrB_NULL, P2, GrB_ALL, dim, i, GrB_NULL);

            // shift columns in P2
            shift_cols_in_matrix(P2, tmp_beta, i, 1, dim);

            // tmp_beta-th column of P2 <--- tmp
            LAGr_assign(P2, GrB_NULL, GrB_NULL, tmp, GrB_ALL, dim, tmp_beta, GrB_NULL);

            // increment stack counter as we now work with s_u without this node
            tmp_beta++;

            (*v)++;

            // *death* <--- tmp
            GrB_free(&tmp);
        }
    }

    //------------------------------------------------------------------------------
    // Matrix P3 creation
    // Moving from bottom of s_p to top of s_p
    //------------------------------------------------------------------------------

    GrB_Matrix P3 = GrB_NULL;
    build_id_permute(&P3, dim);

    // create TMP for copying columns which have to be moved to top of s_p
    GrB_Matrix TMP = GrB_NULL;
    GrB_Index column_amount = tmp_beta - beta;
    LAGr_Matrix_new(&TMP, GrB_UINT64, dim, column_amount);

    // create indexes of columns which have to be moved to top of s_p
    GrB_Index* col_moved_indexes = malloc(sizeof(GrB_Index) * column_amount);
    for (int j = 0; j < column_amount; j++)
        col_moved_indexes[j] = beta + j;

    // extract all columns which have to be moved to top of s_p
    LAGr_extract(TMP, GrB_NULL, GrB_NULL, P3, GrB_ALL, dim, col_moved_indexes, column_amount, GrB_NULL);

    // shift columns in matrix
    shift_cols_in_matrix(P3, tau, beta, column_amount, dim);

    // put the columns which have to be moved to top of s_p to the top of s_p
    for (int j = 0; j < column_amount; j++)
        col_moved_indexes[j] = tau + j;

    LAGr_assign(P3, GrB_NULL, GrB_NULL, TMP, GrB_ALL, dim, col_moved_indexes, column_amount, GrB_NULL);

    //------------------------------------------------------------------------------
    // Matrix P creation
    // P <--- P1 * P2 * P3
    //------------------------------------------------------------------------------

    LAGr_Matrix_new(P, GrB_UINT64, dim, dim);
    LAGr_mxm(*P, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_UINT64, P1, P2, GrB_NULL);
    LAGr_mxm(*P, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_UINT64, *P, P3, GrB_NULL);

    //------------------------------------------------------------------------------
    // Cleanup
    // *death* <--- P1, P2, P3, TMP, col_moved_indexes
    //------------------------------------------------------------------------------

    LAGr_free(&P1);
    LAGr_free(&P2);
    LAGr_free(&P3);
    LAGr_free(&TMP);
    free(col_moved_indexes);
}

void get_vector_e(GrB_Vector *res, GrB_Index size, GrB_Index position)
{
    LAGr_Vector_new(res, GrB_UINT64, size);
    LAGr_assign(*res, NULL, NULL, 0, GrB_ALL, size, NULL);
    LAGr_Vector_setElement(*res, 1, position);
}

// temporary function
void print_vector(GrB_Vector s, GrB_Index dim, GrB_Index tau, GrB_Index beta)
{
    for (int i = 0; i < dim; i++)
    {
        u_int64_t x;
        LAGr_Vector_extractElement(&x, s, i);
        if (i == tau || i == beta) printf("|");
        else printf(" ");
        printf("%ld", x + 1);
    }
    printf("\n");
}

/* Put comments here */
GrB_Info LAGraph_dfs_topsort
(
    GrB_Vector *order_output,
    GrB_Matrix A,
    GrB_Index source
)
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    GrB_Info info;
    GrB_Index nrows, ncols;
    GrB_Vector s = GrB_NULL;

    if (A == GrB_NULL || order_output == GrB_NULL)
    {
        // required argument is missing
        LAGRAPH_ERROR ("required arguments are NULL", GrB_NULL_POINTER);
    }

    *order_output = GrB_NULL;
    LAGRAPH_OK (GrB_Matrix_nrows (&nrows, A));
    LAGRAPH_OK (GrB_Matrix_ncols (&ncols, A));
    if (nrows != ncols)
    {
        // A must be square
        LAGRAPH_ERROR ("A must be square", GrB_INVALID_VALUE);
    }
    GrB_Index n = nrows;

    if (source >= n || source < 0)
    {
        LAGRAPH_ERROR ("invalid value for source vertex s", GrB_INVALID_VALUE) ;
    }

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GrB_Index tau = 0;
    GrB_Index beta = 0;
    GrB_Index child_amount = 0;

    LAGRAPH_OK (GrB_Vector_new(&s, GrB_UINT64, n));
    for (GrB_Index i = 0; i < n; i++)
        LAGRAPH_OK (GrB_Vector_setElement_UINT64(s, i, i));

    GrB_Vector children = GrB_NULL;
    LAGRAPH_OK (GrB_Vector_new(&children, GrB_UINT64, n));

    GrB_Matrix P = GrB_NULL;

    // descriptor for transposing matrix in multiplication
    GrB_Descriptor tr_desc = GrB_NULL;
    LAGr_Descriptor_new(&tr_desc);
    LAGr_Descriptor_set(tr_desc, GrB_INP0, GrB_TRAN);

    //--------------------------------------------------------------------------
    // DFS traversal and topological sorting the nodes
    //--------------------------------------------------------------------------
    // 1st step
    // Root processing...
    //--------------------------------------------------------------------------

    get_vector_e(&children, n, source);

    // child_amount, P <--- *build_permute* <--- children, tau, beta
    build_permute(&P, &child_amount, children, n, tau, beta);

    // A <--- TRANSPOSE(P) * A * P
    LAGr_mxm(A, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_UINT64, P, A, tr_desc);
    LAGr_mxm(A, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_UINT64, A, P, GrB_NULL);

    // s <--- TRANSPOSE(P) * s
    LAGr_mxv(s, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_UINT64, P, s, tr_desc);

    // move the border between s_p and s_u
    beta += child_amount;

    // *death* <--- children, P
    GrB_free(&children);
    GrB_free(&P);

//    print_vector(s, n, tau, beta);

    //--------------------------------------------------------------------------
    // DFS traversal and topological sorting the nodes
    //--------------------------------------------------------------------------
    // 2nd step
    // Loop...
    //--------------------------------------------------------------------------

    while (tau < n)
    {
        get_vector_e(&children, n, tau);

        // children <--- TRANSPOSE(A) * children
        LAGr_mxv(children, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_UINT64, A, children, tr_desc);

        // children, P <--- *build_permute* <--- children, tau, beta
        build_permute(&P, &child_amount, children, n, tau, beta);

        // A <--- TRANSPOSE(P) * A * P
        LAGr_mxm(A, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_UINT64, P, A, tr_desc);
        LAGr_mxm(A, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_UINT64, A, P, GrB_NULL);

        // s <--- TRANSPOSE(P) * s
        LAGr_mxv(s, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_UINT64, P, s, tr_desc);

        // move the border between s_p and s_u
        beta += child_amount;

        // move the border between s_s and s_p (keep beta >= tau)
        if (child_amount == 0) tau++;
        if (beta < tau) beta++;

        // *death* <--- children, P
        GrB_free(&children);
        GrB_free(&P);

//        print_vector(s, n, tau, beta);
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*order_output) = s;
    s = GrB_NULL;
    GrB_free(&tr_desc);
    return (GrB_SUCCESS);
}