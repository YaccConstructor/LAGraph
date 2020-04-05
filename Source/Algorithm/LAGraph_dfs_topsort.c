//------------------------------------------------------------------------------
// LAGraph_dfs_topsort: DFS-based topological sort
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

// LAGraph_dfs_topsort: DFS-based topological sort, returning ordered
// vector of vertices.
// Contributed by Vladislav Myasnikov, SPbU, Saint Petersburg, Russia

// In this algorithm vertices are classified into three different statuses:
//   1. Unordered. The order in which they are visited have not been determined.
// All vertices are assigned the unordered status at the start of the algorithm.
//   2. Partially ordered. An initial visit order for the vertex has been assigned,
// but this order may be changed later. A partially ordered vertex may be visited
// multiple times.
//   3. Ordered. The order in which vertex is processed is fixed. This status is
// assigned to a vertex when the vertex will not be visited again in the future.

// Vector s keeps track of the statuses of all vertices: dim(s) = N = |V|
// Each element of s is the unique label (0..N-1) associated with specific vertex.
// Vector s is partitioned into three sub-vectors as follows:
//
//                             |s_s| - ordered
//                      s ---> |s_p| - partially ordered
//                             |s_u| - unordered
//
// tau and beta variables are tops of s_p and s_u respectively

// LAGraph_dfs_topsort performs a single-source DFS, starting at a source node and
// terminating when all vertices are ordered. This algorithm use matrix and vector
// permutations via pre- and post-multiplying the matrices and vectors with
// a permutation matrix.

// For more information, see https://doi.org/10.1145/3315454.3329962

//------------------------------------------------------------------------------

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_ALL

/* Shift specific matrix columns */
GrB_Info shift_columns(
    GrB_Matrix M,           // source matrix
    GrB_Index start,        // index of 1st column to be shifted
    GrB_Index end,          // index of last column to be shifted
    GrB_Index shift_step,   // shift value
    GrB_Index dim           // dimension of columns
)
{
    if (start <= end && start >= 0 && end >= 0) {
        // TMP := Matrix[dim x (end - start + 1)]
        GrB_Matrix TMP = GrB_NULL;
        GrB_Index column_amount = end - start + 1;
        LAGr_Matrix_new(&TMP, GrB_UINT64, dim, column_amount);

        // create indices of columns to be shifted
        GrB_Index* indices = malloc(sizeof(GrB_Index) * column_amount);
        for (int j = 0; j < column_amount; j++)
            indices[j] = start + j;

        // TMP := Submatrix(M, i = 0..dim-1, j = start..end)
        LAGr_extract(TMP, NULL, NULL, M, GrB_ALL, dim, indices, column_amount, GrB_NULL);

        // shift column indices
        for (int j = 0; j < column_amount; j++)
            indices[j] += shift_step;

        // put specific matrix columns in a new place
        LAGr_assign(M, GrB_NULL, GrB_NULL, TMP, GrB_ALL, dim, indices, column_amount, GrB_NULL);

        // cleanup
        GrB_free(&TMP);
        free(indices);
    }
}

/* Create identity matrix P[dim x dim] */
GrB_Info build_id_permute(
    GrB_Matrix* P,   // permutation matrix to be created (1s are only on the main diagonal)
    GrB_Index dim    // dimension of permutation matrix
)
{
    LAGr_Matrix_new(P, GrB_UINT64, dim, dim);
    LAGr_assign(*P, GrB_NULL, GrB_NULL, 0, GrB_ALL, dim, GrB_ALL, dim, GrB_NULL);
    for (int i = 0; i < dim; i++)
        LAGr_Matrix_setElement(*P, 1, i, i);
}

/* Create permutation matrix P[dim x dim] */
GrB_Info build_permute(
    GrB_Matrix *P,   // permutation matrix to be created
    GrB_Index *v,    // amount of vertices to be moved from s_u to s_p
    GrB_Vector n,    // boolean vector (1s indicate current considered vertices)
    GrB_Index dim,   // dimension of permutation matrix
    GrB_Index tau,   // top of s_p
    GrB_Index beta   // top of s_u
)
{
    //------------------------------------------------------------------------------
    // Matrix P1 creation
    // Moving all considered vertices belonging to s_p to the top of s_p
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

            // tmp := i-th column of P1
            LAGr_extract(tmp, GrB_NULL, GrB_NULL, P1, GrB_ALL, dim, i, GrB_NULL);

            // shift columns in P1
            shift_columns(P1, tau, i - 1, 1, dim);

            // tmp_beta-th column of P1 := tmp
            LAGr_assign(P1, GrB_NULL, GrB_NULL, tmp, GrB_ALL, dim, tau, GrB_NULL);

            // cleanup
            GrB_free(&tmp);
        }
    }

    //------------------------------------------------------------------------------
    // Matrix P2 creation
    // Moving all considered vertices belonging to s_u to the bottom of s_p
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

            // tmp := i-th column of P2
            LAGr_extract(tmp, GrB_NULL, GrB_NULL, P2, GrB_ALL, dim, i, GrB_NULL);

            // shift columns in P2
            shift_columns(P2, tmp_beta, i - 1, 1, dim);

            // tmp_beta-th column of P2 := tmp
            LAGr_assign(P2, GrB_NULL, GrB_NULL, tmp, GrB_ALL, dim, tmp_beta, GrB_NULL);

            // move the temporary border between s_p and s_u
            tmp_beta++;

            // increment amount of vertices to be moved from s_u to s_p
            (*v)++;

            // cleanup
            GrB_free(&tmp);
        }
    }

    //------------------------------------------------------------------------------
    // Matrix P3 creation
    // Moving all shifted vertices in the previous step to the top of s_p
    //------------------------------------------------------------------------------

    GrB_Matrix P3 = GrB_NULL;
    build_id_permute(&P3, dim);

    // TMP := Matrix[dim x (tmp_beta - beta)]
    GrB_Matrix TMP = GrB_NULL;
    GrB_Index column_amount = tmp_beta - beta;
    LAGr_Matrix_new(&TMP, GrB_UINT64, dim, column_amount);

    // create indices of columns to be moved to the top of s_p
    GrB_Index* indices = malloc(sizeof(GrB_Index) * column_amount);
    for (int j = 0; j < column_amount; j++)
        indices[j] = beta + j;

    // TMP := Submatrix(M, i = 0..dim-1, j = beta..tmp_beta-1)
    LAGr_extract(TMP, GrB_NULL, GrB_NULL, P3, GrB_ALL, dim, indices, column_amount, GrB_NULL);

    // shift columns in P3
    shift_columns(P3, tau, beta - 1, column_amount, dim);

    // shift column indices
    for (int j = 0; j < column_amount; j++)
        indices[j] = tau + j;

    // put specific columns of P3 in a new place
    LAGr_assign(P3, GrB_NULL, GrB_NULL, TMP, GrB_ALL, dim, indices, column_amount, GrB_NULL);

    //------------------------------------------------------------------------------
    // Matrix P creation and cleanup
    //------------------------------------------------------------------------------

    // P := P1 * P2 * P3
    LAGr_Matrix_new(P, GrB_UINT64, dim, dim);
    LAGr_mxm(*P, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_UINT64, P1, P2, GrB_NULL);
    LAGr_mxm(*P, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_UINT64, *P, P3, GrB_NULL);

    // cleanup
    LAGr_free(&P1);
    LAGr_free(&P2);
    LAGr_free(&P3);
    LAGr_free(&TMP);
    free(indices);
}

/* Create a vector with a single unit and the remaining zeros */
GrB_Info build_vector_with_single_unit(
    GrB_Vector *v,       // vector to be created
    GrB_Index dim,       // dimension of vector
    GrB_Index position   // position of a unit
)
{
    LAGr_Vector_new(v, GrB_UINT64, dim);
    LAGr_assign(*v, NULL, NULL, 0, GrB_ALL, dim, NULL);
    LAGr_Vector_setElement(*v, 1, position);
}

/* DFS-based topological sort */
GrB_Info LAGraph_dfs_topsort
(
    GrB_Vector *order_output,   // topological order of vertices
    GrB_Matrix A,               // input graph, treated as if boolean in semiring
    GrB_Index source            // starting vertex of the DFS
)
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    GrB_Info info;

    if (A == GrB_NULL || order_output == GrB_NULL)
    {
        LAGRAPH_ERROR ("required arguments are NULL", GrB_NULL_POINTER);
    }

    GrB_Index nrows, ncols;
    LAGRAPH_OK (GrB_Matrix_nrows (&nrows, A));
    LAGRAPH_OK (GrB_Matrix_ncols (&ncols, A));

    if (nrows != ncols)
    {
        LAGRAPH_ERROR ("A must be square", GrB_INVALID_VALUE);
    }

    GrB_Index n = nrows;

    if (source >= n || source < 0)
    {
        LAGRAPH_ERROR ("invalid value for source vertex", GrB_INVALID_VALUE) ;
    }

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GrB_Index tau = 0;
    GrB_Index beta = 0;
    GrB_Index child_amount = 0;
    GrB_Matrix P = GrB_NULL;
    *order_output = GrB_NULL;

    // s := (0,1,2..,n-1)
    GrB_Vector s = GrB_NULL;
    LAGRAPH_OK (GrB_Vector_new(&s, GrB_UINT64, n));
    for (GrB_Index i = 0; i < n; i++)
        LAGRAPH_OK (GrB_Vector_setElement_UINT64(s, i, i));

    // children := (0,..,0)
    GrB_Vector children = GrB_NULL;
    LAGRAPH_OK (GrB_Vector_new(&children, GrB_UINT64, n));

    // create descriptor for transposing matrix in multiplication
    GrB_Descriptor tr_desc = GrB_NULL;
    LAGr_Descriptor_new(&tr_desc);
    LAGr_Descriptor_set(tr_desc, GrB_INP0, GrB_TRAN);

    //--------------------------------------------------------------------------
    // DFS-based topological sort - 1st step
    //--------------------------------------------------------------------------
    // Placing the source vertex at the top of s_p by using a vector permutation
    // and moving the bottom line.
    //--------------------------------------------------------------------------

    // children := (0,..,0,1,0,..,0), where 1 is in source-th position
    build_vector_with_single_unit(&children, n, source);

    // get permutation matrix and amount of vertices moved from s_u to s_p
    build_permute(&P, &child_amount, children, n, tau, beta);

    // A := TRANSPOSE(P) * A * P
    LAGr_mxm(A, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_UINT64, P, A, tr_desc);
    LAGr_mxm(A, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_UINT64, A, P, GrB_NULL);

    // s := TRANSPOSE(P) * s
    LAGr_mxv(s, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_UINT64, P, s, tr_desc);

    // move the border between s_p and s_u
    beta += child_amount;

    // cleanup
    GrB_free(&children);
    GrB_free(&P);

    //--------------------------------------------------------------------------
    // DFS-based topological sort - 2nd step
    //--------------------------------------------------------------------------
    // Until all the vertices are in the s_s:
    // 1. Visiting the top vertex of s_p.
    // 2. Pushing children of this vertex from s_u and s_p at the top of s_p.
    //    This operation is performed via permutation as in 1st step above.
    // 3. Popping the top vertex of s_p. This means moving forward the top line,
    //    which in the first iteration is the source, and adding the vertex to s_s.
    //    This operation is performed only when child_amount = 0, i.e.,
    //    it's a stock vertex or all of its children are already in s_s.
    //--------------------------------------------------------------------------

    while (tau < n)
    {
        // children := (0,..,0,1,0,..,0), where 1 is in tau-th position
        build_vector_with_single_unit(&children, n, tau);

        // children := TRANSPOSE(A) * children
        LAGr_mxv(children, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_UINT64, A, children, tr_desc);

        // get permutation matrix and amount of vertices moved from s_u to s_p
        build_permute(&P, &child_amount, children, n, tau, beta);

        // A := TRANSPOSE(P) * A * P
        LAGr_mxm(A, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_UINT64, P, A, tr_desc);
        LAGr_mxm(A, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_UINT64, A, P, GrB_NULL);

        // s := TRANSPOSE(P) * s
        LAGr_mxv(s, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_UINT64, P, s, tr_desc);

        // move the border between s_p and s_u
        beta += child_amount;

        // move the border between s_s and s_p (keep beta >= tau)
        if (child_amount == 0) tau++;
        if (beta < tau) beta++;

        // cleanup
        GrB_free(&children);
        GrB_free(&P);
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*order_output) = s;
    s = GrB_NULL;
    GrB_free(&tr_desc);
    return (GrB_SUCCESS);
}