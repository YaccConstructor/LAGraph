//------------------------------------------------------------------------------
// dfsTest: test LAGraph_dfs_*.c
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

// Contributed by Dmitrii Iarosh, SPSU

// Usage:
//
// dfs_test matrixmarketfile.mtx
// dfs_test < matrixfile.mtx

#include "LAGraph.h"
#include <time.h>
// #include <sys/time.h>

#define LAGRAPH_FREE_ALL


double to_microsec(struct timespec t1, struct timespec t2)
{
    return
        (t2.tv_sec  - t1.tv_sec ) * 1e6 +
        (t2.tv_nsec - t1.tv_nsec) * 1e-3;
}

int find_root(GrB_Matrix matrix, GrB_Index size) {

//    GrB_Descriptor desc = NULL;
//    LAGr_Descriptor_new(&desc);
//    LAGr_Descriptor_set(desc, GrB_INP0, GrB_TRAN);
//    GrB_Vector sums = NULL;
//    LAGr_Vector_new(&sums, GrB_INT32, size);
//    // it doesn't work with transpose, I don't know why.
//    LAGr_reduce(sums, GrB_NULL, GrB_NULL, GxB_PLUS_INT32_MONOID, matrix, desc);

    // as reduce doesn't work I go straight
    for (GrB_Index j = 0; j < size; j++) {
        GrB_Index i;
        for (i = 0; i < size; i++) {
            int q = 0;
            GrB_Matrix_extractElement(&q, matrix, i, j);
            if (q != 0) {
                break;
            }
        }
        if (i == size) {
            return j;
        }
    }

    return 0;
}


int main (int argc, char **argv)
{
    GrB_Info info ;

    GrB_Matrix graph_matrix = NULL;
    GrB_Vector result = NULL;
    LAGRAPH_OK (LAGraph_init ( ));

    if (argc > 1) {
        char *filename = argv[1];
        printf("matrix will be read from: %s\n", filename);

        FILE *f = fopen(filename, "r");
        if (f == NULL) {
            printf("Matrix file not found: [%s]\n", filename);
            exit(1);
        }
        LAGRAPH_OK(LAGraph_mmread(&graph_matrix, f));
        fclose(f);
    }
    else
    {
        // Usage:  ./dfs_test < matrixfile.mtx
        printf("matrix: from stdin\n");

        // read in the file in Matrix Market format from stdin
        LAGRAPH_OK(LAGraph_mmread(&graph_matrix, stdin));
    }
    GrB_Index n;
    LAGr_Matrix_nrows(&n, graph_matrix);

    int root = find_root(graph_matrix, n);
    printf("Root: %d\n", root);

    struct timespec tstart={0,0}, tend={0,0};
    clock_gettime(CLOCK_MONOTONIC, &tstart);
    LAGraph_dfs_traverse_bin_tree(&result, graph_matrix, root, false);
    clock_gettime(CLOCK_MONOTONIC, &tend);
    double time = to_microsec(tstart, tend);
    printf("Time: %f microsec\n", time);

    int k;
    printf("Result: ");
    for (GrB_Index i = 0; i < n; i++) {
        GrB_Vector_extractElement(&k, result, i);
        printf("%d ", k);
    }
    printf("\n");

    GrB_free(&graph_matrix);
    GrB_free(&result);
    LAGRAPH_FREE_ALL;
    LAGRAPH_OK(GrB_finalize( ));
}

