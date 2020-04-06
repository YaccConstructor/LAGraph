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
// #include <sys/time.h>

#define LAGRAPH_FREE_ALL    \
{                           \
    GrB_free (&result);    \
    GrB_free (&graph_matrix);         \
}

/*
double to_sec(struct timeval t1, struct timeval t2)
{
    return
        (t2.tv_sec  - t1.tv_sec ) + 
        (t2.tv_usec - t1.tv_usec) * 1e-6;
}
*/

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

    LAGraph_dfs_traverse_bin_tree(&result, graph_matrix, 4, true);

    GrB_Index n;
    LAGr_Matrix_nrows(&n, graph_matrix);
    int k;
    printf("Result: ");
    for (GrB_Index i = 0; i < n; i++) {
        GrB_Vector_extractElement(&k, result, i);
        printf("%d ", k);
    }
    printf("\n");

    LAGRAPH_FREE_ALL;
    LAGRAPH_OK(GrB_finalize( ));
}

