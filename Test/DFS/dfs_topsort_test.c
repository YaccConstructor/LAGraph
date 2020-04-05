//------------------------------------------------------------------------------
// dfs_topsort_test: read in (or create) a matrix and test DFS-based topological sort
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

// Contributed by Vladislav Myasnikov, SPbU, Saint Petersburg, Russia

//------------------------------------------------------------------------------

#include "LAGraph.h"

#define LAGRAPH_FREE_ALL    \
{                           \
    GrB_free (&result);     \
    GrB_free (&A);          \
}

int main(int argc, char **argv)
{
    GrB_Info info;

    GrB_Matrix A = NULL;
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
        LAGRAPH_OK(LAGraph_mmread(&A, f));
        fclose(f);
    }
    else
    {
        // Usage:  ./dfs_topsort_test < matrixfile.mtx
        printf("matrix: from stdin\n");

        // read in the file in Matrix Market format from stdin
        LAGRAPH_OK(LAGraph_mmread(&A, stdin));
    }

    GrB_Index n;
    LAGr_Matrix_nrows(&n, A);

    LAGraph_dfs_topsort(&result, A, 3);

    int k;
    printf("Topologically sorted array of vertices:\n");
    for (GrB_Index i = 0; i < n; i++) {
        GrB_Vector_extractElement(&k, result, i);
        printf("%d ", k + 1);
    }
    printf("\n");

    LAGRAPH_FREE_ALL;
    LAGRAPH_OK(GrB_finalize( ));
}