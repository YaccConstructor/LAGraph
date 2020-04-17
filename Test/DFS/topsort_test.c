//------------------------------------------------------------------------------
// topsort_test: read a matrix and test DFS-based topological sort
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
    GrB_free (&A);          \
}

GrB_Info run_test(GrB_Matrix A)
{
    int answer[5][5] = {
            { 1, 0, 2, 4, 3 },
            { 1, 0, 4, 2, 3 },
            { 0, 1, 2, 4, 3 },
            { 0, 1, 4, 2, 3 },
            { 0, 2, 1, 4, 3 }
    };

    GrB_Index n;
    LAGr_Matrix_nrows(&n, A);

    GrB_Vector order = GrB_NULL;
    LAGraph_dfs_topsort(&order, A, 0);

    int isSuccess = 0;
    for (int i = 0; i < 5; i++) {
        int j = 0;
        for (j = 0; j < n; j++) {
            int k;
            GrB_Vector_extractElement(&k, order, j);
            if (k != answer[i][j]) break;
        }
        if (j == n) isSuccess = 1;
    }

    if (isSuccess) printf("Test passed\n");
    else printf("Test failed\n");
}

int main(int argc, char **argv)
{
    GrB_Info info;

    GrB_Matrix A = GrB_NULL;
    LAGRAPH_OK (LAGraph_init ( ));

    // load matrix from test file
    char *filename = "test.mtx";
    printf("Matrix will be read from: %s\n", filename);

    FILE *f = fopen(filename, "r");
    if (f == NULL) {
        printf("Matrix file not found: %s\n", filename);
        exit(1);
    }
    LAGRAPH_OK(LAGraph_mmread(&A, f));
    fclose(f);

    run_test(A);

    LAGRAPH_FREE_ALL;
    LAGRAPH_OK(GrB_finalize( ));
}
