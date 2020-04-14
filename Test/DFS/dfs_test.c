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

#include "LAGraph.h"
#include <time.h>

#define LAGRAPH_FREE_ALL

#define OK(method)                                                          \
{                                                                           \
    GrB_Info this_info = method ;                                           \
    if (! (this_info == GrB_SUCCESS || this_info == GrB_NO_VALUE))          \
    {                                                                       \
        printf ("Test failed: [%d] %s\n", this_info, GrB_error ( )) ;     \
        LAGRAPH_FREE_ALL ;                                                  \
        return (this_info) ;                                                \
    }                                                                       \
}

#define GOLDEN_PATH "../test.golden.txt"
#define TEST_PATH "../test.mtx"


double to_microsec(struct timespec t1, struct timespec t2) {
    return
            (t2.tv_sec - t1.tv_sec) * 1e6 +
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
            OK(GrB_Matrix_extractElement(&q, matrix, i, j));
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

int str_to_int(char *str, int *ptr) {
    int num = 0;
    int i;
    for (i = *ptr; i < strlen(str); i++) {
        if (str[i] > '9' || str[i] < '0') {
            break;
        }
        num = num * 10 + str[i] - '0';
    }
    *ptr = ++i;
    return num;
}


GrB_Info dfs_test(GrB_Matrix matrix, GrB_Index size) {
    GrB_Vector result = NULL;
    int root = find_root(matrix, size);
    GrB_Matrix matrix1 = NULL;
    LAGr_Matrix_dup(&matrix1, matrix);
//    printf("Root: %d\n", root);

    struct timespec tstart = {0, 0}, tend = {0, 0};
    clock_gettime(CLOCK_MONOTONIC, &tstart);
    OK(LAGraph_dfs_traverse_bin_tree(&result, matrix, root, false));
    clock_gettime(CLOCK_MONOTONIC, &tend);
    double time = to_microsec(tstart, tend);
    printf("Pre-order time: %f microsec\n", time);

    char str[256];
    FILE *mf = fopen(GOLDEN_PATH, "r");

    if (mf == NULL) {
        printf("File not found\n");
        return GrB_INVALID_OBJECT;
    }

    int res_val = -1;
    fgets(str, sizeof(str), mf);
    int ptr = 0;
    for (int i = 0; i < size; i++) {
        int check_val = str_to_int(str, &ptr);
        res_val = -1;
        OK(GrB_Vector_extractElement(&res_val, result, i));
        if (res_val != check_val) {
            printf("Not equal %d != %d at position %d\n", res_val, check_val, i);
            return GrB_INVALID_VALUE;
        }
    }

    GrB_free(&result);

    clock_gettime(CLOCK_MONOTONIC, &tstart);
    OK(LAGraph_dfs_traverse_bin_tree(&result, matrix1, root, true));
    clock_gettime(CLOCK_MONOTONIC, &tend);
    time = to_microsec(tstart, tend);
    printf("Post-order time: %f microsec\n", time);

    fgets(str, sizeof(str), mf);
    ptr = 0;
    for (int i = 0; i < size; i++) {
        int check_val = str_to_int(str, &ptr);
        res_val = -1;
        OK(GrB_Vector_extractElement(&res_val, result, i));
        if (res_val != check_val) {
            printf("Not equal %d != %d at position %d\n", res_val, check_val, i);
            return GrB_INVALID_VALUE;
        }
    }

//    int k;
//    printf("Result: ");
//    for (GrB_Index i = 0; i < size; i++) {
//        GrB_Vector_extractElement(&k, result, i);
//        printf("%d ", k);
//    }
//    printf("\n");
    GrB_free(&result);
    GrB_free(&matrix1);
    fclose(mf);
    return GrB_SUCCESS;
}


int main(int argc, char **argv) {
    GrB_Info info;

    GrB_Matrix graph_matrix = NULL;
    LAGRAPH_OK (LAGraph_init());

    char *filename = TEST_PATH;
    printf("matrix will be read from: %s\n", filename);

    FILE *f = fopen(filename, "r");
    if (f == NULL) {
        printf("Matrix file not found: [%s]\n", filename);
        exit(1);
    }
    LAGRAPH_OK(LAGraph_mmread(&graph_matrix, f));
    fclose(f);

    GrB_Index n;
    LAGr_Matrix_nrows(&n, graph_matrix);

    if (dfs_test(graph_matrix, n) == GrB_SUCCESS) {
        printf("Test passed!\n");
    } else {
        printf("Test failed!\n");
    }

    GrB_free(&graph_matrix);
    LAGRAPH_FREE_ALL;
    LAGRAPH_OK(GrB_finalize());
}

