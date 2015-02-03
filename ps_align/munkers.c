/**
 * Implementation of the Hungarian algorithm based on description
 * given at http://csclab.murraystate.edu/bob.pilgrim/445/munkres.html
 *
 * Compile with __STANDALONE__ to include the main and make a standalone
 * program. Compiel with __CHATTY__ to pring out step by step Cost and
 * Mask matrices per iteration.
 */
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

/////////////////////////////////////////////////////////////////////////////////////////////////////
//                           F U N C T I O N   P R O T O T Y P E S				   //
void step_one(int** matrix, int nrow, int ncol, int& step);
void step_two(int** matrix, int nrow, int ncol, int* rowCover, int* colCover, int** m, int& step);
void step_three(int nrow, int ncol, int* colCover, int** m, int& step);
void find_a_zero(int** matrix, int nrow, int ncol, int* rowCover, int* colCover, int& row, int& col);
bool star_in_row(int** m, int ncol, int row);
void find_star_in_row(int** m, int ncol, int row, int& col);
void step_four(int** matrix, int nrow, int ncol, int* rowCover, int* colCover, int** m,
               int& path_row_0, int& path_col_0, int& step);
void find_star_in_col(int nrow, int** m, int c, int& r);
void find_prime_in_row(int nrow, int** m, int r, int& c);
void augment_path(int** m, int path_count, int** path);
void clear_covers(int nrow, int ncol, int* rowCover, int* colCover);
void erase_primes(int nrow, int ncol, int** m);
void step_five(int nrow, int ncol, int* rowCover, int* colCover, int** m, int& path_row_0,
               int& path_col_0, int& path_count, int** path, int& step);
void find_smallest(int** matrix, int nrow, int ncol, int* rowCover, int* colCover,
                   int& minval) ;
void step_six(int** matrix, int nrow, int ncol, int* rowCover, int* colCover, int& step);
int** runMunkers(int** matrix, int nrow, int ncol, bool max);
/////////////////////////////////////////////////////////////////////////////////////////////////////

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

void step_one(int** matrix, int nrow, int ncol, int& step) {
    int min_in_row;

    for (int r = 0; r < nrow; r++) {
        min_in_row = matrix[r][0];
        for (int c = 0; c < ncol; c++)
            min_in_row = MIN(min_in_row, matrix[r][c]);
        for (int c = 0; c < ncol; c++) matrix[r][c] -= min_in_row;
    }
    step = 2;
}

void step_two(int** matrix, int nrow, int ncol, int* rowCover, int* colCover, int** m,
              int& step) {
    for (int r = 0; r < nrow; r++)
        for (int c = 0; c < ncol; c++) {
            if (matrix[r][c] == 0 && rowCover[r] == 0 && colCover[c] == 0) {
                m[r][c] = 1;
                rowCover[r] = 1;
                colCover[c] = 1;
            }
        }
    for (int r = 0; r < nrow; r++) rowCover[r] = 0;
    for (int c = 0; c < ncol; c++) colCover[c] = 0;
    step = 3;
}

void step_three(int nrow, int ncol, int* colCover, int** m, int& step) {
    int colcount;
    for (int r = 0; r < nrow; r++)
        for (int c = 0; c < ncol; c++)
            if (m[r][c] == 1)
                colCover[c] = 1;

    colcount = 0;
    for (int c = 0; c < ncol; c++)
        if (colCover[c] == 1)
            colcount += 1;

    step = (colcount >= nrow || colcount >= ncol) ? 7 : 4;
}

void find_a_zero(int** matrix, int nrow, int ncol, int* rowCover, int* colCover, int& row,
                 int& col) {
    int r = 0;
    int c;
    bool done;
    row = -1;
    col = -1;
    done = false;
    while (!done) {
        c = 0;
        while (true) {
            if (matrix[r][c] == 0 && rowCover[r] == 0 && colCover[c] == 0) {
                row = r;
                col = c;
                done = true;
            }
            c += 1;
            if (c >= ncol || done)
                break;
        }
        r += 1;
        if (r >= nrow)
            done = true;
    }
}

bool star_in_row(int** m, int ncol, int row) {
    for (int c = 0; c < ncol; c++)
        if (m[row][c] == 1)
            return true;
    return false;
}

void find_star_in_row(int** m, int ncol, int row, int& col) {
    col = -1;
    for (int c = 0; c < ncol; c++)
        if (m[row][c] == 1)
            col = c;
}

void step_four(int** matrix, int nrow, int ncol, int* rowCover, int* colCover, int** m,
               int& path_row_0, int& path_col_0, int& step) {
    int row = -1;
    int col = -1;
    bool done;

    done = false;
    while (!done) {
        find_a_zero(matrix, nrow, ncol, rowCover, colCover, row, col);
        if (row == -1) {
            done = true;
            step = 6;
        } else {
            m[row][col] = 2;
            if (star_in_row(m, ncol, row)) {
                find_star_in_row(m, ncol, row, col);
                rowCover[row] = 1;
                colCover[col] = 0;
            } else {
                done = true;
                step = 5;
                path_row_0 = row;
                path_col_0 = col;
            }
        }
    }
}

void find_star_in_col(int nrow, int** m, int c, int& r) {
    r = -1;
    for (int i = 0; i < nrow; i++)
        if (m[i][c] == 1)
            r = i;
}

void find_prime_in_row(int nrow, int** m, int r, int& c) {
    for (int j = 0; j < nrow; j++)
        if (m[r][j] == 2)
            c = j;
}

void augment_path(int** m, int path_count, int** path) {
    for (int p = 0; p < path_count; p++)
        if (m[path[p][0]][path[p][1]] == 1)
            m[path[p][0]][path[p][1]] = 0;
        else
            m[path[p][0]][path[p][1]] = 1;
}

void clear_covers(int nrow, int ncol, int* rowCover, int* colCover) {
    for (int r = 0; r < nrow; r++)
        rowCover[r] = 0;
    for (int c = 0; c < ncol; c++)
        colCover[c] = 0;
}

void erase_primes(int nrow, int ncol, int** m) {
    for (int r = 0; r < nrow; r++)
        for (int c = 0; c < ncol; c++)
            if (m[r][c] == 2)
                m[r][c] = 0;
}

void step_five(int nrow, int ncol, int* rowCover, int* colCover, int** m, int& path_row_0,
               int& path_col_0, int& path_count, int** path, int& step) {
    bool done;
    int r = -1;
    int c = -1;

    path_count = 1;
    path[path_count - 1][0] = path_row_0;
    path[path_count - 1][1] = path_col_0;
    done = false;
    while (!done) {
        find_star_in_col(nrow, m, path[path_count - 1][1], r);
        if (r > -1) {
            path_count += 1;
            path[path_count - 1][0] = r;
            path[path_count - 1][1] = path[path_count - 2][1];
        } else
            done = true;
        if (!done) {
            find_prime_in_row(nrow, m, path[path_count - 1][0], c);
            path_count += 1;
            path[path_count - 1][0] = path[path_count - 2][0];
            path[path_count - 1][1] = c;
        }
    }
    augment_path(m, path_count, path);
    clear_covers(nrow, ncol, rowCover, colCover);
    erase_primes(nrow, ncol, m);
    step = 3;
}

void find_smallest(int** matrix, int nrow, int ncol, int* rowCover, int* colCover,
                   int& minval) {
    for (int r = 0; r < nrow; r++)
        for (int c = 0; c < ncol; c++)
            if (rowCover[r] == 0 && colCover[c] == 0)
                minval = MIN(minval, matrix[r][c]);
}

void step_six(int** matrix, int nrow, int ncol, int* rowCover, int* colCover, int& step) {
    int minval = INT_MAX;
    find_smallest(matrix, nrow, ncol, rowCover, colCover, minval);
    for (int r = 0; r < nrow; r++)
        for (int c = 0; c < ncol; c++) {
            if (rowCover[r] == 1)
                matrix[r][c] += minval;
            if (colCover[c] == 0)
                matrix[r][c] -= minval;
        }
    step = 4;
}

int** runMunkers(int** matrix, int nrow, int ncol, bool max) {

    if (max == true) {
        int maxValue = matrix[0][0];
        for (int i = 0; i < nrow; i++)
            for (int j = 0; j < ncol; j++)
                maxValue = MAX(maxValue, matrix[i][j]);

        for (int i = 0; i < nrow; i++)
            for (int j = 0; j < ncol; j++)
                matrix[i][j] = maxValue - matrix[i][j];
    }

    bool done = false;
    int step = 1;
    int rowCover[nrow];
    int colCover[ncol];
    int **m = (int**)malloc(nrow*(sizeof(int*)));
    int path_row_0 = 0;
    int path_col_0 = 0;
    int path_count = 0;
    for (int i = 0; i < nrow; i++) {
        rowCover[i] = colCover[i] = 0;
        m[i] = (int*)malloc(sizeof(int)*ncol);
        for (int j = 0; j < ncol; j++) m[i][j] = 0;
    }
    int** path = (int**)malloc((nrow + ncol)*(sizeof(int*)));
    for (int i = 0; i < nrow + ncol; i++) path[i] = (int*)malloc(sizeof(int)*ncol);
    while (!done) {
#ifdef __CHATTY__
        printf("--cost\n");
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) printf("%d ", matrix[i][j]);
            printf("\n");
        }
        printf("--mask\n");
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) printf("%d ", m[i][j]);
            printf("\n");
        }
#endif
        switch (step) {
        case 1:
            step_one(matrix, nrow, ncol, step);
            break;
        case 2:
            step_two(matrix, nrow, ncol, rowCover, colCover, m, step);
            break;
        case 3:
            step_three(nrow, ncol, colCover, m, step);
            break;
        case 4:
            step_four(matrix, nrow, ncol, rowCover, colCover, m, path_row_0, path_col_0,
                      step);
            break;
        case 5:
            step_five(nrow, ncol, rowCover, colCover, m, path_row_0, path_col_0,
                      path_count, path, step);
            break;
        case 6:
            step_six(matrix, nrow, ncol, rowCover, colCover, step);
            break;
        case 7:
            done = true;
            break;
        }
    }

    return m;
}

#ifdef __STANDALONE__
int main() {
    int row = 3, col = 3;
    int** c = (int**)malloc(row*(sizeof(int*)));
    for (int i = 0; i < row; i++) {
        c[i] = (int*)malloc(sizeof(int)*col);
        for (int j = 0; j < col; j++)
            c[i][j] = (i + 1) * (j + 1);
    }
    //
    c = runMunkers(c, row, col, true);
    printf("--cost:\n");
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++)
            printf("%d ", c[i][j]);
        printf("\n");
    }
    return 0;
}
#endif
