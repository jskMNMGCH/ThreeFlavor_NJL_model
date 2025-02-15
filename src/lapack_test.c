#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>

int main() {
    // 3x3の行列
    double A[9] = {
        1.0, 2.0, 3.0,
        0.0, 4.0, 5.0,
        1.0, 0.0, 6.0
    };

    // 逆行列の結果を格納するための行列
    double A_inv[9];  // 逆行列

    // 行列のサイズ
    int N = 3;
    int lda = N;  // 行列Aのリーディングディメンション
    int pivots[3]; // ピボットインデックス（LU分解用）
    int info;

    // LU分解 (dgetrf)
    // dgetrf( int *N, int *N, double *A, int *lda, int *ipiv, int *info )
    // A: 行列A、lda: リーディングディメンション、ipiv: ピボットインデックス、info: 実行結果
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, N, N, A, lda, pivots);
    if (info != 0) {
        printf("LU分解に失敗しました。エラーコード: %d\n", info);
        return 1;
    }

    // 逆行列の計算 (dgetri)
    // dgetri( int matrix_layout, int n, double *a, int lda, int *ipiv )
    // A: LU分解後の行列、lda: リーディングディメンション、ipiv: ピボットインデックス、info: 実行結果
    // 逆行列をA_invに保存
    for (int i = 0; i < N * N; i++) {
        A_inv[i] = A[i];  // AをA_invにコピー
    }

    info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, N, A_inv, lda, pivots);
    if (info != 0) {
        printf("逆行列計算に失敗しました。エラーコード: %d\n", info);
        return 1;
    }

    // 結果を表示
    printf("逆行列は以下の通りです:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%lf ", A_inv[i * N + j]);
        }
        printf("\n");
    }

    return 0;
}

/*
lapacke.hなどの場所がわからない時は，
find /opt/homebrew -name lapacke.h

gcc -std=c18 -o lapack_test lapack_test.c -I/opt/homebrew/Cellar/openblas/0.3.25/include -L/opt/homebrew/Cellar/openblas/0.3.25/lib -llapack -lblas -lm
./lapack_test

 */

