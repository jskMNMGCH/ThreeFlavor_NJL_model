#include <stdio.h>
#include <math.h>
#include <gsl/gsl_multiroots.h>

// パラメータを格納する構造体
typedef struct {
    double a;
    double b;
    double c;
} Params;

// 非線形連立方程式の定義
int func(const gsl_vector *v, void *params, gsl_vector *f) {
    double x = gsl_vector_get(v, 0);  // x
    double y = gsl_vector_get(v, 1);  // y
    double z = gsl_vector_get(v, 2);  // z

    Params *p = (Params *) params;  // パラメータを取得
    double a = p->a;
    double b = p->b;
    double c = p->c;

    gsl_vector_set(f, 0, x*x + y*y + a*z - 4);  // f1(x, y, z)
    gsl_vector_set(f, 1, x + b*y - 1);          // f2(x, y, z)
    gsl_vector_set(f, 2, c*x - z);              // f3(x, y, z)

    return GSL_SUCCESS;
}

int main() {
    // パラメータを設定
    Params params = {1.0, 2.0, 0.5};  // a = 1.0, b = 2.0, c = 0.5

    // 初期値の設定
    gsl_vector *x = gsl_vector_alloc(3);
    gsl_vector_set(x, 0, 1.0);  // x = 1.0
    gsl_vector_set(x, 1, 1.0);  // y = 1.0
    gsl_vector_set(x, 2, 1.0);  // z = 1.0

    // ソルバーの設定
    gsl_multiroot_function f = {&func, 3, &params};  // 3つの変数, パラメータ付き

    // 求解器の選択（hybrid method）
    gsl_multiroot_fsolver *solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids, 3);
    gsl_multiroot_fsolver_set(solver, &f, x);

    // 反復計算
    int status, iter = 0;
    do {
        iter++;
        status = gsl_multiroot_fsolver_iterate(solver);
        
        if (status) {
            printf("Solver failed\n");
            break;
        }

        printf("Iteration %d: x = %.7f, y = %.7f, z = %.7f\n", iter,
               gsl_vector_get(solver->x, 0),
               gsl_vector_get(solver->x, 1),
               gsl_vector_get(solver->x, 2));

    } while (gsl_vector_get(solver->f, 0) > 1e-7 || gsl_vector_get(solver->f, 1) > 1e-7 || gsl_vector_get(solver->f, 2) > 1e-7);

    // 結果の表示
    printf("Solution found: x = %.7f, y = %.7f, z = %.7f\n",
           gsl_vector_get(solver->x, 0),
           gsl_vector_get(solver->x, 1),
           gsl_vector_get(solver->x, 2));

    // メモリの解放
    gsl_multiroot_fsolver_free(solver);
    gsl_vector_free(x);

    return 0;
}

