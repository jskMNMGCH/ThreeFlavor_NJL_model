#include "three_flavor_njl.h"
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multiroots.h>


// Global constants (adjust as needed)
const double Nf = 3.0;     // Number of flavors
const double Nc = 3.0;     // Number of colors
const double HBARC = 0.17;
const double RHO0 = 197.3269804;
const double H= 1e-5;
const int WORKSPACE_SIZE = 5000;
const int NUM_SETS = 128;

double Lambda;  // Momentum cutoff (MeV)
double m[3];    // Tormat is {m_u, m_d, m_s}.
double G;       // Coupling const
double K;       // Coupling const

// Function definitions

// --------------------
// Fermi occupation functions
// --------------------
double np_f(double p, double T, double mu, double M) {
    double Ep = sqrt(p*p + M*M);
    return 1.0 / (exp((Ep - mu) / T) + 1.0);
}

double np_f_bar(double p, double T, double mu, double M) {
    double Ep = sqrt(p*p + M*M);
    return 1.0 / (exp((Ep + mu) / T) + 1.0);
}

// Computes the free thermodynamic potential (grand potential) for a Fermi gas.
/*
struct OmegaMf_params{
    double T;
    double mu_f;
    double Mf;
};
*/
double OmegaMf_integrand(double p, void* params_ptr){
    struct OmegaMf_params* params = (struct OmegaMf_params*) params_ptr;
    double T = params->T;
    double mu_f = params->mu_f;
    double Mf = params->Mf;
    double Ep = sqrt(p*p + Mf*Mf);
    double result;
    result = p*p *( Ep
            + T * log(1.0 + exp(-(Ep - mu_f) / T))
            + T * log(1.0 + exp(-(Ep + mu_f) / T))
            );
    return result;
}
double OmegaMf(double T, double mu_f, double Mf){
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(WORKSPACE_SIZE);
    struct OmegaMf_params params = {T, mu_f, Mf};
    double result, error;
    gsl_function F;
    F.function = &OmegaMf_integrand;
    F.params = &params;
    gsl_integration_qag(&F, 0.0, Lambda, 1e-8, 1e-8, WORKSPACE_SIZE, GSL_INTEG_GAUSS15, w, &result, &error);
    gsl_integration_workspace_free(w);
    return -Nc / (M_PI*M_PI) * result;
}
// Computes the total thermodynamic potential including interaction terms.
double Omega_temp(double T, double mu[3], double phi[3]){
    double Omega = 0.0;
    double M[3] = {
        m[0] - 4.0*G*phi[0] + 2.0*K*phi[1]*phi[2],
        m[1] - 4.0*G*phi[1] + 2.0*K*phi[2]*phi[0],
        m[2] - 4.0*G*phi[2] + 2.0*K*phi[0]*phi[1]
    };
    for (int i=0; i<3; ++i){
    Omega += OmegaMf(T, mu[i], M[i]) + 2.0*G*phi[i]*phi[i];
    }
    Omega -= 4.0*K*phi[0]*phi[1]*phi[2];
    return Omega + + 1.9397179193981308e10;
}
// Computes the derivative of Omega_temp with respect to mass M.
/*
struct phi_f_params{
    double T;
    double mu_f;
    double Mf;
};
*/
double phi_f_integrand(double p, void* params_ptr){
    struct phi_f_params* params = (struct phi_f_params*) params_ptr;
    double T = params->T;
    double mu_f = params->mu_f;
    double Mf = params->Mf;
    double Ep = sqrt(p*p + Mf*Mf);
    double result;
    result = p*p * Mf/Ep * ( 1.0
            - np_f(p, T, mu_f, Mf)
            - np_f_bar(p, T, mu_f, Mf)
            );
    return result;
}
double phi_f(double T, double mu_f, double Mf){
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(WORKSPACE_SIZE);
    struct phi_f_params params = {T, mu_f, Mf};
    double result, error;
    gsl_function F;
    F.function = &phi_f_integrand;
    F.params = &params;
    gsl_integration_qag(&F, 0.0, Lambda, 1e-8, 1e-8, WORKSPACE_SIZE, GSL_INTEG_GAUSS15, w, &result, &error);
    gsl_integration_workspace_free(w);
    return -Nc/ (M_PI*M_PI) *result;
}
/*
struct gap_eqs_params{
    double T;
    double mu_u;
    double mu_d;
    double mu_s;
};
*/
// 非線形連立方程式の定義
int gap_eqs(const gsl_vector* v, void* params, gsl_vector* f) {
    double M_u = gsl_vector_get(v, 0);
    double M_d = gsl_vector_get(v, 1);
    double M_s = gsl_vector_get(v, 2);

    struct gap_eqs_params* p = (struct gap_eqs_params*) params;
    double mu[3] = {p->mu_u, p->mu_d, p->mu_s};
    double T = p->T;
    double phi_u = phi_f(T, mu[0], M_u);
    double phi_d = phi_f(T, mu[1], M_d);
    double phi_s = phi_f(T, mu[2], M_s);
    gsl_vector_set(f, 0, m[0] -4.0 *G *phi_u +2.0 *K *phi_d *phi_s);
    gsl_vector_set(f, 1, m[1] -4.0 *G *phi_d +2.0 *K *phi_s *phi_u);
    gsl_vector_set(f, 2, m[2] -4.0 *G *phi_s +2.0 *K *phi_u *phi_d);
    return GSL_SUCCESS;
}

double solv_gap_eq(double T, double mu[3], double(* M)[3]){
    struct gap_eqs_params params = {
        .T = T,
        .mu_u = mu[0],
        .mu_d = mu[1],
        .mu_s = mu[2]
    };
    gsl_vector* M_vec = gsl_vector_alloc(3);
    gsl_vector_set(M_vec, 0, 4e2);
    gsl_vector_set(M_vec, 1, 4e2);
    gsl_vector_set(M_vec, 2, 6e2);
    // ソルバーの設定
    gsl_multiroot_function f = {&gap_eqs, 3, &params};
    // 求解器の選択（hybrid method)
    gsl_multiroot_fsolver *solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids, 3);
    if (!solver) {
        fprintf(stderr, "Failed to allocate solver\n");
        gsl_vector_free(M_vec);
        return -1;
    }
    int status = gsl_multiroot_fsolver_set(solver, &f, M_vec);
    if (status != GSL_SUCCESS) {
        fprintf(stderr, "Failed to set up solver: %s\n", gsl_strerror(status));
        gsl_multiroot_fsolver_free(solver);
        gsl_vector_free(M_vec);
        return -1;
    }
    // 反復計算
    int iter = 0;
    do {
        ++iter;
        status = gsl_multiroot_fsolver_iterate(solver);
        if (status) {
            printf("Solver failed at iteration %d: %s\n",iter,gsl_strerror(status));
            break;
        }
        // Ensure positive values
        for (int i = 0; i < 3; i++) {
            if (gsl_vector_get(solver->x, i) < 0.0) {
                gsl_vector_set(solver->x, i, 1e-8); // Adjust to a small positive value
            }
        }
        status = gsl_multiroot_test_residual(solver->f, 1e-8);
    } while (status == GSL_CONTINUE && iter < 5000);
    (*M)[0] = gsl_vector_get(solver->x, 0);
    (*M)[1] = gsl_vector_get(solver->x, 1);
    (*M)[2] = gsl_vector_get(solver->x, 2);
    gsl_multiroot_fsolver_free(solver);
    gsl_vector_free(M_vec);
    return status;
}

/*
struct phi_eqs_params {
    double T;
    double mu_u;
    double mu_d;
    double mu_s;
};
*/
// Nonlinear system of equations
int phi_eqs(const gsl_vector* v, void* params, gsl_vector* f) {
    double phi_u = gsl_vector_get(v, 0);
    double phi_d = gsl_vector_get(v, 1);
    double phi_s = gsl_vector_get(v, 2);

    struct phi_eqs_params *p = (struct phi_eqs_params*) params;
    double T = p->T;
    double mu[3] = {p->mu_u, p->mu_d, p->mu_s};
    double M_mf[3] = {
        m[0] - 4.0 * G * phi_u + 2.0 * K * phi_d * phi_s,
        m[1] - 4.0 * G * phi_d + 2.0 * K * phi_s * phi_u,
        m[2] - 4.0 * G * phi_s + 2.0 * K * phi_u * phi_d
    };
    gsl_vector_set(f, 0, phi_u - phi_f(T, mu[0], M_mf[0]));
    gsl_vector_set(f, 1, phi_d - phi_f(T, mu[1], M_mf[1]));
    gsl_vector_set(f, 2, phi_s - phi_f(T, mu[2], M_mf[2]));
    return GSL_SUCCESS;
}

double solv_phi_eq(double T, double mu[3], double(* phi)[3]) {
    struct phi_eqs_params params = {
        .T = T,
        .mu_u = mu[0],
        .mu_d = mu[1],
        .mu_s = mu[2]
    };
    gsl_vector* phi_vec = gsl_vector_alloc(3);
    gsl_vector_set(phi_vec, 0, 4e2);
    gsl_vector_set(phi_vec, 1, 4e2);
    gsl_vector_set(phi_vec, 2, 6e2);
    // Solver setup
    gsl_multiroot_function f = {&phi_eqs, 3, &params};
    gsl_multiroot_fsolver *solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids, 3);
    if (!solver) {
        fprintf(stderr, "Failed to allocate solver\n");
        gsl_vector_free(phi_vec);
        return -1;
    }
    gsl_multiroot_fsolver_set(solver, &f, phi_vec);
    // Iterative calculation
    int status, iter = 0;
    do {
        ++iter;
        status = gsl_multiroot_fsolver_iterate(solver);
        if (status) {
            printf("Solver failed at iteration %d: %s\n", iter, gsl_strerror(status));
            break;
        }
        // Ensure positive values
        for (int i = 0; i < 3; i++) {
            if (gsl_vector_get(solver->x, i) < 0.0) {
                gsl_vector_set(solver->x, i, 1e-8); // Adjust to a small positive value
            }
        }
        status = gsl_multiroot_test_residual(solver->f, 1e-8);
    } while (status == GSL_CONTINUE && iter < 1000);

    (*phi)[0] = gsl_vector_get(solver->x, 0);
    (*phi)[1] = gsl_vector_get(solver->x, 1);
    (*phi)[2] = gsl_vector_get(solver->x, 2);

    gsl_multiroot_fsolver_free(solver);
    gsl_vector_free(phi_vec);
    return status;
}

double Omega(double T, double mu[3]){
    double phi_sol[3];
    double mu_vac[3] = {0.0, 0.0, 0.0};
    double phi_vac[3] = {0.0, 0.0, 0.0};
    solv_phi_eq(T, mu, &phi_sol);
    return Omega_temp(T, mu, phi_sol) - Omega_temp(0.0, mu_vac, phi_vac);
}

double pressure(double T, double mu[3]){
    return -Omega(T, mu);
}

void calc_M(double phi_sol[3], double(* M_sol)[3]){
    (*M_sol)[0] = m[0] -4.0 *G *phi_sol[0] +2.0 *K *phi_sol[1] *phi_sol[2];
    (*M_sol)[1] = m[1] -4.0 *G *phi_sol[1] +2.0 *K *phi_sol[2] *phi_sol[0];
    (*M_sol)[2] = m[2] -4.0 *G *phi_sol[2] +2.0 *K *phi_sol[0] *phi_sol[1];
}


