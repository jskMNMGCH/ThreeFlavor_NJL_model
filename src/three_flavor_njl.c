#include "three_flavor_njl.h"
#include <stdio.h>
#include <math.h>
#include <string.h>  // memcpyを使うために必要
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multiroots.h>

// Global constants (adjust as needed)
const double Nf = 3.0;     // Number of flavors
const double Nc = 3.0;     // Number of colors
const double HBARC = 0.17;
const double RHO0 = 197.3269804;
const double H = 1e-4;
const int WORKSPACE_SIZE = 5000;
const int NUM_STEPS = 8000;

double Lambda;  // Momentum cutoff (MeV)
double m[3];    // Tormat is {m_u, m_d, m_s}.
double G;       // Coupling constant
double K;       // Coupling constant

// Function definitions

double np_f(double p, double T, double mu, double M) {
//  --------------------
//  Computes the Fermi-Dirac distribution function for particles
//  --------------------
    double Ep = sqrt(p*p + M*M);
    return 1.0 / (exp((Ep - mu) / T) + 1.0);
}

double np_f_bar(double p, double T, double mu, double M) {
//  --------------------
//  Computes the Fermi-Dirac distribution function for antiparticles
//  --------------------
    double Ep = sqrt(p*p + M*M);
    return 1.0 / (exp((Ep + mu) / T) + 1.0);
}

double OmegaMf_integrand(double p, void* params_ptr) {
//  --------------------
//  Computes the integrand of the free thermodynamic potential for a Fermi gas
//  --------------------
    struct OmegaMf_params* params = (struct OmegaMf_params*) params_ptr;
    double T = params->T;
    double mu_f = params->mu_f;
    double Mf = params->Mf;
    double Ep = sqrt(p*p + Mf*Mf);
    double result;
    result = p*p * (Ep
            + T * log(1.0 + exp(-(Ep - mu_f) / T))
            + T * log(1.0 + exp(-(Ep + mu_f) / T))
            );
    return result;
}

double OmegaMf(double T, double mu_f, double Mf) {
//  --------------------
//  Computes the free thermodynamic potential (grand potential) for a Fermi gas
//  --------------------
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(WORKSPACE_SIZE);
    struct OmegaMf_params params = {T, mu_f, Mf};
    double result, error;
    gsl_function F;
    F.function = &OmegaMf_integrand;
    F.params = &params;
    gsl_integration_qag(&F, 0.0, Lambda, H, H, WORKSPACE_SIZE, GSL_INTEG_GAUSS15, w, &result, &error);
    gsl_integration_workspace_free(w);
    return -Nc / (M_PI*M_PI) * result;
}

double Omega_temp(double T, double mu[3], double phi[3]) {
//  --------------------
//  Computes the total thermodynamic potential including interaction terms
//  --------------------
    double Omega = 0.0;
    double M[3] = {
        m[0] - 4.0*G*phi[0] + 2.0*K*phi[1]*phi[2],
        m[1] - 4.0*G*phi[1] + 2.0*K*phi[2]*phi[0],
        m[2] - 4.0*G*phi[2] + 2.0*K*phi[0]*phi[1]
    };
    for (int i = 0; i < 3; ++i) {
        Omega += OmegaMf(T, mu[i], M[i]) + 2.0*G*phi[i]*phi[i];
    }
    Omega -= 4.0*K*phi[0]*phi[1]*phi[2];
    return Omega + 1.9397179193981308e10;
}

double phi_f_integrand(double p, void* params_ptr) {
//  --------------------
//  Computes the integrand of the derivative of Omega_temp with respect to momentum `p`
//  --------------------
    struct phi_f_params* params = (struct phi_f_params*) params_ptr;
    double T = params->T;
    double mu_f = params->mu_f;
    double Mf = params->Mf;
    double Ep = sqrt(p*p + Mf*Mf);
    double result;
    result = p*p * Mf/Ep * (1.0
            - np_f(p, T, mu_f, Mf)
            - np_f_bar(p, T, mu_f, Mf)
            );
    return result;
}

double phi_f(double T, double mu_f, double Mf) {
//  --------------------
//  Computes the derivative of Omega_temp with respect to mass M
//  --------------------
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(WORKSPACE_SIZE);
    struct phi_f_params params = {T, mu_f, Mf};
    double result, error;
    gsl_function F;
    F.function = &phi_f_integrand;
    F.params = &params;
    gsl_integration_qag(&F, 0.0, Lambda, H, H, WORKSPACE_SIZE, GSL_INTEG_GAUSS15, w, &result, &error);
    gsl_integration_workspace_free(w);
    return -Nc / (M_PI*M_PI) * result;
}

int gap_eqs(const gsl_vector* v, void* params, gsl_vector* f) {
//  --------------------
//  Defines the nonlinear system of equations for the gap equations
//  --------------------
    double M_u = gsl_vector_get(v, 0);
    double M_d = gsl_vector_get(v, 1);
    double M_s = gsl_vector_get(v, 2);

    struct gap_eqs_params* p = (struct gap_eqs_params*) params;
    double mu[3] = {p->mu_u, p->mu_d, p->mu_s};
    double T = p->T;
    double phi_u = phi_f(T, mu[0], M_u);
    double phi_d = phi_f(T, mu[1], M_d);
    double phi_s = phi_f(T, mu[2], M_s);
    gsl_vector_set(f, 0, m[0] - 4.0*G*phi_u + 2.0*K*phi_d*phi_s -M_u);
    gsl_vector_set(f, 1, m[1] - 4.0*G*phi_d + 2.0*K*phi_s*phi_u -M_d);
    gsl_vector_set(f, 2, m[2] - 4.0*G*phi_s + 2.0*K*phi_u*phi_d -M_s);
    return GSL_SUCCESS;
}

double solv_gap_eq(double T, double mu[3], double(* M)[3]) {
//  --------------------
//  Solves the gap equations using a nonlinear solver
//  --------------------
    struct gap_eqs_params params = {
        .T = T,
        .mu_u = mu[0],
        .mu_d = mu[1],
        .mu_s = mu[2]
    };
    gsl_vector* M_vec = gsl_vector_alloc(3);
    gsl_vector_set(M_vec, 0, 3.7e2);
    gsl_vector_set(M_vec, 1, 4e2);
    gsl_vector_set(M_vec, 2, 5.5e2);
    // Solver setup
    gsl_multiroot_function f = {&gap_eqs, 3, &params};
    // Select solver (hybrid method)
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
    // Iterative calculation
    int iter = 0;
    do {
        ++iter;
        status = gsl_multiroot_fsolver_iterate(solver);
        if (status) {
            printf("Solver failed at iteration %d: %s\n", iter, gsl_strerror(status));
            break;
        }
        // Ensure positive values
        for (int i = 0; i < 3; ++i) {
            if (gsl_vector_get(solver->x, i) < 0.0) {
                gsl_vector_set(solver->x, i, H); // Adjust to a small positive value
            }
        }
        status = gsl_multiroot_test_residual(solver->f, H);
    } while (status == GSL_CONTINUE && iter < NUM_STEPS);
    (*M)[0] = gsl_vector_get(solver->x, 0);
    (*M)[1] = gsl_vector_get(solver->x, 1);
    (*M)[2] = gsl_vector_get(solver->x, 2);
    gsl_multiroot_fsolver_free(solver);
    gsl_vector_free(M_vec);
    return status;
}

int phi_eqs(const gsl_vector* v, void* params, gsl_vector* f) {
//  --------------------
//  Defines the nonlinear system of equations for the phi equations
//  --------------------
    double phi_u = gsl_vector_get(v, 0);
    double phi_d = gsl_vector_get(v, 1);
    double phi_s = gsl_vector_get(v, 2);

    struct phi_eqs_params* p = (struct phi_eqs_params*) params;
    double T = p->T;
    double mu[3] = {p->mu_u, p->mu_d, p->mu_s};
    double M_mf[3] = {
        m[0] - 4.0*G*phi_u + 2.0*K*phi_d*phi_s,
        m[1] - 4.0*G*phi_d + 2.0*K*phi_s*phi_u,
        m[2] - 4.0*G*phi_s + 2.0*K*phi_u*phi_d
    };
    gsl_vector_set(f, 0, phi_u - phi_f(T, mu[0], M_mf[0]));
    gsl_vector_set(f, 1, phi_d - phi_f(T, mu[1], M_mf[1]));
    gsl_vector_set(f, 2, phi_s - phi_f(T, mu[2], M_mf[2]));
    return GSL_SUCCESS;
}

double solv_phi_eq(double T, double mu[3], double(* phi)[3]) {
//  --------------------
//  Solves the phi equations using a nonlinear solver
//  --------------------
    struct phi_eqs_params params = {
        .T = T,
        .mu_u = mu[0],
        .mu_d = mu[1],
        .mu_s = mu[2]
    };
    gsl_vector* phi_vec = gsl_vector_alloc(3);
    gsl_vector_set(phi_vec, 0, 8e7);
    gsl_vector_set(phi_vec, 1, 8e7);
    gsl_vector_set(phi_vec, 2, 8e7);
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
        status = gsl_multiroot_test_residual(solver->f, H);
    } while (status == GSL_CONTINUE && iter < NUM_STEPS);

    (*phi)[0] = gsl_vector_get(solver->x, 0);
    (*phi)[1] = gsl_vector_get(solver->x, 1);
    (*phi)[2] = gsl_vector_get(solver->x, 2);

    gsl_multiroot_fsolver_free(solver);
    gsl_vector_free(phi_vec);
    return status;
}

double Omega(double T, double mu[3]) {
//  --------------------
//  Computes the thermodynamic potential Omega at a given temperature and chemical potentials
//  --------------------
    double phi_sol[3];
    double mu_vac[3] = {0.0, 0.0, 0.0};
    double phi_vac[3] = {0.0, 0.0, 0.0};
    solv_phi_eq(T, mu, &phi_sol);
    return Omega_temp(T, mu, phi_sol) - Omega_temp(0.0, mu_vac, phi_vac);
}

double pressure(double T, double mu[3]) {
//  --------------------
//  Computes the pressure at a given temperature and chemical potentials
//  --------------------
    return -Omega(T, mu);
}

int calc_M(double phi_sol[3], double(* M_sol)[3]) {
//  --------------------
//  Calculates the constituent quark masses given the solution for the chiral condensates
//  --------------------
    (*M_sol)[0] = m[0] - 4.0*G*phi_sol[0] + 2.0*K*phi_sol[1]*phi_sol[2];
    (*M_sol)[1] = m[1] - 4.0*G*phi_sol[1] + 2.0*K*phi_sol[2]*phi_sol[0];
    (*M_sol)[2] = m[2] - 4.0*G*phi_sol[2] + 2.0*K*phi_sol[0]*phi_sol[1];
    return 0;
}

int M_phi_eqs(const gsl_vector* v, void* params, gsl_vector* f){
    gsl_vector* phi_vec = gsl_vector_alloc(3);
    double phi_u = gsl_vector_get(v, 0);
    double phi_d = gsl_vector_get(v, 1);
    double phi_s = gsl_vector_get(v, 2);

    struct M_phi_eqs_params* p = (struct M_phi_eqs_params*) params;
    double M_u = p->M_u;
    double M_d = p->M_d;
    double M_s = p->M_s;
    // Solver setup
    gsl_vector_set(f, 0, m[0] - 4.0*G*phi_u + 2.0*K*phi_d*phi_s -M_u);
    gsl_vector_set(f, 1, m[1] - 4.0*G*phi_d + 2.0*K*phi_s*phi_u -M_d);
    gsl_vector_set(f, 2, m[2] - 4.0*G*phi_s + 2.0*K*phi_u*phi_d -M_s);
    return GSL_SUCCESS;
}

int calc_phi(double M_sol[3], double(* phi_sol)[3]) {
//  --------------------
//  Calculates the quark condensates given the solution for the constituent quark masses
//  --------------------
    struct M_phi_eqs_params params = {
        .M_u = M_sol[0],
        .M_d = M_sol[1],
        .M_s = M_sol[2]
    };
    // Solver setup
    gsl_vector* phi_vec = gsl_vector_alloc(3);
    gsl_vector_set(phi_vec, 0, 7.226774675685889e7);
    gsl_vector_set(phi_vec, 1, 7.226774675685889e7);
    gsl_vector_set(phi_vec, 2, 7.226774675685889e7);
    gsl_multiroot_function f = {&M_phi_eqs, 3, &params};
    // Select solver (hybrid method)
    gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids, 3);
    if (!solver) {
        fprintf(stderr, "Failed to allcocate solver\n");
       gsl_vector_free(phi_vec);
      return -1; 
    }
    int status = gsl_multiroot_fsolver_set(solver, &f, phi_vec);
    if (status != GSL_SUCCESS) {
        fprintf(stderr, "Failed to set up solver: %s\n", gsl_strerror(status));
        gsl_multiroot_fsolver_free(solver);
        gsl_vector_free(phi_vec);
        return -1;
    }
    // Iterative calculation
    int iter = 0;
    do {
        ++iter;
        status = gsl_multiroot_fsolver_iterate(solver);
        if (status) {
            printf("Solver failed at iteration %d: %s\n",iter,gsl_strerror(status));
            break;
        }  
        // Ensure positive values
        for (int i = 0; i < 3; ++i) {
            if (gsl_vector_get(solver->x, i) < 0.0) {
                gsl_vector_set(solver->x, i, H); // Adjust to a small positive value
            }
        }
        status = gsl_multiroot_test_residual(solver->f, H);
    } while (status == GSL_CONTINUE && iter < NUM_STEPS);
    (*phi_sol)[0] = gsl_vector_get(solver->x, 0);
    (*phi_sol)[1] = gsl_vector_get(solver->x, 1);
    (*phi_sol)[2] = gsl_vector_get(solver->x, 2);
    gsl_multiroot_fsolver_free(solver);
    gsl_vector_free(phi_vec);
    return status;
}


double solv_gap_eq_multi_guess(double T, double mu[3], double(* M)[3]) {
//  --------------------
//  Solves the gap equations using a nonlinear solver
//  --------------------
    struct gap_eqs_params params = {
        .T = T,
        .mu_u = mu[0],
        .mu_d = mu[1],
        .mu_s = mu[2]
    };
    double M_u_guess[3] = {367.6, 52.5, m[0]};
    double M_d_guess[3] = {367.6, 52.5, m[1]};
    double M_s_guess[3] = {549.5, 464.4, m[2]};
    
    int guess_best = 100;
    double M_best[3] = {M_u_guess[0], M_d_guess[0], M_s_guess[0]};
    double Omega_best = HUGE_VAL;
    double M_tmp[3];
    double phi_tmp[3];
    double Omega_tmp;
    
    // Solver setup
    int status;
    gsl_multiroot_function f = {&gap_eqs, 3, &params};
    // Select solver (hybrid method)
    gsl_multiroot_fsolver *solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_dnewton, 3);
    if (!solver) {
        fprintf(stderr, "Failed to allocate solver at gap_eq_multi_guess\n");
        return -1;
    }
    gsl_vector* M_vec = gsl_vector_alloc(3);
    for (int i=0; i < sizeof(M_u_guess)/sizeof(M_u_guess[0]); ++i) {
        gsl_vector_set(M_vec, 0, M_u_guess[i]);
        gsl_vector_set(M_vec, 1, M_d_guess[i]);
        gsl_vector_set(M_vec, 2, M_s_guess[i]);
        status = gsl_multiroot_fsolver_set(solver, &f, M_vec);
        if (status != GSL_SUCCESS) {
            fprintf(stderr, "Failed to set up solver at gap_eq_multi_guess i=%d: %s\n",i,gsl_strerror(status));
            gsl_multiroot_fsolver_free(solver);
            gsl_vector_free(M_vec);
            return -1;
        }
        // Iterative calculation
        int iter = 0;
        do {
            ++iter;
            status = gsl_multiroot_fsolver_iterate(solver);
            if (status) {
                //printf("Stop when params. are T=%lf, mu[0]=%lf, i=%d\n", T,mu[0],i);
                break;
            }
            status = gsl_multiroot_test_delta(solver->dx, solver->x, H, H);
        } while (status == GSL_CONTINUE && iter < NUM_STEPS);
        M_tmp[0] = gsl_vector_get(solver->x, 0);
        M_tmp[1] = gsl_vector_get(solver->x, 1);
        M_tmp[2] = gsl_vector_get(solver->x, 2);
        calc_phi(M_tmp, &phi_tmp);
        // calc_phi_multi_guess(T, mu, M_tmp, &phi_tmp);
        // Check consistency.
        /*
        double M_check[3];
        calc_M(phi_tmp, &M_check);
        if (fabs(M_tmp[2] - M_check[2]) > 1.0) {
            printf("Phi is inconsistent: Ms_tmp=%lf, Ms_check=%lf\n", M_tmp[2], M_check[2]);
        }
        */
        // Check the Omega
        Omega_tmp = Omega_temp(T, mu, phi_tmp);
        //printf("Compare Omega: Omega_tmp=%lf, phi_tmp[0]=%lf, for M_u=%lf at mu[0]=%lf\n",Omega_tmp,phi_tmp[0],M_tmp[0],mu[0]);
        //printf("Strange Mass at mu[2]=%lf, id=%d: %lf\n", mu[2], i, M_tmp[2]);
        if (Omega_tmp <= Omega_best) {
            Omega_best = Omega_tmp;
            memcpy(M_best, M_tmp, sizeof(M_tmp));
            guess_best = i;
        }
    }
    gsl_vector_free(M_vec);
    //printf("Solve end T=%lf, mu[0]=%lf, guess_besl=%d\n", T,mu[0],guess_best);
    memcpy((*M), M_best, sizeof(M_best));
    gsl_multiroot_fsolver_free(solver);
    return status;
}

int calc_phi_multi_guess(double T, double mu[3], double M_sol[3], double(* phi_sol)[3]) {
//  --------------------
//  Calculates the quark condensates given the solution for the constituent quark masses
//  --------------------
    struct M_phi_eqs_params params = {
        .M_u = M_sol[0],
        .M_d = M_sol[1],
        .M_s = M_sol[2]
    };
    // Select solver (hybrid method)
    gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids, 3);
    if (!solver) {
        fprintf(stderr, "Failed to allcocate solver\n");
        return -1; 
    }
    int status;
    double phi_u_guess[2] = {8e7, 6.5e7};
    double phi_d_guess[2] = {8e7, 6.5e7};
    double phi_s_guess[2] = {8e7, 6.4e7};
    double phi_best[3] = {phi_u_guess[0], phi_d_guess[0], phi_s_guess[0]};
    double Omega_best = HUGE_VAL;
    double phi_tmp[3];
    double Omega_tmp;
    // Solver setup
    gsl_vector* phi_vec = gsl_vector_alloc(3);

    for (int i = 0; i < 2; ++i) {
    gsl_vector_set(phi_vec, 0, phi_u_guess[i]);
    gsl_vector_set(phi_vec, 1, phi_d_guess[i]);
    gsl_vector_set(phi_vec, 2, phi_s_guess[i]);
    gsl_multiroot_function f = {&M_phi_eqs, 3, &params};

    status = gsl_multiroot_fsolver_set(solver, &f, phi_vec);
    if (status != GSL_SUCCESS) {
        fprintf(stderr, "Failed to set up solver at calc_phi_multi_guess: %s\n", gsl_strerror(status));
        gsl_multiroot_fsolver_free(solver);
        gsl_vector_free(phi_vec);
        return -1;
    }
    // Iterative calculation
    int iter = 0;
    do {
        ++iter;
        status = gsl_multiroot_fsolver_iterate(solver);
        if (status) {
            break;
        }  
        status = gsl_multiroot_test_delta(solver->dx, solver->x, H, H);
    } while (status == GSL_CONTINUE && iter < NUM_STEPS);
    phi_tmp[0] = gsl_vector_get(solver->x, 0);
    phi_tmp[1] = gsl_vector_get(solver->x, 1);
    phi_tmp[2] = gsl_vector_get(solver->x, 2);
    Omega_tmp = Omega_temp(T, mu, phi_tmp);
    if (Omega_tmp < Omega_best) {
        Omega_best = Omega_tmp;
        memcpy(phi_best, phi_tmp, sizeof(phi_tmp));
    }
    }
    memcpy((*phi_sol), phi_best, sizeof(phi_best));
    gsl_multiroot_fsolver_free(solver);
    gsl_vector_free(phi_vec);
    return status;
}

double solv_gap_eq_multi_guess_solver(double T, double mu[3], double(* M)[3]) {
    //  --------------------
    //  Solves the gap equations using a nonlinear solver
    //  --------------------
    struct gap_eqs_params params = {
        .T = T,
        .mu_u = mu[0],
        .mu_d = mu[1],
        .mu_s = mu[2]
    };
    double M_u_guess[] = {367.6, 52.5,  m[0]};
    double M_d_guess[] = {367.6, 52.5,  m[1]};
    double M_s_guess[] = {549.5, 464.4, m[2]};
    
    int guess_best = 100;
    double M_best[3] = {M_u_guess[0], M_d_guess[0], M_s_guess[0]};
    double Omega_best = HUGE_VAL;
    double M_tmp[3];
    double phi_tmp[3];
    double Omega_tmp;
    
    // Solver setup
    int status;
    gsl_multiroot_function f = {&gap_eqs, 3, &params};
    // Select solver (hybrid method)
    gsl_multiroot_fsolver *solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_dnewton, 3);
    if (!solver) {
        fprintf(stderr, "Failed to allocate solver at gap_eq_multi_guess\n");
        return -1;
    }
    gsl_vector* M_vec = gsl_vector_alloc(3);
    for (int i=0; i < sizeof(M_u_guess)/sizeof(M_u_guess[0]); ++i) {
        gsl_vector_set(M_vec, 0, M_u_guess[i]);
        gsl_vector_set(M_vec, 1, M_d_guess[i]);
        gsl_vector_set(M_vec, 2, M_s_guess[i]);
        status = gsl_multiroot_fsolver_set(solver, &f, M_vec);
        if (status != GSL_SUCCESS) {
            fprintf(stderr, "Failed to set up solver at gap_eq_multi_guess i=%d: %s\n",i,gsl_strerror(status));
            gsl_multiroot_fsolver_free(solver);
            gsl_vector_free(M_vec);
            return -1;
        }
        // Iterative calculation
        int iter = 0;
        do {
            ++iter;
            status = gsl_multiroot_fsolver_iterate(solver);
            if (status) {
                //printf("Stop when params. are T=%lf, mu[0]=%lf, i=%d\n", T,mu[0],i);
                break;
            }
            //status = gsl_multiroot_test_delta(solver->dx, solver->x, sqrt(H), sqrt(H));
            status = gsl_multiroot_test_residual(solver->f, H);
        } while (status == GSL_CONTINUE && iter < NUM_STEPS);
        
        // If iteration count reaches maximum, switch to gsl_multiroot_fsolver_hybrids
        if (iter >= NUM_STEPS) {
            gsl_vector_set(M_vec, 0, M_u_guess[i]);
            gsl_vector_set(M_vec, 1, M_d_guess[i]);
            gsl_vector_set(M_vec, 2, M_s_guess[i]);
            gsl_multiroot_fsolver_free(solver);
            solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids, 3);
            if (!solver) {
                fprintf(stderr, "Failed to allocate gsl_multiroot_fsolver_hybrids solver at gap_eq_multi_guess\n");
                gsl_vector_free(M_vec);
                return -1;
            }
            status = gsl_multiroot_fsolver_set(solver, &f, M_vec);
            if (status != GSL_SUCCESS) {
                fprintf(stderr, "Failed to set up gsl_multiroot_fsolver_hybrids solver at gap_eq_multi_guess i=%d: %s\n",i,gsl_strerror(status));
                gsl_multiroot_fsolver_free(solver);
                gsl_vector_free(M_vec);
                return -1;
            }
            iter = 0;
            do {
                ++iter;
                status = gsl_multiroot_fsolver_iterate(solver);
                if (status) {
                    break;
                }
                status = gsl_multiroot_test_delta(solver->dx, solver->x, H, H);
            } while (status == GSL_CONTINUE && iter < NUM_STEPS);
        }

        M_tmp[0] = gsl_vector_get(solver->x, 0);
        M_tmp[1] = gsl_vector_get(solver->x, 1);
        M_tmp[2] = gsl_vector_get(solver->x, 2);
        calc_phi(M_tmp, &phi_tmp);
        // Check the Omega
        Omega_tmp = Omega_temp(T, mu, phi_tmp);
        if (Omega_tmp <= Omega_best) {
            Omega_best = Omega_tmp;
            memcpy(M_best, M_tmp, sizeof(M_tmp));
            guess_best = i;
        }
    }
    gsl_vector_free(M_vec);
    //printf("Solve end T=%lf, mu[0]=%lf, guess_besl=%d\n", T,mu[0],guess_best);
    memcpy((*M), M_best, sizeof(M_best));
    gsl_multiroot_fsolver_free(solver);
    return status;
}



double solv_gap_eq_multi_guess_solver_restrict(double T, double mu[3], double(* M)[3]) {
    //  --------------------
    //  Solves the gap equations using a nonlinear solver
    //  --------------------
    struct gap_eqs_params params = {
        .T = T,
        .mu_u = mu[0],
        .mu_d = mu[1],
        .mu_s = mu[2]
    };

    // Upper bounds when mu={0.0, 0.0, 0.0}
    double upper_bound[3];
    double mu_vac[3] = {0.0, 0.0, 0.0};
    solv_gap_eq_multi_guess_solver(T, mu_vac, &upper_bound);
    // Middle masses when mu={3.8e2, 3.8e2, 3.8e2}
    double middle_mass[3];
    double mu_middle[3] = {4e2, 4e2, 4e2};
    solv_gap_eq_multi_guess_solver(T, mu_middle, &middle_mass);
    
    double M_u_guess[] = {upper_bound[0], middle_mass[0], 2.0*m[0]};
    double M_d_guess[] = {upper_bound[1], middle_mass[1], 2.0*m[1]};
    double M_s_guess[] = {upper_bound[2], middle_mass[2], 2.0*m[2]};

    int guess_best = 100;
    double M_best[3] = {M_u_guess[0], M_d_guess[0], M_s_guess[0]};
    double Omega_best = HUGE_VAL;
    double M_tmp[3];
    double phi_tmp[3];
    double Omega_tmp;

    // Solver setup
    int status;
    gsl_multiroot_function f = {&gap_eqs, 3, &params};
    // Select solver (hybrid method)
    gsl_multiroot_fsolver *solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_dnewton, 3);
    if (!solver) {
        fprintf(stderr, "Failed to allocate solver at gap_eq_multi_guess\n");
        return -1;
    }
    gsl_vector* M_vec = gsl_vector_alloc(3);
    for (int i=0; i < sizeof(M_u_guess)/sizeof(M_u_guess[0]); ++i) {
        gsl_vector_set(M_vec, 0, M_u_guess[i]);
        gsl_vector_set(M_vec, 1, M_d_guess[i]);
        gsl_vector_set(M_vec, 2, M_s_guess[i]);
        status = gsl_multiroot_fsolver_set(solver, &f, M_vec);
        if (status != GSL_SUCCESS) {
            fprintf(stderr, "Failed to set up solver at gap_eq_multi_guess i=%d: %s\n",i,gsl_strerror(status));
            gsl_multiroot_fsolver_free(solver);
            gsl_vector_free(M_vec);
            return -1;
        }
        // Iterative calculation
        int iter = 0;
        do {
            ++iter;
            status = gsl_multiroot_fsolver_iterate(solver);
            if (status) {
                //printf("Stop when params. are T=%lf, mu[0]=%lf, i=%d\n", T,mu[0],i);
                break;
            }
            // Ensure values are within bounds
            /*
            for (int j = 0; j < 3; ++j) {
                if (gsl_vector_get(solver->x, j) < 0.0) {
                    gsl_vector_set(solver->x, j, H); // Adjust to a small positive value
                }
                else if (gsl_vector_get(solver->x, j) > upper_bound[j]+sqrt(H)) {
                    gsl_vector_set(solver->x, j, upper_bound[j]); // Adjust to the upper bound
                }
            }
            */
            status = gsl_multiroot_test_residual(solver->f, H);
        } while (status == GSL_CONTINUE && iter < NUM_STEPS);

        // If iteration count reaches maximum, switch to gsl_multiroot_fsolver_hybrids
        if (iter >= NUM_STEPS) {
            gsl_vector_set(M_vec, 0, M_u_guess[i]);
            gsl_vector_set(M_vec, 1, M_d_guess[i]);
            gsl_vector_set(M_vec, 2, M_s_guess[i]);
            gsl_multiroot_fsolver_free(solver);
            solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids, 3);
            if (!solver) {
                fprintf(stderr, "Failed to allocate gsl_multiroot_fsolver_hybrids solver at gap_eq_multi_guess\n");
                gsl_vector_free(M_vec);
                return -1;
            }
            status = gsl_multiroot_fsolver_set(solver, &f, M_vec);
            if (status != GSL_SUCCESS) {
                fprintf(stderr, "Failed to set up gsl_multiroot_fsolver_hybrids solver at gap_eq_multi_guess i=%d: %s\n",i,gsl_strerror(status));
                gsl_multiroot_fsolver_free(solver);
                gsl_vector_free(M_vec);
                return -1;
            }
            iter = 0;
            do {
                ++iter;
                status = gsl_multiroot_fsolver_iterate(solver);
                if (status) {
                    break;
                }
                // Ensure values are within bounds
                if (gsl_vector_get(solver->x, 0) < 0.0
                    || gsl_vector_get(solver->x, 1) < 0.0
                    || gsl_vector_get(solver->x, 2) < 0.0) {
                    gsl_vector_set(solver->x, 0, m[0]);
                    gsl_vector_set(solver->x, 1, m[1]);
                    gsl_vector_set(solver->x, 2, m[2]);
                }else if (gsl_vector_get(solver->x, 0) > upper_bound[0]+H
                          || gsl_vector_get(solver->x, 1) > upper_bound[1]+H
                          || gsl_vector_get(solver->x, 2) > upper_bound[2]+H) {
                    gsl_vector_set(solver->x, 0, middle_mass[0]);
                    gsl_vector_set(solver->x, 1, middle_mass[1]);
                    gsl_vector_set(solver->x, 2, middle_mass[2]);
                }
                status = gsl_multiroot_test_delta(solver->dx, solver->x, H, H);
            } while (status == GSL_CONTINUE && iter < NUM_STEPS);
        }

        M_tmp[0] = gsl_vector_get(solver->x, 0);
        M_tmp[1] = gsl_vector_get(solver->x, 1);
        M_tmp[2] = gsl_vector_get(solver->x, 2);
        calc_phi(M_tmp, &phi_tmp);
        // Check the Omega
        Omega_tmp = Omega_temp(T, mu, phi_tmp);
        if (Omega_tmp <= Omega_best) {
            Omega_best = Omega_tmp;
            memcpy(M_best, M_tmp, sizeof(M_tmp));
            guess_best = i;
        }
    }
    gsl_vector_free(M_vec);
    //printf("Solve end T=%lf, mu[0]=%lf, guess_besl=%d\n", T,mu[0],guess_best);
    memcpy((*M), M_best, sizeof(M_best));
    gsl_multiroot_fsolver_free(solver);
    return status;
}
