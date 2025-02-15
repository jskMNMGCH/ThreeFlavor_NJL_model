#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

// Define global constants
const double Nf = 2.0;
const double Nc = 3.0;
// For coupling constants and current mass:
double Lambda;
double m;  // will be given in main.
double Gs; // will be computed as 2.44 / (Lambda^2)
double Gv; // set equal to Gs

// Helper: GSL integration workspace size
#define HBARC 197.3269804  // hbar * c in MeV * fm
#define RHO0 0.17  // Nuclear saturation density
#define WORKSPACE_SIZE 1000
#define NUM_SETS 128
#define H 1e-5

// ------------------
// Fermi occupation functions
// ------------------
double np(double p, double T, double mu, double M) {
    double Ep = sqrt(p*p + M*M);
    return 1.0 / (exp((Ep - mu) / T) + 1.0);
}

double np_bar(double p, double T, double mu, double M) {
    double Ep = sqrt(p*p + M*M);
    return 1.0 / (exp((Ep + mu) / T) + 1.0);
}

// ------------------
// OmegaM: free thermodynamic potential
// ------------------
struct OmegaM_params {
    double Temp;
    double mu_t;
    double Mass;
};

double OmegaM_integrand(double p, void* params_ptr) {
    struct OmegaM_params *params = (struct OmegaM_params *) params_ptr;
    double Temp = params->Temp;
    double mu_t = params->mu_t;
    double Mass = params->Mass;
    double Ep = sqrt(p*p + Mass*Mass);
    double result;
    if (Temp > 0.0) {
        result = p*p * (Ep 
            + Temp * log(1.0 + exp(-(Ep - mu_t) / Temp))
            + Temp * log(1.0 + exp(-(Ep + mu_t) / Temp)));
    } else {
        double step = (mu_t - Ep >= 0.0) ? 1.0 : 0.0;
        result = p*p * (Ep + (mu_t - Ep)*step);
    }
    return result;
}

double OmegaM(double Temp, double mu_t, double Mass) {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(WORKSPACE_SIZE);
    struct OmegaM_params params = {Temp, mu_t, Mass};
    double result, error;
    gsl_function F;
    F.function = &OmegaM_integrand;
    F.params = &params;
    gsl_integration_qag(&F, 0.0, Lambda, 1e-8, 1e-8, WORKSPACE_SIZE, GSL_INTEG_GAUSS15, w, &result, &error);
    gsl_integration_workspace_free(w);
    return - Nf * Nc / (M_PI * M_PI) * result;
}

// ------------------
// Omega_temp: total thermodynamic potential
// ------------------
double Omega_temp(double T, double mu, double M, double mu_tilde) {
    return OmegaM(T, mu_tilde, M) + (M - m)*(M - m) / (4.0 * Gs) 
           - (mu - mu_tilde)*(mu - mu_tilde) / (4.0 * Gv)
           + 1.9397179193981308e10;
}

// ------------------
// dOmegadM: derivative with respect to Mass
// ------------------
struct dOmegaM_params {
    double T;
    double mu;
    double Mass;
    double mu_tilde;
};

double dOmegadM_integrand(double p, void *params_ptr) {
    struct dOmegaM_params *params = (struct dOmegaM_params *) params_ptr;
    double T = params->T;
    double mu = params->mu;
    double Mass = params->Mass;
    double mu_tilde = params->mu_tilde;
    double Ep = sqrt(p*p + Mass*Mass);
    double result;
    if (T > 0.0) {
        result = p*p * Mass / Ep * (1.0 - np(p, T, mu_tilde, Mass) - np_bar(p, T, mu_tilde, Mass));
    } else {
        double step = (mu_tilde - Ep >= 0.0) ? 1.0 : 0.0;
        result = p*p * Mass / Ep * (1.0 - step);
    }
    return result;
}

double dOmegadM(double T, double mu, double Mass, double mu_tilde) {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(WORKSPACE_SIZE);
    struct dOmegaM_params params = {T, mu, Mass, mu_tilde};
    double result, error;
    gsl_function F;
    F.function = &dOmegadM_integrand;
    F.params = &params;
    gsl_integration_qag(&F, 0.0, Lambda, 1e-8, 1e-8, WORKSPACE_SIZE, GSL_INTEG_GAUSS15, w, &result, &error);
    gsl_integration_workspace_free(w);
    return (Mass - m) / (2.0 * Gs) - Nf * Nc / (M_PI * M_PI) * result;
}

// ------------------
// dOmegadmutilde: derivative with respect to mu_tilde
// ------------------
double dOmegadmutilde_integrand(double p, void *params_ptr) {
    struct dOmegaM_params *params = (struct dOmegaM_params *) params_ptr;
    double T = params->T;
    double mu = params->mu;
    double Mass = params->Mass;
    double mu_tilde = params->mu_tilde;
    double Ep = sqrt(p*p + Mass*Mass);
    double result;
    if (T > 0.0) {
        result = p*p * (np(p, T, mu_tilde, Mass) - np_bar(p, T, mu_tilde, Mass));
    } else {
        result = p*p * ((mu_tilde - Ep) >= 0.0 ? 1.0 : 0.0);
    }
    return result;
}

double dOmegadmutilde(double T, double mu, double Mass, double mu_tilde) {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(WORKSPACE_SIZE);
    struct dOmegaM_params params = {T, mu, Mass, mu_tilde};
    double result, error;
    gsl_function F;
    F.function = &dOmegadmutilde_integrand;
    F.params = &params;
    gsl_integration_qag(&F, 0.0, Lambda, 1e-8, 1e-8, WORKSPACE_SIZE, GSL_INTEG_GAUSS15, w, &result, &error);
    gsl_integration_workspace_free(w);
    return (mu - mu_tilde) / (2.0 * Gv) - Nf * Nc / (M_PI * M_PI) * result;
}

// ------------------
// solve_system: Solve dOmega/dM = 0 and dOmega/dmu_tilde = 0
// ------------------
// For simplicity, we implement a basic Newton-Raphson method for a 2x2 system.
int solve_system(double T, double mu, double* M_sol, double* mu_tilde_sol) {
    // Two initial guesses
    double guesses[2][2] = { {400.0, mu}, {m, 0.6 * mu} };
    double best_M = guesses[0][0], best_mu_tilde = guesses[0][1];
    double best_val = Omega_temp(T, mu, best_M, best_mu_tilde);
    
    for (int k = 0; k < 2; ++k) {
        double M = guesses[k][0];
        double mu_tilde = guesses[k][1];
        double d = H; // step size for finite differences
        for (int iter = 0; iter < 100; iter++) {
            double F[2];
            F[0] = dOmegadM(T, mu, M, mu_tilde);
            F[1] = dOmegadmutilde(T, mu, M, mu_tilde);
            if (fabs(F[0]) < d && fabs(F[1]) < d)
                break;
            // Estimate Jacobian using finite differences:
            double J[2][2];
            J[0][0] = (dOmegadM(T, mu, M + d, mu_tilde) - F[0]) / d;
            J[0][1] = (dOmegadM(T, mu, M, mu_tilde + d) - F[0]) / d;
            J[1][0] = (dOmegadmutilde(T, mu, M + d, mu_tilde) - F[1]) / d;
            J[1][1] = (dOmegadmutilde(T, mu, M, mu_tilde + d) - F[1]) / d;
            // Solve the 2x2 system J * dx = -F:
            double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
            if (fabs(det) < 1e-12) break;
            double dx = (-J[1][1] * F[0] + J[0][1] * F[1]) / det;
            double dy = ( J[1][0] * F[0] - J[0][0] * F[1]) / det;
            M += dx;
            mu_tilde += dy;
        }
        double current_val = Omega_temp(T, mu, M, mu_tilde);
        if (current_val < best_val) {
            best_val = current_val;
            best_M = M;
            best_mu_tilde = mu_tilde;
        }
    }
    *M_sol = best_M;
    *mu_tilde_sol = best_mu_tilde;
    return 0;
}


// Define struct for integration parameters
struct num_dens_params {
    double T;
    double mu;
    double M;
};

// Function to integrate
double num_dens_integrand(double p, void *p_ptr) {
    struct num_dens_params *params = (struct num_dens_params *) p_ptr;
    double T = params->T;
    double mu = params->mu;
    double M = params->M;
    double Ep = sqrt(p*p + M*M);
    double term;

    if (T > 0.0) {
        term = (np(p, T, mu, M) - np_bar(p, T, mu, M)) * p*p;
    } else {
        term = ((mu - Ep) >= 0.0 ? 1.0 : 0.0) * p*p;
    }
    return term;
}

// Compute quark number density using GSL integration
double num_dens(double T, double mu, double M) {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(WORKSPACE_SIZE);
    struct num_dens_params params = {T, mu, M};
    double result, error;

    gsl_function F;
    F.function = &num_dens_integrand; // Correct function pointer
    F.params = &params;

    gsl_integration_qag(&F, 0.0, Lambda, 1e-8, 1e-8, WORKSPACE_SIZE, GSL_INTEG_GAUSS15, w, &result, &error);
    gsl_integration_workspace_free(w);

    return Nf * Nc / (M_PI * M_PI) * result;
}

double Omega(double T, double mu){
    double M_sol, mu_tilde_sol;
    double M_vac, mu_tilde_vac;
    solve_system(T, mu, &M_sol, &mu_tilde_sol);
    solve_system(0.0, 0.0, &M_vac, &mu_tilde_vac);
    return Omega_temp(T, mu, M_sol, mu_tilde_sol) - Omega_temp(0.0, 0.0, M_vac, mu_tilde_vac);
}


// ------------------
// main: Example usage
// ------------------
int main() {
    // Initialize coupling constants
    Lambda = 587.9;
    Gs = 2.44 / pow(Lambda, 2.0);
    Gv = Gs * 0.5;
    m = 5.6;

    // Define temperature and chemical potential ranges
    double T_start = 0.0, T_end = 250.0;
    double T_step = (T_end - T_start) / (NUM_SETS - 1);
    double mu_start = 0.0, mu_end = 500.0;
    double mu_step = (mu_end - mu_start) / (NUM_SETS - 1);
    double M_sol, mu_tilde_sol;

    char filename[100];
    char formatted_GvGs[50];
    char formatted_m[50];

    // Gv / Gs と m を文字列にフォーマット
    snprintf(formatted_GvGs, sizeof(formatted_GvGs), "%.1e", Gv / Gs);
    snprintf(formatted_m, sizeof(formatted_m), "%.1f", m);

    // 少数点をアンダースコアに置き換え
    for (int i = 0; formatted_GvGs[i] != '\0'; i++) {
        if (formatted_GvGs[i] == '.') {
           formatted_GvGs[i] = '_';
        }
    }
    for (int i = 0; formatted_m[i] != '\0'; i++) {
        if (formatted_m[i] == '.') {
            formatted_m[i] = '_';
        }
    }

    // 最終的なファイル名を生成
    snprintf(filename, sizeof(filename), "njl_data_GvGs=%s_m=%s_gsl.txt", formatted_GvGs, formatted_m);

    // Open file for writing
    FILE *outfile = fopen(filename, "w");
    if (outfile == NULL) {
        perror("Error opening file");
        return 1;
    }
    fprintf(outfile, "T (MeV), mu (MeV), M solution (MeV), mu_tilde solution (MeV), Baryon density (RHO0)\n");

    // Loop over temperature and chemical potential values
    for (int k = 0; k < NUM_SETS; ++k) {
        double T = T_start + k * T_step;
        for (int i = 0; i < NUM_SETS; ++i) {
            double mu = mu_start + i * mu_step;

            solve_system(T, mu, &M_sol, &mu_tilde_sol);
            double density = num_dens(T, mu_tilde_sol, M_sol) * pow(HBARC, -3.0) / (3.0 * RHO0);

            // Write data to file
            fprintf(outfile, "%lf\t%lf\t%lf\t%lf\t%lf\n", T, mu, M_sol, mu_tilde_sol, density);
        }
    }

    // Close file and print confirmation message
    printf("Results written to %s\n", filename);
    fclose(outfile);

    return 0;
}
