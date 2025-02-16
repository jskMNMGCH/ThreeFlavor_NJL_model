#include "three_flavor_njl.h"
#include <stdio.h>
#include <gsl/gsl_errno.h>

void initialize_parameters(int param_set) {
    if (param_set == 0) {
        Lambda = 602.3;  // MeV
        m[0] = 5.5;      // Current quark mass for u (MeV)
        m[1] = 5.5;      // Current quark mass for d (MeV)
        m[2] = 140.7;    // Current quark mass for s (MeV)
        G = 1.835*pow(Lambda, -2.0);       // Coupling constant G
        K = 12.36*pow(Lambda, -5.0);       // Coupling constant K
    }
    else if (param_set == 1) {
        Lambda = 631.4;  // MeV
        m[0] = 5.5;      // Current quark mass for u (MeV)
        m[1] = 5.5;      // Current quark mass for d (MeV)
        m[2] = 135.7;    // Current quark mass for s (MeV)
        G = 1.835*pow(Lambda, -2.0);       // Coupling constant G
        K = 9.29*pow(Lambda, -5.0);       // Coupling constant K
    }
}

int main() {
    initialize_parameters(1);

    double T = 1.0;  // Temperature (MeV)
    double mu_start = 2e2;  // Starting chemical potential (MeV)
    double mu_end = 6e2;   // Ending chemical potential (MeV)
    double mu_step = 4.0;    // Step size for chemical potential (MeV)

    double phi_sol[3];
    double M_sol[3];

    char filename[256];
    char formatted_GL2[50];
    char formatted_KL5[50];

    // G, K を文字列にフォーマット
    snprintf(formatted_GL2, sizeof(formatted_GL2), "%.3f", G*pow(Lambda, 2.0));
    snprintf(formatted_KL5, sizeof(formatted_KL5), "%.2f", K*pow(Lambda, 5.0));

    // 少数点をアンダースコアに置き換え
    for (int i = 0; formatted_GL2[i] != '\0'; ++i) {
        if (formatted_GL2[i] == '.') {
           formatted_GL2[i] = '_';
        }
    }
    for (int i = 0; formatted_KL5[i] != '\0'; ++i) {
        if (formatted_KL5[i] == '.') {
            formatted_KL5[i] = '_';
        }
    }

    snprintf(filename, sizeof(filename), "data_three_flavor_njl_GL2=%s_KL5=%s.txt", formatted_GL2, formatted_KL5);

    FILE *output_file = fopen(filename, "w");
    if (output_file == NULL) {
        perror("Error opening file");
        return 1;
    }
    fprintf(output_file, "T(MeV)\t\tmu(MeV)\t\tM_u(MeV)\t\tM_d(MeV)\t\tM_s(MeV)\t\tphi_u/L^3\t\tphi_d/L^3\t\tphi_s/L^3\n");

    int status;
    for (double mu = mu_start; mu <= mu_end; mu += mu_step) {
        double mu_arr[3] = {mu, mu, mu};

        status = solv_gap_eq_multi_guess_solver_restrict(T, mu_arr, &M_sol);
        calc_phi_multi_guess(T, mu_arr, M_sol, &phi_sol);
        fprintf(output_file, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",T,mu,M_sol[0],M_sol[1],M_sol[2],phi_sol[0],phi_sol[1],phi_sol[2]);
    }
    printf("Output file: %s\n", filename);
    fclose(output_file);
    return 0;
}
