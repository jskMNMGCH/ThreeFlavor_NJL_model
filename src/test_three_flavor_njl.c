#include "three_flavor_njl.h"
#include <stdio.h>
#include <gsl/gsl_errno.h>

void initialize_parameters() {
    Lambda = 602.3;  // MeV
    m[0] = 5.5;      // Current quark mass for u (MeV)
    m[1] = 5.5;      // Current quark mass for d (MeV)
    m[2] = 140.7;    // Current quark mass for s (MeV)
    G = 1.835*pow(Lambda, -2.0);       // Coupling constant G
    K = 12.36*pow(Lambda, -5.0);       // Coupling constant K
}

int main() {
    initialize_parameters();

    double T = 1e-3;  // Temperature (MeV)
    double mu[3] = {300.0, 300.0, 300.0};  // Chemical potentials for u, d, s (MeV)

    int status;
    // Test solv_phi_eq
    double phi_sol[3];
    status = solv_phi_eq(T, mu, &phi_sol);
    
    double M_sol[3];
    calc_M(phi_sol, &M_sol);
    printf("solv_phi_eq: M_u = %lf, M_d = %lf, M_s = %lf\n",M_sol[0],M_sol[1],M_sol[2]);

    solv_gap_eq(T, mu, &M_sol);
    printf("solv_phi_eq: M_u = %lf, M_d = %lf, M_s = %lf\n",M_sol[0],M_sol[1],M_sol[2]);

    return 0;
}
