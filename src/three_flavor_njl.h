#ifndef THREE_FLAVOR_NJL_H
#define THREE_FLAVOR_NJL_H

#include <gsl/gsl_integration.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_errno.h>

// Global constants (adjust as needed)
extern const double Nf;     // Number of flavors
extern const double Nc;     // Number of colors
extern const double HBARC;
extern const double RHO0;
extern const double H;
extern const int WORKSPACE_SIZE;
extern const int NUM_SETS;

extern double Lambda;  // Momentum cutoff (MeV)
extern double m[3];    // Tormat is {m_u, m_d, m_s}.
extern double G;       // Coupling const
extern double K;       // Coupling const

// Fermi occupation functions
double np_f(double p, double T, double mu, double M);
double np_f_bar(double p, double T, double mu, double M);

// Struct definitions
struct OmegaMf_params {
    double T;
    double mu_f;
    double Mf;
};

struct phi_f_params {
    double T;
    double mu_f;
    double Mf;
};

struct gap_eqs_params {
    double T;
    double mu_u;
    double mu_d;
    double mu_s;
};

struct phi_eqs_params {
    double T;
    double mu_u;
    double mu_d;
    double mu_s;
};

struct M_phi_eqs_params {
    double M_u;
    double M_d;
    double M_s;
};

// Function declarations
double OmegaMf_integrand(double p, void* params_ptr);
double OmegaMf(double T, double mu_f, double Mf);
double Omega_temp(double T, double mu[3], double phi[3]);
double phi_f_integrand(double p, void* params_ptr);
double phi_f(double T, double mu_f, double Mf);
int gap_eqs(const gsl_vector* v, void* params, gsl_vector* f);
double solv_gap_eq(double T, double mu[3], double (*M)[3]);
int phi_eqs(const gsl_vector* v, void* params, gsl_vector* f);
double solv_phi_eq(double T, double mu[3], double (*phi)[3]);
double Omega(double T, double mu[3]);
double pressure(double T, double mu[3]);
int calc_M(double phi_sol[3], double(* M_sol)[3]);
int calc_phi(double M_sol[3], double(* phi_sol)[3]);
double solv_gap_eq_multi_guess(double T, double mu[3], double(* M)[3]);
int calc_phi_multi_guess(double T, double mu[3], double M_sol[3], double(* phi_sol)[3]);

#endif /* THREE_FLAVOR_NJL_H */
