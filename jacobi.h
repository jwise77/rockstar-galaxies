#ifndef _JACOBI_H_
#define _JACOBI_H_

void calc_deviations(double corr[][6], double *sig_x, double *sig_v, double *axis_x, double *axis_v);
void inv_matrix_multiply(double m[][3], double *in, double *out);
void jacobi_decompose(double cov_matrix[][3], double *eigenvalues, double orth_matrix[][3]);
void identity(double val, double *eigenvalues, double orth_matrix[][3]);
#endif /* _JACOBI_H_ */
