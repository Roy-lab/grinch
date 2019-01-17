#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_rng.h>


#define min(x,y) (((x) < (y)) ? (x) : (y))
#define max(x,y) (((x) > (y)) ? (x) : (y))


/* write matrix to file 
format:
% comment
num_rows num_columns num_nonzeros
row col nnz
.....
row col nnz
*/
void matrix_write_to_file(gsl_matrix *M, char *fname);


/* load matrix from file 
format:
% comment
num_rows num_columns num_nonzeros
row col nnz
.....
row col nnz
*/
gsl_matrix * matrix_load_from_text_file(char *fname);



/* load matrix from binary file 
 * the nonzeros are in order of double loop over rows and columns
format:
num_rows (int) 
num_columns (int)
nnz (double)
...
nnz (double)
*/
gsl_matrix * matrix_load_from_binary_file(char *fname);


void matrix_write_to_binary_file(gsl_matrix *M, char *fname);


/* load vector from file 
format:
% comment
num_rows
value
.....
value
*/
gsl_vector * vector_load_from_file(char *fname);



/* build up a random matrix R */
void initialize_random_matrix_old(gsl_matrix *M);

void initialize_random_matrix(gsl_matrix *M);

void initialize_random_matrix_with_seed(gsl_matrix *M, int seed);
/*
% project v in direction of u
function p=project_vec(v,u)
p = (dot(v,u)/norm(u)^2)*u;
*/
void project_vector(gsl_vector *v, gsl_vector *u, gsl_vector *p);



/* build orthonormal basis matrix
Q = Y;
for j=1:k
    vj = Q(:,j);
    for i=1:(j-1)
        vi = Q(:,i);
        vj = vj - project_vec(vj,vi);
    end
    vj = vj/norm(vj);
    Q(:,j) = vj;
end
*/
void build_orthonormal_basis_from_mat(gsl_matrix *A, gsl_matrix *Q);



/* C = A*B */
void matrix_matrix_mult(gsl_matrix *A, gsl_matrix *B, gsl_matrix *C);


/* C = A^T*B */
void matrix_transpose_matrix_mult(gsl_matrix *A, gsl_matrix *B, gsl_matrix *C);


/* y = M*x */
void matrix_vector_mult(gsl_matrix *M, gsl_vector *x, gsl_vector *y);


/* y = M^T*x */
void matrix_transpose_vector_mult(gsl_matrix *M, gsl_vector *x, gsl_vector *y);


/* compute evals and evecs of symmetric matrix */
void compute_evals_and_evecs_of_symm_matrix(gsl_matrix *M, gsl_vector *eval, gsl_matrix *evec);


/* compute QR factorization 
M is mxn; Q is mxm and R is mxn
this is slow
*/
void compute_QR_factorization(gsl_matrix *M, gsl_matrix *Q, gsl_matrix *R);


/* compute compact QR factorization 
M is mxn; Q is mxk and R is kxk
*/
void compute_QR_compact_factorization(gsl_matrix *M, gsl_matrix *Q, gsl_matrix *R);



/* compute compact QR factorization and get Q 
M is mxn; Q is mxk and R is kxk (not computed)
*/
void QR_factorization_getQ(gsl_matrix *M, gsl_matrix *Q);



/* build diagonal matrix from vector elements */
void build_diagonal_matrix(gsl_vector *dvals, int n, gsl_matrix *D);



/* invert diagonal matrix */
void invert_diagonal_matrix(gsl_matrix *Dinv, gsl_matrix *D);


/* frobenius norm */
double matrix_frobenius_norm(gsl_matrix *M);


/* print matrix */
void matrix_print(gsl_matrix *M);


/* P = U*S*V^T */
void form_svd_product_matrix(gsl_matrix *U, gsl_matrix *S, gsl_matrix *V, gsl_matrix *P);


/* calculate percent error between A and B: 100*norm(A - B)/norm(A) */
double get_percent_error_between_two_mats(gsl_matrix *A, gsl_matrix *B);


/* copy the first k columns of M into M_out where k = M_out->ncols (M_out pre-initialized) */
void matrix_copy_first_columns(gsl_matrix *M_out, gsl_matrix *M);

/* C = [A, B]; */
void append_matrices_horizontally(gsl_matrix *A, gsl_matrix *B, gsl_matrix *C);

/* for rsvd autorank1 */
void estimate_rank_and_buildQ(gsl_matrix *M, double frac_of_max_rank, double TOL, gsl_matrix **Q, int *good_rank);

/* for rsvd autorank2 */
void estimate_rank_and_buildQ2(gsl_matrix *M, int kblock, double TOL, gsl_matrix **Y, gsl_matrix **Q, int *good_rank);

