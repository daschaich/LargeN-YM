float *vector(int nl, int nh);
int *ivector(int nl,int nh);
double *dvector(int nl,int nh);
float **matrix(int nrl,int nrh,int ncl,int nch);
double **dmatrix(int nrl,int nrh,int ncl,int nch);
int **imatrix(int nrl,int nrh,int ncl,int nch);
float **submatrix(float **a,int oldrl,int oldrh,int oldcl,int oldch,int newrl,int newcl);
void free_vector(float *v,int nl,int nh);
void free_ivector(int *v,int nl,int nh);
void free_dvector(double *v,int nl,int nh);
void free_matrix(float **m,int nrl,int nrh,int ncl,int nch);
void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch);
void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch);
void free_submatrix(float **b,int nrl,int nrh,int ncl,int nch);
float **convert_matrix(float *a,int nrl,int nrh,int ncl,int nch);
void free_convert_matrix(float **b,int nrl,int nrh,int ncl,int nch);
void ludcmp(float **a,int n,int *indx, float *d);
void lubksb(float **a, int n, int *indx, float *b);
void polcoe(double *x,double *y,int n,double *cof);
void polcof(double *xa,double *ya,int n,double *cof);
void polint(double *xa,double *ya,int n,double x,double *y,double *dy);


