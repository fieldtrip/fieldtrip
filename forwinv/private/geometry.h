double dot(double* a, double* b);
void cross(double* a, double* b, double* r);
double determinant(double* a, double* b, double* c);
double pdist(double* v1);
double ppdist(double* v1, double* v2);
double plinproj(double* l1, double* l2, double* r, double *proj, int flag);
double ptriproj(double* v1, double* v2, double* v3, double* r, double* proj, int flag);
void ltrisect(double* v1, double* v2, double* v3, double* l1, double* l2, double* proj);
void lmoutr(double* v1, double* v2, double* v3, double* r, double *la, double *mu, double *d);
void routlm(double* v1, double* v2, double* v3, double la, double mu, double* r);
int ptriside(double* v1, double* v2, double* v3, double *r);
double solang(double* v1, double* v2, double* v3, int *on_triangle);

