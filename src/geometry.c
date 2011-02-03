/*
 * dot		dot product
 * cross	cross product
 * determinant	determinant of matrix built from three vectors
 * pdist	distance from point towards origin
 * ppdist	distance between two points
 * plinproj	projection of point onto line spanned by two points 
 * ptriproj	projection of point onto plane spanned by three points
 * ltrisect	intersection of line with plane spanned by three points
 * lmoutr	compute la/mu parameters for projection on triangle
 * routlm	compute projection on triangle from la/mu parameters
 * ptriside	determines on which side of a triangle a point lies
 * solang       computes solid angle of triangle as seen from origin
 * 
 *  (c) 2002 Robert Oostenveld
 *
 */

#include <math.h>
#include "matrix.h"
#include "geometry.h"

/****************************************************************************/
double dot(double* a, double* b)
{
  double ss=0;
  ss += a[0]*b[0];
  ss += a[1]*b[1];
  ss += a[2]*b[2];
  return ss;
}

/****************************************************************************/
void cross(double* a, double* b, double* r)
{
  r[0] = a[1]*b[2] - a[2]*b[1];
  r[1] = a[2]*b[0] - a[0]*b[2];
  r[2] = a[0]*b[1] - a[1]*b[0];
  return;
}

/****************************************************************************/
double determinant(double* a, double* b, double* c)
{
  return a[0]*(b[1]*c[2]-c[1]*b[2]) - a[1]*(b[0]*c[2]-c[0]*b[2]) + a[2]*(b[0]*c[1]-c[0]*b[1]);
}

/****************************************************************************/
double pdist(double* v1)
{
  double ss=0;
  ss += v1[0]*v1[0];
  ss += v1[1]*v1[1];
  ss += v1[2]*v1[2];
  return sqrt(ss);
}

/****************************************************************************/
double ppdist(double* v1, double* v2)
{
  double ss=0;
  ss += (v1[0]-v2[0])*(v1[0]-v2[0]);
  ss += (v1[1]-v2[1])*(v1[1]-v2[1]);
  ss += (v1[2]-v2[2])*(v1[2]-v2[2]);
  return sqrt(ss);
}

/****************************************************************************/
double plinproj(double* l1, double* l2, double* r, double *proj, int flag)
{
  double l12[3], l1r[3];
  double l12_l=0, l1r_l=0, la=0;
  int k;
  for (k=0; k<3; k++)
  {
    l12[k]=l2[k]-l1[k];
    l1r[k]=r[k]-l1[k];
    l12_l += l12[k]*l12[k];
    l1r_l += l1r[k]*l1r[k];
  }
  l12_l = sqrt(l12_l);
  l1r_l = sqrt(l1r_l);
  if (l12_l==0)
  {
    /* begin and endpoint of linepiece are the same */
    proj[0] = l1[0];
    proj[1] = l1[1];
    proj[2] = l1[2];
    return l1r_l;
  }
  if (l1r_l==0)
  {
    /* point lies exactly at beginpoint of linepiece */
    proj[0] = l1[0];
    proj[1] = l1[1];
    proj[2] = l1[2];
    return 0;
  }
  /* compute relative distance along linepiece */
  la = dot(l12,l1r)/(l12_l*l12_l);
  if (flag && la<0)
    /* intended projection point would lie before begin of linepiece */
    la=0;
  else if (flag && la>1)
    /* intended projection point would lie after end of linepiece */
    la=1;
  /* compute projection point along line through l1 and l2 */
  proj[0] = l1[0] + la*l12[0];
  proj[1] = l1[1] + la*l12[1];
  proj[2] = l1[2] + la*l12[2];
  return ppdist(proj,r);
}

/****************************************************************************/
double ptriproj(double* v1, double* v2, double* v3, double* r, double* proj, int flag)
{
  double la, mu, d;
  lmoutr(v1, v2, v3, r, &la, &mu, &d);
  if (flag)
  {
    /* projection of r on triangle */
    if (la>=0 && mu>=0 && (la+mu)<=1)
      routlm(v1, v2, v3, la, mu, proj);
    else if (la<0)
      d = plinproj(v1, v3, r, proj, flag);
    else if (mu<0)
      d = plinproj(v1, v2, r, proj, flag);
    else
      d = plinproj(v2, v3, r, proj, flag);
  }
  else
  {
    /* projection of r on plane */
    routlm(v1, v2, v3, la, mu, proj);
  }
  return d;
}

/****************************************************************************/
void ltrisect(double* v1, double* v2, double* v3, double* l1, double* l2, double* proj)
{
  double la1, mu1, d1;
  double la2, mu2, d2;
  double p1[3];
  double p2[3];
  int s1, s2;
  char str[256];
  s1 = ptriside(v1, v2, v3, l1);
  s2 = ptriside(v1, v2, v3, l2);
  if (s1==0)
  {
    /* point l1 lies on plane spanned by triangle */
    proj[0] = l1[0];
    proj[1] = l1[1];
    proj[2] = l1[2];
    return;
  }
  else if (s2==0)
  {
    /* point l2 lies on plane spanned by triangle */
    proj[0] = l2[0];
    proj[1] = l2[1];
    proj[2] = l2[2];
    return;
  }
  lmoutr(v1, v2, v3, l1, &la1, &mu1, &d1);
  lmoutr(v1, v2, v3, l2, &la2, &mu2, &d2);
  routlm(v1, v2, v3, la1, mu1, p1);	/* projection of l1 on plane */
  routlm(v1, v2, v3, la2, mu2, p2);	/* projection of l2 on plane */
  /* d1 and d2 are the distances from l1 and l2 to the plane spanned by v1, v2, v3 */
  /* these are used to weigh both projected points to their weighed geometric mean */
  d1 *= s1;
  d2 *= s2;
  if (d1==d2)
  {
    /* the line is parallel to the plane and does not intersect */
    proj[0] = mxGetNaN();
    proj[1] = mxGetNaN();
    proj[2] = mxGetNaN();
    return;
  }
  proj[0] = (d1*p2[0] - d2*p1[0])/(d1-d2);
  proj[1] = (d1*p2[1] - d2*p1[1])/(d1-d2);
  proj[2] = (d1*p2[2] - d2*p1[2])/(d1-d2);
  return;
}

/****************************************************************************/
void lmoutr(double* v1, double* v2, double* v3, double* r, double *la, double *mu, double *ze)
{
  double a[3], b[3], c[3], d[3];
  double a_l=0, b_l=0, c_l=0, d_l=0;
  double det, det_la, det_mu, det_ze;
  int k;

  for (k=0; k<3; k++)
  {
    a[k] = r[k]-v1[k];
    b[k] = v2[k]-v1[k];
    c[k] = v3[k]-v1[k];
    a_l += a[k]*a[k];
    b_l += b[k]*b[k];
    c_l += c[k]*c[k];
  }
  a_l = sqrt(a_l);
  b_l = sqrt(b_l);
  c_l = sqrt(c_l);
  
  if (a_l==0)
  {
    /* point lies exactly on the first vertex */
    *la = 0;
    *mu = 0;
    *ze = 0;
    return;
  }
  else if (b_l==0 || c_l==0)
  {
    /* this is a degenerate triangle, not possible to compute la/mu */
    *la = mxGetNaN();
    *mu = mxGetNaN();
    *ze = mxGetNaN();
    return;
  }
  else
  {
    /* compute vector d orthogonal to the triangle */
    cross(b, c, d); 
    d_l = pdist(d);
    /* normalize vector d */
    d[0] = d[0]/d_l;
    d[1] = d[1]/d_l;
    d[2] = d[2]/d_l;
    /* solve the system of equations using Cramer's method (method of determinants) */
    /* c = la*a + mu*b + ze*o */
    det = determinant(b, c, d);
    det_la = determinant(a, c, d);
    det_mu = determinant(b, a, d);
    det_ze = determinant(b, c, a);
    *la = det_la/det;
    *mu = det_mu/det;
    *ze = det_ze/det;
    /* since d is an orthonormal vector to the plane, the distance of r */
    /* to the plane equals the absolume value of parameter ze */
    *ze = ((*ze)<0 ? -(*ze) : +(*ze));
  }
  return;
}

/****************************************************************************/
void routlm(double* v1, double* v2, double* v3, double la, double mu, double* r)
{
  r[0] = (1-la-mu)*v1[0]+la*v2[0]+mu*v3[0];
  r[1] = (1-la-mu)*v1[1]+la*v2[1]+mu*v3[1];
  r[2] = (1-la-mu)*v1[2]+la*v2[2]+mu*v3[2];
  return;
}

/****************************************************************************/
int ptriside(double* v1, double* v2, double* v3, double *r)
{
  double val;
  double a[3], b[3], c[3], d[3];
  int k;
  for (k=0; k<3; k++)
  {
    a[k] = r[k]-v1[k];
    b[k] = v2[k]-v1[k];
    c[k] = v3[k]-v1[k];
  }
  cross(b, c, d);
  val = dot(a, d);
  if (val>0)
    return 1;
  else if (val<0)
    return -1;
  else
    return 0;
}

/****************************************************************************/
double solang(double *r1, double *r2, double *r3, int *on_triangle)
{
  double cp23_x,cp23_y,cp23_z,n1,n2,n3,ip12,ip23,ip13,nom,den;
  *on_triangle=0;
  cp23_x = r2[1] * r3[2] - r2[2] * r3[1];
  cp23_y = r2[2] * r3[0] - r2[0] * r3[2];
  cp23_z = r2[0] * r3[1] - r2[1] * r3[0];
  nom = cp23_x*r1[0] + cp23_y*r1[1] + cp23_z*r1[2];
  n1 = sqrt(r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2]);
  n2 = sqrt(r2[0]*r2[0] + r2[1]*r2[1] + r2[2]*r2[2]);
  n3 = sqrt(r3[0]*r3[0] + r3[1]*r3[1] + r3[2]*r3[2]);
  ip12 = r1[0]*r2[0] + r1[1]*r2[1] + r1[2]*r2[2];
  ip23 = r2[0]*r3[0] + r2[1]*r3[1] + r2[2]*r3[2];
  ip13 = r1[0]*r3[0] + r1[1]*r3[1] + r1[2]*r3[2];
  den = n1*n2*n3 + ip12*n3 + ip23*n1 + ip13*n2;
  if (nom==0 & den<=0)
  {
    *on_triangle=1;
    return 0;
  }
  return -2.*atan2 (nom,den);
}
