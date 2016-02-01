/* ARS.C - Procedure for performing Adaptive Rejection Sampling. */

/* Copyright (c) 1995-2003 by Carl Edward Rasmussen and Radford M. Neal
 *
 * Permission is granted for anyone to copy, use, modify, or distribute this
 * program and accompanying programs and documents for any purpose, provided 
 * this copyright notice is retained and prominently displayed, along with
 * a note saying that the original programs are available from Radford Neal's
 * web page, and note is made of any changes made to the programs.  The
 * programs and documents are distributed without any warranty, express or
 * implied.  As the programs were written for research purposes only, they have
 * not been tested to the degree that would be advisable in any important
 * application.  All use of these programs is entirely at the user's own risk.
 */


/* This module implements the "Adaptive Rejection Sampling" scheme due to
 * Gilks and Wild.  See "Adaptive Rejection Sampling for Gibbs Sampling",
 * Applied Statistics, vol. 41, no. 2, pp. 337-348 (1992).  This is not
 * the most sophisticated possible implementation of the method, however. */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "rand.h"
#include "ars.h"


#define MIN_DROP  0.1  /* minimum drop in log probability for initial points */
#define TINY      1e-9         /* minimum tolerance for inserting new points */
#define MAX_LIST  100     /* max number of segments in piecewise exponential */

struct segtype { double a;           /* in the log domain, the piece wise... */
                 double b;     /* exponential is a piece wise linear: y=ax+b */
                 double x;                        /* interior point in piece */
                 double xmax;                   /* upper limit of this piece */
                 double mass;          /* the probability mass of this piece */
                 struct segtype *prv;               /* ptr to previous piece */
                 struct segtype *nxt; };                /* ptr to next piece */

static struct segtype *root; 

static double rpwed(struct segtype **);


/* ADAPTIVE REJECTION SAMPLING PROCEDURE.
 *
 * The ars function returns a single point sampled at random from the univariate
 * distribution specified by the function logp, provided by the user.  The
 * logp function takes a point as its first argument, and returns the log 
 * of the probability density at that point, plus any arbitrary constant (not
 * depending on the point).  The logp function takes two additional arguments,
 * the first a pointer to a place where it must store the derivative of the
 * log probability density, the second a pointer to additional information
 * describing the distribution, which is passed on unchanged from the last 
 * argument of ars.
 *
 * The logp function passed MUST be log-concave. It is assumed that any real 
 * number is legal as input to logp. 
 *
 * The user must also supply an initial guess, "init", and a typical "scale" 
 * of variation.  It is not essential that these values be very accurate, but 
 * performance will generally depend on their accuracy.
 *
 * The ars function first tries to locate points on either side of the mode;
 * the derivative must have the right sign and be non-negligible, and the drop 
 * from the max. probability seen must be of at least moderate size to qualify.
 * Then a piece-wise exponential distribution is iteratively improved using 
 * knowledge from rejected points, until a sample is accepted. At most MAX_LIST
 * pieces are introduced in the approximation. Before pieces are inserted 
 * various checks are made in order to prevent numerical problems. If new points
 * don't qualify, the piece-wise exponential is simply not updated. A warning 
 * will be issued (one time only) when 10000 rejections are exceeded */

double ars(
  double init,                                              /* initial guess */
  double scale,                         /* guess for scale of variation in x */
  double (*logp)(double, double *, void *),       /* function to sample from */
  void   *extra                   /* any extra information to pass to logp() */
)
{
  struct segtype seg[MAX_LIST], *prv, *cur, *nxt;
  int    i, no_seg = 2;
  static int warning = 0;
  double x, max, f, df;

  root = &seg[0]; nxt = &seg[1];
  root->prv = NULL; root->nxt = nxt; nxt->prv = root; nxt->nxt = NULL;

  x = init; f = logp(x, &df, extra);       /* find point to the left of mode */
  max = f;
  while (df<TINY || f>max-MIN_DROP) {
    x -= scale+(init-x); 
    f = logp(x, &df, extra);
    if (f>max) max = f;
  }
  root->x = x; root->a = df; root->b = f-x*df;

  x = init; f = logp(x, &df, extra);      /* find point to the right of mode */
  while (df>-TINY || f>max-MIN_DROP) {
    x += scale+(x-init); 
    f = logp(x, &df, extra);
    if (f>max) max = f;
  }
  nxt->x = x; nxt->a = df; nxt->b = f-x*df;
  root->xmax = (nxt->b-root->b)/(root->a-nxt->a);
 
  for (i=0; ; i++) {                     /* repeat until a point is accepted */

    if (i==10000 && !(warning)) {
      fprintf(stderr, "WARNING: More than 10000 rejections in ars\n");
      warning = 1;
    } 

    cur = root;      /* find max y-value; needed to avoid numerical problems */
    max = cur->a*cur->xmax+cur->b;
    while (cur = cur->nxt, cur->nxt) 
      if ((x = cur->a*cur->xmax+cur->b) > max) max = x;

    cur = root;                                            /* compute masses */
    cur->mass = exp(cur->a*cur->xmax+cur->b-max)/cur->a;
    while (prv = cur, cur = cur->nxt, cur->nxt) 
      cur->mass = (exp(cur->a*cur->xmax+cur->b-max)-
                   exp(cur->a*prv->xmax+cur->b-max))/cur->a;
    cur->mass = -exp(cur->a*prv->xmax+cur->b-max)/cur->a;

    x = rpwed(&cur);                               /* this is the new sample */
    f = logp(x, &df, extra);
    if (rand_uniform() <= exp(f-cur->a*x-cur->b)) return x;      /* success! */

/* Now, insert a new piece in the piece-wise approximation if the situation is
 * appropriate. Eg, if we have enough memory, if the slope at the new x isn't
 * too small (the exponential distribution will degenerate), and if the slope
 * isn't too close to the slope of current, previous and next segment (or
 * which ever of these may exist), since this may cause numerical problems. */

    if (no_seg < MAX_LIST && fabs(df) > TINY && fabs(df-cur->a) > TINY &&
      (!(cur->prv) || fabs(df-cur->prv->a) > TINY) &&
      (!(cur->nxt) || fabs(df-cur->nxt->a) > TINY)) {      

      if (x<cur->x) cur = cur->prv;               /* now, insert *after* cur */
      prv = cur; cur = &seg[no_seg++]; cur->prv = prv;
      if (prv)
        { cur->nxt = prv->nxt; prv->nxt = cur; }
      else
        { cur->nxt = root; root = cur; }
      nxt = cur->nxt; if (nxt) nxt->prv = cur;
      cur->x = x; cur->a = df; cur->b = f-x*df;

      if (prv) prv->xmax = (cur->b-prv->b)/(prv->a-cur->a);
      if (nxt) cur->xmax = (nxt->b-cur->b)/(cur->a-nxt->a);
    }
  }
}


/* Private function to sample from piece-wise exponential distribution. First a
 * piece under the piece-wise distribution is sampled at random. Then a random
 * sample is drawn from this piece.  A pointer to the segment, or piece which
 * was used is returned in q, and the function returns the random sample. Care
 * is taken to avoid numerical over and underflow.  */

static double rpwed(struct segtype **q)
{
  double mass = 0.0, t, u;
  
  *q = root; while (*q) { mass += (*q)->mass;  *q = (*q)->nxt; }
  t = mass*rand_uniform();
  *q = root; while ((*q)->nxt && ((t -= (*q)->mass) >= 0.0)) *q = (*q)->nxt; 

  u = rand_uniopen();
  if ((*q)->prv == NULL)
    return (*q)->xmax+log(u)/(*q)->a;
  if ((*q)->nxt == NULL) 
    return (*q)->prv->xmax+log(u)/(*q)->a;
  t = log(u+(1.0-u)*exp(fabs((*q)->a)*((*q)->prv->xmax-(*q)->xmax)))/(*q)->a;
  return ((*q)->a > 0) ? (*q)->xmax+t : (*q)->prv->xmax+t;
}
