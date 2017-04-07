#ifndef __TemplateVectorMath_h
#define __TemplateVectorMath_h

template <typename Tout, typename Ts, typename Tin>
void tvmSetScaledVector(Tout *y, Ts a, const Tin *x, int n) {
   /* for (i=0;i<n;i++) y[i] = a*x[i]; */
   while (n>=8) {
      y[0] = a*x[0];
      y[1] = a*x[1];
      y[2] = a*x[2];
      y[3] = a*x[3];
      y[4] = a*x[4];
      y[5] = a*x[5];
      y[6] = a*x[6];
      y[7] = a*x[7];
      n-=8;
      y+=8;
      x+=8;
   }
   switch(n) {
      case 7: y[6] = a*x[6];
      case 6: y[5] = a*x[5];
      case 5: y[4] = a*x[4];
      case 4: y[3] = a*x[3];
      case 3: y[2] = a*x[2];
      case 2: y[1] = a*x[1];
      case 1: y[0] = a*x[0];
   }
}

template <typename Tout, typename Ts, typename Tin>
void tvmAddScaledVector(Tout *y, Ts a, const Tin *x, int n) {
   /* for (i=0;i<n;i++) y[i] += a*x[i]; */
   while (n>=8) {
      y[0] += a*x[0];
      y[1] += a*x[1];
      y[2] += a*x[2];
      y[3] += a*x[3];
      y[4] += a*x[4];
      y[5] += a*x[5];
      y[6] += a*x[6];
      y[7] += a*x[7];
      n-=8;
      y+=8;
      x+=8;
   }
   switch(n) {
      case 7: y[6] += a*x[6];
      case 6: y[5] += a*x[5];
      case 5: y[4] += a*x[4];
      case 4: y[3] += a*x[3];
      case 3: y[2] += a*x[2];
      case 2: y[1] += a*x[1];
      case 1: y[0] += a*x[0];
   }
}

#endif
