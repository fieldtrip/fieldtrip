#ifndef __MultiChannelFilter_h
#define __MultiChannelFilter_h

#include <math.h>

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

template <typename T> 
class MultiChannelFilter {
	public:
	
	// Create filter model
	MultiChannelFilter(int nChans, int order);
	~MultiChannelFilter();
	
	// Calculate Butterworth lowpass filter coefficients
	// from normalised cutoff frequency (0<cutoff<1=Fnyquist)
	void setButterLP(double cutoff); 
	
	// Copy coefficients. A and B must point to 1+order doubles
	void setCoefficients(const double *B, const double *A) {
		if (A[0]==1.0) {
			for (int i=0;i<=order;i++) {
				this->B[i] = B[i];
				this->A[i] = A[i];
			}
		} else {
			for (int i=0;i<=order;i++) {
				this->B[i] = B[i]/A[0];
				this->A[i] = A[i]/A[0];
			}
		}
	}
	
	// Run single input sample through filter without computing
	// output. This is mostly useful for downsampling purposes.
	void process(const T *source) {
		for (int j=0;j<=order;j++) {
			z[j] = states + pos*nChans;
			if (++pos>order) pos=0;
		}
	
		for (int i=0;i<nChans;i++) z[0][i] = source[i];
		for (int j=1;j<=order;j++) {
			add_scalar_vector(z[0], -A[j], z[j], nChans);
		}
		if (--pos < 0) pos=order;
	}
	
	// Run single input sample through filter and write output
	// to "dest". Both source and dest must point to an nChans-array of T.
	void process(T *dest, const T *source) {
		process(source);
		set_scalar_vector(dest, B[0], z[0], nChans);
		for (int j=1;j<=order;j++) {
			add_scalar_vector(dest, B[j], z[j], nChans);
		}
	}
	
	// Process multiple samples and write output to "dest"
	void process(int nSamples, T *dest, const T *source) {
		for (int i=0;i<nSamples;i++) {
			process(dest + i*nChans, source + i*nChans);
		}
	}
	
	// Clear internal filter states (all zero)
	void clear();
	
	// helper functions with loop-unrolling
	static void set_scalar_vector(T *y, double a,const T *x,int n);
	static void add_scalar_vector(T *y, double a,const T *x,int n);
	
	protected:
	
	T **z;
	T *states;
	double *B, *A;
	int order, nChans, pos;
};

template <typename T> void MultiChannelFilter<T>::set_scalar_vector(T *y, double a,const T *x,int n) {
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

template <typename T> void MultiChannelFilter<T>::add_scalar_vector(T *y, double a,const T *x,int n) {
   /*
   DAXPY_SSE2(X,n,a,x,y);
   */
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

template <typename T> 
MultiChannelFilter<T>::MultiChannelFilter(int nChans, int order) {
	this->nChans = nChans;
	this->order = order;
	this->pos = 0;

	states = new T[nChans*(1+order)];
	clear();
	B = new double[1+order];
	A = new double[1+order];
	z = new T*[1+order];
	B[0] = 1.0;
	A[0] = 1.0;
}

template <typename T> 
MultiChannelFilter<T>::~MultiChannelFilter() {
	delete[] states;
	delete[] A;
	delete[] B;
	delete[] z;
}

template <typename T> 
void MultiChannelFilter<T>::clear() {
	for (int i=0;i<nChans*(1+order);i++) states[i]=0.0;
	pos = 0;
}

template <typename T> 
void MultiChannelFilter<T>::setButterLP(double cutoff) {
	int n;
		
	// warping factor for bilinear transform
	// number inside tan( ) is pi/2
	double f = 1.0/tan(0.5*M_PI*cutoff);

	// nominator is sth. like prodB0*[1 4 6 4 1] later
	double prodB0;

	// init coefficients with unit response
	A[0] = 1.0;
	B[0] = 1.0;
	for (int i=1;i<=order;i++) A[i]=B[i]=0.0;
	
	// for safety: very high cutoff-frequency (almost Nyquist) -> unit response
	if (cutoff > 0.95) return;
	
	// if odd order, handle 1. pole separately
	if (order & 1) {
		prodB0 = 1.0/(1.0+f);
		A[1] = (1-f)/(1.0+f);
		B[1] = 1.0;
		n=1;
	} else {
		prodB0 = 1.0;
		n=0;
	}	

	// add 2 poles at a time (complex conjugates)
	for (int i=n;i<order;i+=2) {
		// location of pole on analog unit circle
		double ang = M_PI * (1.0 - (double)(i+1)/(double)(2*order));
		// analog twopole denominator => (1 + q*s + s^2)
		double q = -2.0*cos(ang);
	
		// bilinear transformation 
		// s -> f*(z+1)/(z-1)
		// 
		// 1 + 2z^-1 + z-^2
		// --------------------------------------------
		// (1+qf+f^2) + (2-2f^2)*z-^1 + (1-qf+f^2)*z-^2
		double b0 = 1.0/(1.0 + q*f + f*f);
		//double b1 = 2.0*b0; 
		//double b2 = b0;
		//double a0 = 1.0;
		double a1 = (2.0-2.0*f*f)*b0;
		double a2 = (1.0-q*f+f*f)*b0;
		// convolve A by [a0 a1 a2] and B by b0*[1 2 1]
		// we can do this in place if we start at the back
		for (int j=i;j>=0;j--) {
			A[j+2] += a2*A[j];
			A[j+1] += a1*A[j];
			B[j+2] += B[j];
			B[j+1] += 2*B[j];
		}
		prodB0 *= b0;
	}
	for (int i=0;i<=order;i++) B[i] *= prodB0;
}

#endif
