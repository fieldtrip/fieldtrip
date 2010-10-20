#ifndef __MultiChannelFilter_h
#define __MultiChannelFilter_h

template <typename T> 
class MultiChannelFilter {
	public:
	
	// Create filter model
	MultiChannelFilter(int nChans, int order);
	~MultiChannelFilter();
	
	// Calculate Butterworth lowpass filter coefficients
	// from normalised cutoff frequency (TODO)
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
	
}

#endif
