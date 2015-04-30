/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */
#ifndef __Brain3dWindow_h
#define __Brain3dWindow_h

#include <math.h>

#include <FL/Fl_Gl_Window.H>

struct LightProps {
   int active;
   float position[4];
   float specular[4];
   float ambient[4];
   float diffuse[4];
   
   LightProps() {
      active = 0;
      position[0]=position[1]=position[2]=position[3]=0.0f;
      specular[0]=ambient[0]=diffuse[0]=0.0f;
      specular[1]=ambient[1]=diffuse[1]=0.0f;
      specular[2]=ambient[2]=diffuse[2]=0.0f;
      specular[3]=ambient[3]=diffuse[3]=1.0f;            
   }
};

class Brain3dWindow : public Fl_Gl_Window {
	public:
	
	Brain3dWindow(int x, int y, int w, int h, const char *label = NULL);
	
	virtual void draw();
	virtual int handle(int event);
	
	void drawBrain();
   
	void setDefault() {
		center[0] = center[1] = center[2] = 0.0;
		theta = 0.6*M_PI;
		phi = 0; // 0.5*M_PI;
		psi = 0;
		dist = 3.0;
	}

	void calcEyeUpRight();   

	const double *getEye() const { return eye; }
	const double *getCenter() const { return center; }
	const double *getUp() const { return up; }
	const double *getRight() const { return right; }
	const double *getBackground() const { return bgcolor; }   

	double getTheta() const { return theta; }
	double getPhi() const { return phi; }
	double getPsi() const { return psi; }   
	double getDist() const { return dist; }
	double getFanSize() const { return fanSize; }
	
	void setWarp(float factor, float offset) {
		warpFactor = factor;
		warpOffset = offset;
	}

	void setDist(double dist) {
		this->dist = dist;
		calcEyeUpRight();
	}

	void setTheta(double theta) {
		this->theta = theta;
		calcEyeUpRight();
	}

	void setPhi(double phi) {
		this->phi = phi;
		calcEyeUpRight();
	}

	void setPsi(double psi) {
		this->psi = psi;
		calcEyeUpRight();
	}

	void setCenter(double x, double y, double z) {
		center[0] = x;
		center[1] = y;
		center[2] = z;
		calcEyeUpRight();
	}

	void setFanSize(double size) {
		if (size>0) {
			fanSize = size; 
		} else {
			fanSize = 0.0;
		}
	}

	void setBackground(double r = 0.0, double g = 0.0, double b = 0.0) {
		bgcolor[0] = r;
		bgcolor[1] = g;
		bgcolor[2] = b;
	}
	
	
	int getLighting() const { return useLight; }	
	void setLighting(int yesno) { useLight = yesno; }
   
	const LightProps *getLight(int i) {
		if (i<0 || i>7) return 0;
		return &light[i]; 
	}
   
	void setLight(int i, const LightProps *l) {
		if (i<0 || i>7 || l==0) return;
		light[i] = *l;
	}
	
	void setSliceTextures(int numSlices, const GLuint *textures) {
		if (numSlices<=0 || textures == NULL) {
			numSlices = 0;
			textures = NULL;
		} else {
			this->numSlices = numSlices;
			this->texture = textures;
		}
	}

	protected:
	
	double eye[3];
	double center[3];
	double up[3];
	double right[3];
	double bgcolor[3];
	double theta,phi,psi;
	double dist;
	double fanSize;
	double warpFactor, warpOffset;
	bool useLight;
	LightProps light[8];
	
	int numSlices;
	const GLuint *texture;	

	private:

	int mouse_x;
	int mouse_y;
	int button;
};

#endif
