/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */

#include <platform.h>
#if defined (PLATFORM_OSX)
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <FL/Fl.H>

#include <Brain3dWindow.h>

Brain3dWindow::Brain3dWindow(int X, int Y, int W, int H, const char *L) : Fl_Gl_Window(X,Y,W,H,L) {

	bgcolor[0] = bgcolor[1] = bgcolor[2] = 0.0;

	light[0].active = 1;
	light[0].position[0] = 1000.0;
	light[0].position[1] = -800.0;
	light[0].position[2] = 1000.0;
	light[0].ambient[0] = 0.3;         
	light[0].ambient[1] = 0.3;         
	light[0].ambient[2] = 0.3;               
	light[0].diffuse[0] = 1.0;
	light[0].diffuse[1] = 1.0;
	light[0].diffuse[2] = 1.0;      

	light[1].active = 1;
	light[1].position[0] = -1000.0;
	light[1].position[1] =   800.0;
	light[1].position[2] =  1000.0;
	light[1].diffuse[0] = 0.7;
	light[1].diffuse[1] = 0.7;
	light[1].diffuse[2] = 0.4;      


	fanSize = 0.0;

	setDefault();
	calcEyeUpRight();

	mode(FL_DOUBLE|FL_DEPTH|FL_ALPHA);
	resizable(this);
	
	numSlices = 0;
	texture = NULL;
}



void Brain3dWindow::drawBrain() {
	if (numSlices<=0) return;
	
	float z_min = 0;
	float dz = 1.0/numSlices;
		
	glColor3f(1,1,1);
	for (int i=0;i<numSlices;i++) {
		float z_org = z_min + i*dz;
		
		float z = -0.7 + 1.4*((1-warpFactor)*z_org + warpFactor/(1+exp(-50*(z_org - warpOffset))));
		
		glBindTexture(GL_TEXTURE_2D, texture[i]);
		
		glBegin(GL_QUADS);
		glTexCoord2f(0,0); glVertex3f(-1,-1,z);
		glTexCoord2f(1,0); glVertex3f(+1,-1,z);
		glTexCoord2f(1,1); glVertex3f(+1,+1,z);
		glTexCoord2f(0,1); glVertex3f(-1,+1,z);
		glEnd();
	}
}
   
   
void Brain3dWindow::draw() {
	int ww = w();
	int wh = h();

	#ifndef MESA
	// glDrawBuffer(GL_FRONT_AND_BACK);
	#endif // !MESA

	glViewport(0, 0, ww, wh);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0f, (double)ww/(double)wh, 0.01, 10000.0);

	glFrontFace(GL_CCW);
	glEnable(GL_TEXTURE_2D);

	if (useLight) {
		glEnable(GL_LIGHTING);      
		glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
		// glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientLight);      

		for (int i=0;i<8;i++) {
			if (light[i].active) {
				glEnable(GL_LIGHT0 + i);
			} else {
				glDisable(GL_LIGHT0 + i);
				continue;
			}
			glLightfv(GL_LIGHT0 + i, GL_POSITION, light[i].position);
			glLightfv(GL_LIGHT0 + i, GL_SPECULAR, light[i].specular);
			glLightfv(GL_LIGHT0 + i, GL_AMBIENT, light[i].ambient);
			glLightfv(GL_LIGHT0 + i, GL_DIFFUSE, light[i].diffuse);
		}
	} else {
		glDisable(GL_LIGHTING);
		glDisable(GL_COLOR_MATERIAL);
	}

	glShadeModel(GL_SMOOTH);
	glPolygonMode(GL_FRONT, GL_FILL);

	glDisable(GL_DEPTH_TEST);
	//glDepthMask(GL_TRUE);

	glEnable (GL_BLEND); 
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glClearColor(bgcolor[0],bgcolor[1],bgcolor[2],0); // clear the window to black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // clear the window

	gluLookAt(eye[0], eye[1], eye[2], 
				center[0], center[1], center[2], 
				up[0], up[1], up[2]);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	if (fanSize>0) {
		glBegin(GL_TRIANGLE_FAN);      
		glColor3f(1,1,1);   glVertex3f(0,0,0);
		glColor3f(1,0,0);   glVertex3f(fanSize,0,0);
		glColor3f(0,1,0);   glVertex3f(0,fanSize,0);
		glColor3f(0,0,1);   glVertex3f(0,0,fanSize);      
		glColor3f(1,0,0);   glVertex3f(fanSize,0,0);      
		glEnd();

		glBegin(GL_TRIANGLES);            
		glColor3f(1,0,0);   glVertex3f(fanSize,0,0);
		glColor3f(0,1,0);   glVertex3f(0,fanSize,0);
		glColor3f(0,0,1);   glVertex3f(0,0,fanSize);      
		glEnd();  
	}

	drawBrain();

	#ifndef MESA
	// glDrawBuffer(GL_BACK);
	#endif // !MESA
}


int Brain3dWindow::handle(int event) {
   if (event == FL_CLOSE) {
      return 0;
   } else if (event == FL_MOUSEWHEEL) {
      if (Fl::event_dy()<0) {
         dist-=0.1;
         if (dist<0.1) dist = 0.1;
      } else {
         dist+=0.1;
      }            
      calcEyeUpRight();         
      redraw();
      return 1;
   } else if (event == FL_PUSH) {
      mouse_x = Fl::event_x();
      mouse_y = Fl::event_y();
      return 1;
   } else if (event == FL_RELEASE) {
      return 1;
   } else if (event == FL_DRAG) {
      double dx = 0.01*(double) (Fl::event_x() - mouse_x);
      double dy = 0.01*(double) (Fl::event_y() - mouse_y);            

      if (Fl::event_state() & FL_BUTTON1) {
         if (Fl::event_state() & FL_CTRL) {
            center[0]+=dy*(eye[0]-center[0]);
            center[1]+=dy*(eye[1]-center[1]);
            center[2]+=dy*(eye[2]-center[2]);
         } else {
            center[0]+=dx*right[0] + dy*up[0];
            center[1]+=dx*right[1] + dy*up[1];
            center[2]+=dx*right[2] + dy*up[2];
         }
         mouse_x = Fl::event_x();
         mouse_y = Fl::event_y();
         calcEyeUpRight();
         redraw();
         return 1;
      } else if (Fl::event_state() & FL_BUTTON3) {
         if (Fl::event_state() & FL_CTRL) {
            psi+=dx;
         } else {
            phi-=dx;
            theta+=dy;

            if (theta<0) theta=0;
            if (theta>M_PI) theta=M_PI;
         } 
         mouse_x = Fl::event_x();
         mouse_y = Fl::event_y();
         calcEyeUpRight();
         redraw();
         return 1;
      }
   } else if (event == FL_KEYDOWN) {
      switch (Fl::event_key()) {
         case FL_Home:
            setDefault();
            calcEyeUpRight();
            redraw();
            return 1;
          default:
            return 0;
      }
   } 
   return 0;  
}


void Brain3dWindow::calcEyeUpRight() {
	double ct = cos(theta);
	double st = sin(theta);
	double cph = cos(phi);
	double sph = sin(phi);
	double cps = cos(psi);
	double sps = sin(psi);


	eye[0] = center[0] + dist*st*cph;
	eye[1] = center[1] + dist*st*sph;
	eye[2] = center[2] - dist*ct;

	/* without psi 
	up[0] = ct*cph;
	up[1] = ct*sph;      
	up[2] = st;

	right[0] = sph;
	right[1] = -cph;
	right[2] = 0;  
	*/   

	up[0] = cps*ct*cph - sps*sph;
	up[1] = cps*ct*sph + sps*cph;      
	up[2] = cps*st;

	right[0] = cps*sph  + sps*ct*cph;
	right[1] = -cps*cph + sps*ct*sph;
	right[2] = 0        + sps*st;
}



