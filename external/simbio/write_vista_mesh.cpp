#include <mex.h>

#include "FIL_Vista.h"
#include <stdio.h>
#include <string>
#include <assert.h>

unsigned int countNodes(VGraph graph)
{
  unsigned int n = 0;
  for (int i=1; i<=graph->lastUsed; i++)
    if (VGraphGetNode(graph, i)) n++;
  return n;
}

class primitive : public VNodeBaseRec  {       
public:
        VLong        vcnt;
        VLong        id[8];
        primitive()     { hops = 0; visited = 0; head = 0; weight = 0; }
};

class vertex : public VNodeBaseRec  {
public:
	VFloat		type;
	VFloat		x,y,z;
    VUByte      color;

 	vertex()	{ hops = 0; visited = 0; head = 0; weight = 0; };

 	vertex(double a, double b, double c)
       			{ hops = 0; visited = 0; head = 0; weight = 0;
			  type = 1; x = (VFloat)a; y = (VFloat)b; z = (VFloat)c; };
};

// main function
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   VUByte flag=0;
   
/*   static VOptionDescRec options[] = {
            { "conv", VUByteRepn, 1, &flag, VOptionalOpt, NULL, 
   		 	"Input node numbering in C-convention: 0 -- " }
   };*/
    
   //VParseFilterCmd( VNumber(options), options, argc, argv, &inf, &outf );
   
   // check number of arguments
   int nargin = nrhs;
   /*if((nargin != 4)&(nargin != 5))
   {
       mexErrMsgTxt("Not enough arguments. Usage: <filename>, <nodes>, <elements> (,<labels)\n");
   }
   else if(!mxIsChar(prhs[0]))
   {
       mexErrMsgTxt("Wrong type of argument. <filename> must be of type 'char'.");
   }
   else if(!mxIsDouble(prhs[1]))
   {
       mexErrMsgTxt("Wrong type of argument. <nodes> must be of type 'double'.");
   }
   else if(!mxIsDouble(prhs[2]))
   {
       mexErrMsgTxt("Wrong type of argument. <elements> must be of type 'double'.");
   }
   else if(!mxIsDouble(prhs[3]))
   {
       mexErrMsgTxt("Wrong type of argument. <labels> must be of type 'double'.");
   }
   else if((!mxIsDouble(prhs[4]))&(nargin==5))
   {
       mexErrMsgTxt("Wrong type of argument. <tensors> must be of type 'double'.");
   }*/
   
   const mxArray *mxelements = prhs[2];
   const mxArray *mxnodes = prhs[1];
   const mxArray *mxlabels = prhs[3];
   const mxArray *mxtensors = NULL;
   double *elements = mxGetPr(mxelements);
   double *nodes = mxGetPr(mxnodes);
   double *labels = mxGetPr(mxlabels);
   double *tensors = NULL;
   int rowLab = mxGetN(mxlabels);
   int colLab = mxGetM(mxlabels);
   int rowEle = mxGetN(mxelements); //number of elements
   int colEle = mxGetM(mxelements);
   int rowNod = mxGetN(mxnodes); //number of nodes
   int colNod = mxGetM(mxnodes);
   int rowTen = 0;
   int colTen = 0;
     
   if(nargin == 5){
      mxtensors = prhs[4];
      tensors = mxGetPr(mxtensors);
      rowTen = mxGetN(mxtensors);
      colTen = mxGetM(mxtensors);
   }
   
   if(colLab != colEle)
       mexErrMsgTxt("Number of elements must fit number of labels.");
   
   if((colTen != colEle)&&(nargin == 5))
       mexErrMsgTxt("Number of elements must fit number of tensors.");
   
   VGraph vertices, primitives;
   VImage mp = 0, AttributeImage = 0;                     // 1D image storing the material properties
   int id[8];
    
   // get filename from arguments
   std::string filename(mxArrayToString(prhs[0]));
   
   
   FILE *outf;
   outf = fopen (filename.c_str(),"w");

   vertices = VCreateGraph(colNod, 5, VFloatRepn, false);
    if (vertices == NULL)  {
      mexErrMsgTxt("Cannot create vertices graph.\n");
    }

    VSetAttr(VGraphAttrList(vertices), "patient", NULL, VStringRepn, "Matlab to Vista");
    VSetAttr(VGraphAttrList(vertices), "date", NULL, VStringRepn, "27.4.2011");
    VSetAttr(VGraphAttrList(vertices), "convention", NULL, VStringRepn, "natural");
//    VSetAttr(VGraphAttrList(vertices), "voxel", NULL, VStringRepn, "1.000000 1.000000 1.000000");
    VSetAttr(VGraphAttrList(vertices), "component_interp", NULL, VStringRepn, "vertex"   );

    for(int i=0; i<colNod; i++){

		double x=(double)nodes[i]; //static_cast<double>(nodes[i])
		double y=(double)nodes[colNod + i]; //static_cast<double>
		double z=(double)nodes[2*colNod + i]; //static_cast<double>
        
        //if(i<10) mexPrintf("coord: %f %f %f \n",x,y,z);

		vertex u(x, y, z);
        u.hops = i+1;
        u.color = 0;

    	VGraphAddAndGrow(vertices, (VNode)&u, i+1);	  
    }

    //
    // elements
    //

   
    primitives = VCreateGraph(colEle, 9, VLongRepn,  false);

    if (primitives == NULL)  {
      mexErrMsgTxt("cannot create primitives graph.\n");
    }

    mp = VCreateImage(1, 1, colEle, VUByteRepn);
    VImageNFrames( mp )     = 1;
    VImageNComponents( mp ) = 1;
    if ( mp == 0 ) {
      mexErrMsgTxt("generateCells: cannot create attribute images\n");
    }
    
    VAppendAttr(VImageAttrList(mp), "component_repn",   0, VStringRepn, "scalar" );
    VAppendAttr(VImageAttrList(mp), "component_interp", 0, VStringRepn, "element label");

    if(nargin == 5){
        AttributeImage = VCreateImage(6, 1, colEle, VDoubleRepn);
        if (AttributeImage == 0)
            mexErrMsgTxt("Allocate: cannot create attribute image\n");
    
        VSetAttr(VImageAttrList(AttributeImage),"component_repn",NULL,VStringRepn,"tensor6");
        VSetAttr(VImageAttrList(AttributeImage),"component_interp",NULL,VStringRepn,"conductivity S/mm");
    }
    
    primitive p;
    
	for(int i=0; i<colEle; i++) {

      // FE-type: m = 303 -> tet    m = 323 -> hex 
      // Surface type m = 302 -< triangle    m = 312 -> quatriliteral
    
      if (rowEle == 4) {

        p.vcnt  = 4;
        p.id[0] = (int)elements[i]; //static_cast<int>(elements[i]);
        p.id[1] = (int)elements[colEle + i];
        p.id[2] = (int)elements[2*colEle + i];
        p.id[3] = (int)elements[3*colEle + i];

	  linkNodes(vertices, p.id[0], p.id[1]);
	  linkNodes(vertices, p.id[0], p.id[2]);
	  linkNodes(vertices, p.id[0], p.id[3]);
	  linkNodes(vertices, p.id[1], p.id[2]);
	  linkNodes(vertices, p.id[1], p.id[3]);
	  linkNodes(vertices, p.id[2], p.id[3]);

        VGraphAddAndGrow(primitives, (VNode)&p, i+1);
      }
      
      else if (rowEle == 8){                     // internal connectivity: Hughes convention
//        fprintf(stdout," %i %i %i %i %i %i %i %i\n",id[0], id[1], id[2], id[3], id[4], id[5], id[6], id[7]); fflush(stdout);

        p.vcnt  = 8;
        p.id[0] = (int)elements[i]; //static_cast<int>(elements[i]);
        p.id[1] = (int)elements[colEle + i];
        p.id[2] = (int)elements[2*colEle + i];
        p.id[3] = (int)elements[3*colEle + i];
        p.id[4] = (int)elements[4*colEle + i];
        p.id[5] = (int)elements[5*colEle + i];
        p.id[6] = (int)elements[6*colEle + i];
        p.id[7] = (int)elements[7*colEle + i];


		linkNodes(vertices, p.id[0], p.id[1]);
		linkNodes(vertices, p.id[1], p.id[2]);
		linkNodes(vertices, p.id[2], p.id[3]);
		linkNodes(vertices, p.id[3], p.id[0]);
		linkNodes(vertices, p.id[0], p.id[4]);
		linkNodes(vertices, p.id[1], p.id[5]);
		linkNodes(vertices, p.id[2], p.id[6]);
		linkNodes(vertices, p.id[3], p.id[7]);
		linkNodes(vertices, p.id[4], p.id[5]);
		linkNodes(vertices, p.id[5], p.id[6]);
		linkNodes(vertices, p.id[6], p.id[7]);
		linkNodes(vertices, p.id[7], p.id[4]);

        VGraphAddAndGrow(primitives, (VNode)&p, i+1);
      }
      else if (rowEle == 3) {
        p.vcnt  = 3;
        p.id[0] = (int)elements[i];
        p.id[1] = (int)elements[colEle + i];
        p.id[2] = (int)elements[2*colEle + i];

	    linkNodes(vertices, p.id[0], p.id[1]);
	    linkNodes(vertices, p.id[0], p.id[2]);
	    linkNodes(vertices, p.id[1], p.id[2]);

		VGraphAddAndGrow(primitives, (VNode)&p, i+1);
      }
      
      //if(i<10) mexPrintf("primitive: %d %d %d %d %d %d %d %d \n",p.id[0],p.id[1],p.id[2],p.id[3],p.id[4],p.id[5],p.id[6],p.id[7]);
      
      if ((rowEle == 4)||(rowEle == 8)) {
          int m = (int)labels[i];
          //if(i<10) mexPrintf("label: %d \n",m);
          // set mp image
          VSetPixel(mp, 0, 0, i, m);
          
          if( nargin == 5){
            for(int j=0; j < 6; j++ ){		// for all dimensions, sequence xx-xy-xz-yy-yz-zz
                double tens = (double)tensors[j*colEle + i];
                VSetPixel(AttributeImage, j, 0, i, tens);
            }
        }
      }

    
    }

    VAppendAttr(VGraphAttrList(primitives), "matprops",   NULL, VImageRepn, mp);
    if(nargin == 5)
        VSetAttr(VGraphAttrList(primitives),"condtensor",NULL,VImageRepn,AttributeImage);
    
    VSetAttr(VGraphAttrList(primitives), "component_interp", NULL, VStringRepn, "primitive");
    if ( rowEle == 4 || rowEle == 8 )
        VSetAttr(VGraphAttrList(primitives), "primitive_interp", NULL, VStringRepn, "volume");
    else if ( rowEle == 3 )
        VSetAttr(VGraphAttrList(primitives), "primitive_interp", NULL, VStringRepn, "surface");
    VSetAttr(VGraphAttrList(primitives), "implicit_links", NULL, VStringRepn, "false");
    
    //trim graphs
    assert(vertices->lastUsed   == (int)countNodes(vertices));
    vertices->size     = vertices->lastUsed;
    vertices->nnodes   = vertices->lastUsed;
    
    assert(primitives->lastUsed   == (int)countNodes(primitives));
    primitives->size     = primitives->lastUsed;
    primitives->nnodes   = primitives->lastUsed;
        
    //save graphs
    VGraph g[2];
   
/*    vertices->size = colNod;
    vertices->nnodes = colNod;
    primitives->size = colEle;
    primitives->nnodes = colEle;*/
   
   
    g[0] = vertices;
    g[1] = primitives;
    VWriteGraphs(outf, NULL, 2, g);
    fclose(outf);
    if (vertices) 		VDestroyGraph(vertices);   vertices   = 0;
    if (primitives) 	VDestroyGraph(primitives); primitives = 0;
    
}
