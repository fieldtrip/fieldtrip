/*

 REPOP.C
 Generalized arithmetic operators which implicity replicate their inputs when
 they are not large enough for the input dimensions.

 This file contains the generic type-independent c-code utility functions
 which can be used stand-alone (i.e. w/o reference to matlab)
 
 This should be compilied and linked with with repop.c.

 Copyright 2006-     by Jason Farquhar (jdrf@zepler.org)
 Permission is granted for anyone to copy, use, or modify this
 software and accompanying documents for any uncommercial purposes,
 provided this copyright notice is retained, and note is made of
 any changes that have been made. This software and documents are
 distributed without any warranty, express or implied.

 This code is based upon ideas by Douglas M. Schwarz (schwarz@servtech.com)
 and Adi Vehtari -- but heavily modifed to make it more general, easier to
 understand and faster!

*/
#include "string.h"
#include "mxInfo.h"
#include "repop.h"

/*---------------------------------------------------------------------------*/

/* compute the size of the output z result & id dims to replicate    */
MxInfo initzinfo(const MxInfo xinfo, const MxInfo yinfo, int repNonUnitDim){
  int i;
  int znd=xinfo.nd;
  MxInfo zinfo=mkemptymxInfo(znd);
  zinfo.stride[0]=1; 
  zinfo.numel    =1;

  /* set the result type to the least precise of the input types */
  if ( xinfo.dtype == DOUBLE_DTYPE && yinfo.dtype == DOUBLE_DTYPE ) {
	 zinfo.dtype    = DOUBLE_DTYPE ;
  } else if ( xinfo.dtype == SINGLE_DTYPE || yinfo.dtype == SINGLE_DTYPE ){
	 zinfo.dtype    = SINGLE_DTYPE;
  } else if ( xinfo.dtype == INT32_DTYPE && yinfo.dtype == INT32_DTYPE ) {
	 zinfo.dtype    = INT32_DTYPE ;
  } else { /* set to invalid type */
	 zinfo.dtype = 0 ; 
  }

  for( i=0; i<znd; i++ ){
	 /* z size is max of x or y size */
	 zinfo.sz[i]       = max(xinfo.sz[i],yinfo.sz[i]);
	 zinfo.stride[i+1] = zinfo.stride[i]*zinfo.sz[i];
	 zinfo.numel      *= zinfo.sz[i];
	 
	 /* check this type of replication is allowed */
	 if( xinfo.sz[i]>1 && yinfo.sz[i]>1 && xinfo.sz[i]!=yinfo.sz[i]) {
		if ( repNonUnitDim==0 && 
			  ((xinfo.sz[i]>yinfo.sz[i])?yinfo.sz[i]>1:xinfo.sz[i]>1) )
		  ERROR("REPOP:Replicating a non-unit dimension! -- if you *really* meant this use the 'm' or 'n' options");
		if ( repNonUnitDim>=1 && /* warning if not integer multiples */
			  ((xinfo.sz[i]>yinfo.sz[i])?xinfo.sz[i]%yinfo.sz[i]!=0:yinfo.sz[i]%xinfo.sz[i]!=0) ) {
		  if ( repNonUnitDim==1 ) /* non-integer is an error */
			 ERROR("REPOP:X and Y size not integer multiples of each other! -- if you *really* meant this use the 'n' option");
		}
	 }
  }
  return zinfo;
}


/*-------------------------------------------------------------------*/
/* compress output size to minimise the dimension loop overhead */
void repopqueryOptimise(MxInfo *zinfo, MxInfo *xinfo, MxInfo *yinfo, 
								int *repxy){
  int znd=zinfo->nd;
  int i,resnd=0;
  for( i=1; i<znd; i++ ){
	 if ( repxy[i-1]==0 && xinfo->stride[i-1] > 0 && yinfo->stride[i-1] > 0 && 
			repxy[i]==0 && xinfo->stride[i] > 0 && yinfo->stride[i] > 0 ) { 
		/* non-rep pair dims -- COMPRESS*/
		zinfo->sz[resnd]=zinfo->sz[resnd]*zinfo->sz[i];
		xinfo->sz[resnd]=xinfo->sz[resnd]*xinfo->sz[i];
		yinfo->sz[resnd]=yinfo->sz[resnd]*yinfo->sz[i];

	 } else if ( repxy[i-1]==0 && repxy[i]==0 &&
					 ( (xinfo->stride[i-1]==0 && xinfo->stride[i]==0) 
						|| (yinfo->stride[i-1]==0 && yinfo->stride[i]==0) ) ){
		/* unit sized rep in X or Y -- COMPRESS */
		zinfo->sz[resnd]=zinfo->sz[resnd]*zinfo->sz[i];
		xinfo->sz[resnd]=xinfo->sz[resnd]*xinfo->sz[i];
		yinfo->sz[resnd]=yinfo->sz[resnd]*yinfo->sz[i];		

	 } else {/* normally replicated dim -- can't compress so move on as normal*/
		resnd++;
		zinfo->sz[resnd]=zinfo->sz[i]; zinfo->stride[resnd]=zinfo->stride[i];
		xinfo->sz[resnd]=xinfo->sz[i]; xinfo->stride[resnd]=xinfo->stride[i];
		yinfo->sz[resnd]=yinfo->sz[i]; yinfo->stride[resnd]=yinfo->stride[i];
		repxy[resnd]=repxy[i];
	 }
  }
  resnd++;
  zinfo->nd = xinfo->nd = yinfo->nd = resnd;
  zinfo->stride[resnd]=zinfo->stride[i];
  xinfo->stride[resnd]=xinfo->stride[i];
  yinfo->stride[resnd]=yinfo->stride[i];
}

/* char to int lookup of the operator name */
int getOpid(const char *str){
  if ( strcmp(str,"+")==0  || strcmp(str,".+")==0 || 
		 strcmp(str,"plus")==0 || strcmp(str,"PLUS")==0 )
	 return PLUS;
  if ( strcmp(str,"-")==0  || strcmp(str,".-")==0 || 
		 strcmp(str,"minus")==0 || strcmp(str,"MINUS")==0)
	 return MINUS;
  if ( strcmp(str,"*")==0  || strcmp(str,".*")==0 || 
		 strcmp(str,"times")==0 || strcmp(str,"TIMES")==0 )
	 return TIMES;
  if ( strcmp(str,"/")==0  || strcmp(str,"./")==0 || 
		 strcmp(str,"rdivide")==0 || strcmp(str,"RDIVIDE")==0)
	 return RDIVIDE;
  if ( strcmp(str,"\\")==0 || strcmp(str,".\\")==0|| 
		 strcmp(str,"ldivide")==0 || strcmp(str,"LDIVIDE")==0 )
	 return LDIVIDE;
  if ( strcmp(str,"^")==0  || strcmp(str,".^")==0 || 
		 strcmp(str,"power")==0 || strcmp(str,"POWER")==0)
	 return POWER;
  if ( strcmp(str,"%")==0  ||strcmp(str,"mod")==0 || strcmp(str,"MOD")==0 )
	 return MOD;
  if ( strcmp(str,"==")==0 || strcmp(str,"eq")==0 || strcmp(str,"EQ")==0 )
	 return EQ;
  if ( strcmp(str,"~=")==0 || strcmp(str,"ne")==0 || strcmp(str,"NE")==0 )
	 return NE;
  if ( strcmp(str,"<")==0  || strcmp(str,"lt")==0 || strcmp(str,"LT")==0 )
	 return LT;
  if ( strcmp(str,">")==0  || strcmp(str,"gt")==0 || strcmp(str,"GT")==0 )
	 return GT;
  if ( strcmp(str,"<=")==0 || strcmp(str,"le")==0 || strcmp(str,"LE")==0 )
	 return LE;
  if ( strcmp(str,">=")==0 || strcmp(str,"ge")==0 || strcmp(str,"GE")==0 )
	 return GE;
  if ( strcmp(str,"min")==0 || strcmp(str,"MIN")==0 )
	 return MINOP;
  if ( strcmp(str,"max")==0 || strcmp(str,"MAX")==0 ) 
	 return MAXOP;
  return -1;
}


/******************************************************************************
*   Tests to see if every imaginary element is identically zero and, if so,   *
*   removes the complex part.                                                 *
******************************************************************************/
char removeZeroImag(MxInfo *zinf) {
        
  /* check if we have any non-zero imaginary coefficients */
  /* only necessary for double/single types */
  if ( zinf->dtype == DOUBLE_DTYPE ) {
	 double *endp, *ip;
	 for(ip=zinf->ip,endp=ip+zinf->numel; ip<endp; ip++) if( *ip!=0.0 ) break;
	 if ( ip!=endp ) return 1; /* if any complex coeff do nothing */
  } else if( zinf->dtype == SINGLE_DTYPE ) {
	 float *endp, *ip;
	 for(ip=(float*)zinf->ip,endp=ip+zinf->numel; ip<endp; ip++) if( *ip!=0.0 ) break;
	 if ( ip!=endp ) return 1; /* if any complex coeff do nothing */  
  } else {
	 /* Not a valid complex type so do nothing */
	 return 0 ; 
  }

  /* if we got here then we have complex pointer but no complex info */
  /*simpler and faster version is to simply set the arrays pi pointer to null*/
  FREE(zinf->ip);  /* free the old memory and return */
  zinf->ip=0; 
  return 0;
}
