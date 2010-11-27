function [M_V] = vox2ras_ksolve(M_R, Vc_Ps, varargin)
%%
%% NAME
%%
%%     vox2ras_ksolve.m (vox2ras_k{-space column}solve) 
%%
%% AUTHOR
%%
%%	Rudolph Pienaar
%%
%% VERSION
%%
%%	$Id$
%%
%% SYNOPSIS
%%
%%     [M_V] = vox2ras_ksolve(M_R, Vc_Ps, Vr_logicalSpace)
%%
%% ARGUMENTS
%%
%%	M_R		in	4x4 vox2ras candidate matrix with correct
%%					direction cosines, but invalid
%%					k-space column
%%      Vc_Ps		in      column vector defining the slice position
%%					(typically read from meas.asc)
%%	Vr_logicalSpace	in/opt	row vector that defines the logical dimension
%%					size. If omitted, this is assumed
%%					to be [256 256 128]
%%	M_V		out	complete 4x4 vox2ras matrix 
%%
%% DESCRIPTION
%%
%%	"vox2ras_ksolve" uses the given rotation components of M_R and the
%%	starting point, Vc_Ps, to determine the center of k-space and thus
%%	create a fully-qualified vox2ras matrix.
%%
%% PRECONDITIONS
%%
%%	o M_R is 4x4 where the first 3x3 submatrix contains the direction 
%%		cosines. The 4th column is ignored, as is the 4th row. 
%%	o The vector Vc_Ps is typically read from a Siemens meas.asc file such 
%%	  that
%%		Vc_Ps(1)	= sSliceArray.asSlice[0].sPosition.dSag
%%		Vc_Ps(2)	= sSliceArray.asSlice[0].sPosition.dCor
%%		Vc_ps(3)	= sSliceArray.asSlice[0].sPosition.dTra
%%
%% POSTCONDITIONS
%%
%%	o M_V(1:3, 1:3) = M_R(1:3, 1:3)
%%	o M_V(:,4)	= k-space center
%%
%% SEE ALSO
%%
%%	vox2ras_rsolveAA- determine the rotational component of a vox2ras matrix
%%				using Siemens reference orientations directly
%%	vox2ras_rsolve	- determine the rotational component of a vox2ras matrix
%%				using Siemens reference orientations indirectly
%%	vox2ras_dfmeas	- main function: determines the vox2ras matrix from a
%%			  Siemens meas.asc file.
%% 
%% HISTORY
%%
%% 26 May 2004
%% o Initial design and coding.
%%


%
% vox2ras_ksolve.m
%
% Original Author: Rudolph Pienaar
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:10 $
%    $Revision$
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

M_V		= M_R;
M_V(4,:)	= M_R(4,:)*0;
M_V(:,4)	= M_R(:,4)*0;
M_V(4,4)	= 1;

Vr_logicalSpace	= [ 256 256 128 ];
if length(varargin)
	Vr_logicalSpace = varargin{1};
end

%% First read the direction cosines
Vc_x 		= M_R(1:3,1);
Vc_y 		= M_R(1:3,2);
Vc_z 		= M_R(1:3,3);

xoff		= Vr_logicalSpace(1)/2;
yoff		= Vr_logicalSpace(2)/2;
zoff		= Vr_logicalSpace(3)/2;

%% Solve for the k-space center:
%% The first two components are found by moving along the direction cosines
Vc_Pe1		=  - (Vc_Ps + (xoff*Vc_x + yoff*Vc_y + zoff*Vc_z));
%% and the last component is found by moving against the direction cosines
Vc_Pe2		= (Vc_Ps - (xoff*Vc_x + yoff*Vc_y + zoff*Vc_z));
Vc_Pe(1)	= Vc_Pe1(1) -0.67;		%% a strange correction?
Vc_Pe(2)	= Vc_Pe1(2);
Vc_Pe(3)	= Vc_Pe2(3);
M_V(1:3, 4)	= Vc_Pe';

%% All done!

