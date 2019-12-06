function [M_V] = vox2ras_dfmeas(str_filename, varargin)
%%
%% NAME
%%
%%     vox2ras_dfmeas.m (vox2ras_d{etermine}f{rom}meas)
%%
%% AUTHOR 
%%
%%	Rudolph Pienaar
%%
%% VERSION
%%
%%	$Id: vox2ras_dfmeas.m,v 1.10 2011/03/02 00:04:13 nicks Exp $
%%
%% SYNOPSIS
%%
%%      [M_V] = vox2ras_dfmeas(str_filename, ch_override, ...
%%				voxelDimension spaceSize)
%%
%% ARGUMENTS
%%
%%	str_filename	in      string specifying meas.asc to parse
%%	ch_override	in/opt	character string specifying override behaviour
%%	voxelDimension	in/opt	row vector specifying the data set's voxel mm
%%					dimensions. This should be a 1x3 vector
%%					of [ReadOut PhaseEncode SliceSelect] 
%%	spaceSize	in/opt	row vector specifying the logical space
%%					dimensions. This should be a 1x3 vector
%%					of [rows columns slices]
%%	M_V		out	vox2ras matrix
%%
%% DESCRIPTION
%%
%%	"vox2ras_dfmeas"  examines a Siemens meas.asc file and returns a 
%%	vox2ras matrix, M_V.
%%
%%	There are several possible mechanisms by which M_V can be determined:
%%
%%	The most natural way is to simply parse M_V from the meas.asc file
%%	itself. Since the meas.asc does not record complete vox2ras information
%%	this method depends on "enhancements" to raw data saving sequences
%%	and may not be available on all meas.asc files. In the NMR center, all
%%	raw data saving sequences from 2004 onwards are "enhanced" and record
%%	complete vox2ras information.
%%
%%	In the absence of complete vox2ras information in the meas.asc file, we
%%	can attempt a vox2ras solution using only one direction cosine. By 
%%	default, the "standard" meas.asc file records a single direction cosine
%%	(normal to the coronal plane). Additional direction cosines for the 
%%	remaining dimensions can be calculated, using the vox2ras_rsolve(...) 
%%	function (see additional documentation for more details on this method). 
%%	This calculation makes several assumptions on the structure of the final
%%	vox2ras matrix in order to address numerous ambiguities inherent in 
%%	trying to solve for a 3D space given only one direction cosine. This 
%%	calculation may fail in some cases, or simply return an incorrect 
%%	vox2ras matrix. Visual sanity-check of the volume associated with the 
%%	given vox2ras matrix is highly recommended.
%%
%%	Finally, M_V can be simply assigned a "default" value and returned back 
%%	out of this function without any further processing. This is useful if 
%%	a suitable vox2ras was not found using a direct calculation.
%%
%%	In all cases, the center of k-space is also calculated and returned as 
%%	the 4th column in M_V.
%%
%%	By default, this function will
%%
%%		o First, parse for an embedded vox2ras in the given meas.asc
%%			If found, use this as the first 3x3 submatrix of
%%			vox2ras, and calculate the position of the k-space
%%			center.
%%
%%		o Second, if no embedded vox2ras is found, attempt to solve for
%%			the first 3x3 submatrix of vox2ras using the existing 
%%			direction cosine in meas.asc, and calculate the 
%%			position of the k-space center.
%%
%%		o Thirdly, if specified by ch_override, simply return a "default"
%%			first 3x3 submatrix and calculate the position of the
%%			k-space center.
%%
%%	ch_override values:
%%		o 'n'		No override - used if additional variable 
%%					arguments are pending.
%%		o 'c'		Force a vox2ras calculation even if an embedded 
%%					matrix was found.
%%		o 'd'		Force a "default" return (i.e. no parsing or 
%%					calculation for the first 3x3 submatrix. 
%%					Still solve for center of k-space).
%%		
%% PRECONDITIONS
%%
%%	o Siemens meas.asc scanner header file
%%	o The parsing of meas.asc is performed by dropping to the OS and doing 
%%	  some quick text parsing. This assumes a *nix environment!
%%	o The default calculations assume voxel dimensions of 
%%	  [1.0 1.0 1.33], i.e. standard GLEEK-type structure.
%%	o The default calculation also assume a logical space size of
%%	  [256 256 128] = [rows cols slices].
%%
%% POSTCONDITIONS
%%
%%	o a vox2ras matrix, M_V, is returned.
%%
%% SEE ALSO
%%
%%	vox2ras_rsolveAA- determine the rotational component of a vox2ras matrix
%%				using Siemens reference orientations directly
%%	vox2ras_rsolve	- determine the rotational component of a vox2ras matrix
%%				using Siemens reference orientations indirectly
%%	vox2ras_ksolve	- determine the k-space col in RAS of a vox2ras matrix
%%
%% HISTORY
%%
%% 25 May 2004
%% o Initial design and coding.
%%
%% 02 June 2004
%% o vox2ras_rsolveAA
%%


%
% vox2ras_dfmeas.m
%
% Original Author: Rudolph Pienaar
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:13 $
%    $Revision: 1.10 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

%% Check for override specifications and set default values
ch_override		= '0';
Vr_voxelDimension 	= [ 1 1 1.33 ];
Vr_logicalSpace		= [ 256 256 128 ];
if length(varargin)
	ch_override		= varargin{1}(1);
	Vr_voxelDimension = [ 1 1 1.33];
end
if length(varargin) == 2
	Vr_voxelDimension = varargin{2};
end
if length(varargin) == 3
	Vr_voxelDimension = varargin{2};
	Vr_logicalSpace	 = varargin{3};
end


M_Vdefault = [
 0	 		 0			-Vr_voxelDimension(3)	0
-Vr_voxelDimension(1)	 0	 		 0			0
 0			-Vr_voxelDimension(2)	 0			0
 0	 		 0	 		 0			1
];
M_voxScale 	= eye(4);
M_voxScale(1,1)	= Vr_voxelDimension(1);
M_voxScale(2,2)	= Vr_voxelDimension(2);
M_voxScale(3,3)	= Vr_voxelDimension(3);
M_W = 0;

%% Check if the file exists and is readable
fid = fopen(str_filename, 'r', 'l');
if fid < 0
    fprintf(1, 'Error - could not open file %s for reading.\n', str_filename);
    fprintf(1, 'Returning default vox2ras with no zero k-space center.\n');
    M_V = M_Vdefault;
    return;
end


%% Parse the meas.asc file for an embedded vox2ras matrix.
%%	This parsing relies on an underlying *nix environment. The
%% 	somewhat convoluted construction is necessary to embed the 
%%	single quote character in the commmand string
cmd 	= ['cat ' str_filename ' | grep -a adRM | awk ' char(39) ...
		'{print $4 " " $7 " " $10}' char(39)];
[s, v2r]= system(cmd);	%% Execute the command, capturing vox2ras in v2r
[r c]	= size(v2r);	
if r~=0
	%% If the meas.asc acutally had an embedded vox2ras, do some
	%%	fine tuning.
	M_mask1	= [
			 1	0	 0	0
			 0	1	 0	0
			 0	0	-1	0
			 0	0	 0	1
	];
	M_mask2	= [
		 	-1	0	 0	0
			 0	1	 0	0
			 0	0	-1	0
			 0	0	 0	1
	];
	M_Wn 		= str2num(v2r);
	M_W  		= eye(4);
	M_W(1:3, 1:3)	= M_Wn;	
	M_W  		= M_mask1 * M_W' * M_mask2 * M_voxScale;
end

%% Parse meas.asc file for sNormal direction cosine
Vc_N = zeros(3,1);
for i=1:3,
    switch i
	case 1
	    str_search = 'dSag';
	case 2
	    str_search = 'dCor';
	case 3
	    str_search = 'dTra';
    end   
    cmd = ['cat ' str_filename ' | grep -a sSliceArray.asSlice | grep sNormal.' ... 
    		str_search ' | awk ' char(39) '{print $3}' char(39)];
    [s, f]   = system(cmd);
    Vc_N(i)  = str2num(f);
end
%% Swap the direction sense of the first two normals (LPS to RAS)
Vc_N(1)	= -Vc_N(1);	
Vc_N(2)	= -Vc_N(2);
%% and scale with the Slice Select voxel dimension
Vc_N	= Vr_voxelDimension(3) * Vc_N;
%% Parse meas.asc file for InPlaneRotation value
cmd = ['cat ' str_filename ' | grep -a sSliceArray.asSlice | grep dInPlaneRot' ... 
    		' | awk ' char(39) '{print $3}' char(39)];
[s, inPlaneRot]    		= system(cmd);
inPlaneRotation			= str2num(inPlaneRot);

[M_R]				= vox2ras_rsolveAA(Vc_N, inPlaneRotation);

%% Parse meas.asc file for sPosition vector
Vc_Ps = zeros(3,1);
for i=1:3,
    switch i
	case 1
	    str_search = 'dSag';
	case 2
	    str_search = 'dCor';
	case 3
	    str_search = 'dTra';
    end   
    cmd = ['cat ' str_filename ' | grep -a sSliceArray.asSlice | grep sPosition.' ... 
    		str_search ' | awk ' char(39) '{print $3}' char(39)];
    [s, f]    = system(cmd);
    Vc_Ps(i)  = str2num(f);
end

%% Depending on file contents and/or override flags, fix the vox2ras matrix
%%	to process further.
if r==0
	M_V	= M_R;
else
	M_V	= M_W;
end

if ch_override == 'c'
	M_V	= M_R;
elseif ch_override == 'd'
	M_V	= M_Vdefault;
end

%% Finally, calculate the center of k-space
M_V	= vox2ras_ksolve(M_V, Vc_Ps, Vr_logicalSpace);

%% All done!
