function [err,XYZdic] = AFNI_Index2XYZcontinuous (Indx, Info, CoordCode)
%
%   [err,XYZdic] = AFNI_Index2XYZcontinuous (Indx, Info, [CoordCode])
%
%Purpose:
%   Change from voxel XYZindex (called Voxel Coords in AFNI) to XYZ in mm 
%   The mm and voxel coordinates refer to the values displayed 
%   on the top left corner of AFNI controller.
%   CoordCode is the one you'd set from the Coord Order plugin
%   
%   
%Input Parameters:
%   Indx an Mx3 matrix or an  Mx1 vector containing the voxel indices to be
%        transformed to voxel coordinates.  (indices start at 0)
%   Info is the output of BrikInfo
%   CoordCode is an optional parameter used to specify the coordinates system of the output
%      if empty or not specified, the default is 'RAI'. The code can be either a string or a vector 
%      of numbers (see AFNI_CoordChange for more on that)
%
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%   XYZdic : The continuous coordinates corresponding to Indx
%       The coordnate system output is in RAI (DICOM) 
%       unless otherwise specified by CoordCode
%   
%      
%Key Terms:
%   
%More Info :
%   BrikInfo
%   Test_AFNI_Index2XYZcontinuous
%   AFNI_XYZcontinuous2Index
%   Test_AFNI_XYZcontinuous2Index
%
%     Author : Ziad Saad
%     Date : Tue Sep 5 21:48:06 PDT 2000           Latest Modification: Feb 18 04
%     LBC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'AFNI_Index2XYZcontinuous';

%Debug Flag
DBG = 1;

ChangeCoord = 0;
if (nargin == 3)
	if (~isempty(CoordCode)),
		ChangeCoord = 1;
	end   
end


%initailize return variables
err = 1;
XYZmm = [];

%make sure Indx is the right size 
switch size(Indx,2),
	case 1, %change 1D index to XYZ index
		[err, Indx] = AfniIndex2AfniXYZ (Indx, Info.DATASET_DIMENSIONS(1), Info.DATASET_DIMENSIONS(2))
	case 3, %OK
	otherwise,
		err = ErrEval(FuncName,'Err_Bad dimension for Indx');
		return
end

XYZmm = Indx;

	%The equations that would change the indices to coordinate system result in a coordinate system that 
	% may be any permutation of RAI (like IRA or AIR or IAR or RIA or ARI) so one only needs to find the 
	%dimension permutation needed to bring the final result to RAI.

	%determine the ordering map to go from any permutation of RAI to RAI
		%[maploc(1),jnk] = find(Info.Orientation == 'R');
		%[maploc(2),jnk] = find(Info.Orientation == 'A');
		%[maploc(3),jnk] = find(Info.Orientation == 'I');
	
	%pre - Wed May 23 18:20:56 PDT 2001 - WRONG !
		%XYZmm(:, maploc(1)) = Info.ORIGIN(1) + Indx(:,1) .* Info.DELTA(1);
		%XYZmm(:, maploc(2)) = Info.ORIGIN(2) + Indx(:,2) .* Info.DELTA(2);
		%XYZmm(:, maploc(3)) = Info.ORIGIN(3) + Indx(:,3) .* Info.DELTA(3);

	%post - Wed May 23 18:20:56 PDT 2001 - WRONG! 
		%XYZmm(:, 1) = Info.ORIGIN(maploc(1)) + Indx(:,maploc(1)) .* Info.DELTA(maploc(1));
		%XYZmm(:, 2) = Info.ORIGIN(maploc(2)) + Indx(:,maploc(2)) .* Info.DELTA(maploc(2));
		%XYZmm(:, 3) = Info.ORIGIN(maploc(3)) + Indx(:,maploc(3)) .* Info.DELTA(maploc(3));
	
   %Feb 18 04, back to the original
      XYZmm(:, 1) = Info.ORIGIN(1) + Indx(:,1) .* Info.DELTA(1);
		XYZmm(:, 2) = Info.ORIGIN(2) + Indx(:,2) .* Info.DELTA(2);
		XYZmm(:, 3) = Info.ORIGIN(3) + Indx(:,3) .* Info.DELTA(3);
      %Now this is in the axis orientation which is Info.Orientation(:,1)' called 3dmm in thd_coords.c
      [err,XYZdic, map] = THD_3dmm_to_dicomm (Info, XYZmm);

if (ChangeCoord),
	[err, maplocation, mapsign, XYZdic] = AFNI_CoordChange ('RAI', CoordCode, XYZdic);
end

err = 0;
return;

