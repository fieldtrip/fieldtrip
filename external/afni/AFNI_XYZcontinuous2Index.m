function [err,Indx] = AFNI_XYZcontinuous2Index (XYZmm, Info, CoordCode, IndexDim, retNaNoutside)
%
%   [err,Indx] = AFNI_XYZcontinuous2Index (XYZmm, Info, [CoordCode], [IndexDim], [retNaNoutside)
%
%Purpose:
%   Change from voxel XYZ in mm to XYZindex (called Voxel Coords in AFNI)
%   The mm and voxel coordinates refer to the values displayed
%   on the top left corner of AFNI controller.
%   CoordCode is the one you'd set from the Coord Order plugin
%
%
%Input Parameters:
%   XYZmm : The continuous coordinates corresponding to Indx
%       The coordnate system is assumed to be RAI (DICOM)
%       unless otherwise specified by CoordCode
%   Info is the output of BrikInfo
%   CoordCode is an optional parameter used to specify the
%            coordinates system of XYZmm. If empty or not specified,
%            the default is 'RAI'. The code can be either a string or a vector
%            of numbers (see AFNI_CoordChange for more on that)
%   IndexDim (3 or 1) is an optional parameter used to specify if Indx
%                     is Mx3 or Mx1 vector
%      (see AfniIndex2AfniXYZ for more info)
%   	The default is 3 . If you choose to specify IndexDim, you must
%     specify CoordCode
%      (you could use an empty string to leave CoordCode to the default)
%
%   [NNO Nov 2009 added <n.oosterhof@bangor.ac.uk>]
%
%   retNaNoutside (true or false) specifies if NaN should be returned if
%   a coordinate is outsize the volume specified by Info. The default is
%   false and yields valid indices in Indx of the nearest voxel. If set to
%   true, then this yields NaN values in Indx for coordinates outside the
%   volume
%
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%   Indx an Mx3 matrix or an  Mx1 vector (depending on IndexDim)
%        containing the voxel indices to be
%        transformed to voxel coordinates.  (indices start at 0)
%
%
%
%Key Terms:
%
%More Info :
%   BrikInfo
%   Test_AFNI_XYZcontinuous2Index
%   AFNI_Index2XYZcontinuous
%   Test_AFNI_Index2XYZcontinuous
%
%   and usage of  IJK_TO_DICOM_REAL illustrated in help for
%   function AFNI_Index2XYZcontinuous
%
%
%     Author : Ziad Saad
%     Date : Thu Sep 7 16:50:38 PDT 2000        Latest Modification: Feb 18 04
%     LBC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'AFNI_XYZcontinuous2Index';

%Debug Flag
DBG = 1;

ChangeCoord = 0;
if (nargin > 2)
	if (~isempty(CoordCode)),
		ChangeCoord = 1;
	end
end

ChangeDim = 0;
if (nargin >= 4), %NNO was '=='
	if (~isempty(IndexDim)),
		ChangeDim = 1;
	end
end


%initailize return variables
err = 1;
Indx = [];

% NNO added
if nargin < 5 %check last input argument
    retNaNoutside=false; % default: use old behaviour, do not return NaNs
elseif ~islogical(retNaNoutside) || numel(retNaNoutside) ~= 1
    fprintf(2,'%s: Illegal value for retNaNoutside\n\n', FuncName);
    return
end




Indx = XYZmm;

%make sure coordinate system is RAI
if (ChangeCoord),
	[err, maplocation, mapsign, XYZdic] = AFNI_CoordChange (CoordCode, 'RAI', XYZmm);
else
   XYZdic = XYZmm;
end

% now change dicomm to 3dmm
[err, XYZmm, map] = THD_dicomm_to_3dmm (Info, XYZdic);

	%The equations that would change the coordinate system to indices must take the indces in the same
	%RAI permutation that the slices are entered into to3d in (No need to worry about R versus L or A versus P)
	%determine the ordering map to go from any permutation of RAI to RAI
		%[maploc(1),jnk] = find(Info.Orientation == 'R');
		%[maploc(2),jnk] = find(Info.Orientation == 'A');
		%[maploc(3),jnk] = find(Info.Orientation == 'I');

		%pre - Wed May 23 18:20:56 PDT 2001 - WRONG !
		%Indx(:,1) = round( ( XYZmm(:, maploc(1)) - Info.ORIGIN(1) ) ./ Info.DELTA(1) );
		%Indx(:,2) = round( ( XYZmm(:, maploc(2)) - Info.ORIGIN(2) ) ./ Info.DELTA(2) );
		%Indx(:,3) = round( ( XYZmm(:, maploc(3)) - Info.ORIGIN(3) ) ./ Info.DELTA(3) );
		
		%post - Wed May 23 18:20:56 PDT 2001 - WRONG !
		%Indx(:,maploc(1)) = round( ( XYZmm(:, 1) - Info.ORIGIN(maploc(1)) ) ./ Info.DELTA(maploc(1)) );
		%Indx(:,maploc(2)) = round( ( XYZmm(:, 2) - Info.ORIGIN(maploc(2)) ) ./ Info.DELTA(maploc(2)) );
		%Indx(:,maploc(3)) = round( ( XYZmm(:, 3) - Info.ORIGIN(maploc(3)) ) ./ Info.DELTA(maploc(3)) );
		
      %Feb 18 04

    % NNO added
    IndxCount=size(Indx,1); % number of indices
    settoNaNmask=false(IndxCount,1); % values to be set to NaN
    outsidevalues=[0 0 0;Info.DATASET_DIMENSIONS(1:3)-1]; % if not set to NaN, values outside the volume.

    % for x, y, z dimension see which indices fall outside the volume
    for dim=1:3
        IndxDim=round( ( XYZmm(:, dim) - Info.ORIGIN(dim) ) ./ Info.DELTA(dim) ); % transform to voxel (ijk) space

        mask1=IndxDim<0; % values too low
        mask2=IndxDim>Info.DATASET_DIMENSIONS(dim); % values too high

        if retNaNoutside
            settoNaNmask=settoNaNmask | mask1 | mask2; % add to mask
        else
            IndxDim(mask1)=outsidevalues(1,dim); % } set outside
            IndxDim(mask2)=outsidevalues(2,dim); % } value
        end

        Indx(:,dim)=IndxDim;
    end

    %Now, if needed, change the Index dimension from Mx3 to Mx1
    if ChangeDim && IndexDim == 1
        Indx1=zeros(IndxCount,1);
        [err, Indx1(~settoNaNmask)] = AfniXYZ2AfniIndex (Indx(~settoNaNmask,:), Info.DATASET_DIMENSIONS(1), Info.DATASET_DIMENSIONS(2));
        Indx=Indx1;
    end

    % set to NaN if requested
    if retNaNoutside
        Indx(settoNaNmask,:)=NaN;
    end
	
err = 0;
return;

