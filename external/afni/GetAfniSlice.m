function [err, slc, slcDisp] = GetAfniSlice (Brik, BrikInfo, Opt)
%
%   [err, slc, slcDisp] = GetAfniSlice (Brik, BrikInfo, Opt)
%
%Purpose:
%   Display a slice as it appears in the AFNI window
%
%
%Input Parameters:
%   Brik: is a vector or a matrix containing the brik values
%   BrikInfo is a structure containing Header info
%      (both are output by BrikLoad function)
%   Opt is a an options structure with the following fields
%    .plane is the slicing plane. Choose from 'Ax', 'Sa' and 'Co' for axial, sagittal and coronal
%    .iSlc is an Nx1 vector with slice indices to show (a la AFNI, first slice is indexed 0)
%    .index is the sub-brick index, default is 0
%
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%   slc is a structure containing the slices
%    .M   [N x M x k] , k being the number of slices requested in Opt.i*
%    .orientax a 2x2 character matrix. The 1st and 2nd rows indicates the orientation
%          of the voxels in the 1st  and 2nd dimension of .M, respectively.
%          for example: if orientax is [RL; AP] or
%          RL
%          AP
%          This means that the pixels in slc.M go from Right to Left along the 1st dimension of the
%          matrix (ie along columns, top to bottom). Pixels go from Anterior to Posterior along the second
%          dimension of the matrix (ie along rows, left to right).
%
%   slcDisp is a structure similar to slc and containing the slices oriented as they would be displayed in AFNI's window
%
%
%Key Terms:
%
%More Info :
%   Test_GetAfniSliceTriplet
%   GetAfniSlice
%	 AFNI_SliceDispManip
%	GetAfniSliceTriplet
%
%
%     Author : Ziad Saad
%     Date : Tue Aug 22 18:36:04 PDT 2000
%     LBC/NIMH/ National Institutes of Health, Bethesda Maryland
%    last modified: Mon Aug 27 13:29:42 PDT 2001


%Define the function name for easy referencing
FuncName = 'GetAfniSlice';

%Debug Flag
DBG = 1;

%turn Brick into a matrix if it isn't already
if (ndims(Brik) < 3),
	Brik = reshape(Brik, BrikInfo.DATASET_DIMENSIONS(1), BrikInfo.DATASET_DIMENSIONS(2),...
						BrikInfo.DATASET_DIMENSIONS(3), BrikInfo.DATASET_RANK(2));
end

%initailize return variables
err = 1;
	switch Opt.plane,
		case 'Ax',
			iplane = 1;
		case 'Sa',
			iplane = 2;
		case 'Co',
			iplane = 3;
	end

if (~isfield(Opt,'index') |isempty(Opt.index)), Opt.index = 0; end

indx = Opt.index+1;

DispManip = AFNI_SliceDispManip (BrikInfo);		

if (DispManip(iplane).zflip),
			%added -1 on Mon Aug 27 13:29:42 PDT 2001
	Opt.iSlc = BrikInfo.DATASET_DIMENSIONS(DispManip(iplane).SliceDim) - Opt.iSlc -1;
end

	nz = length(Opt.iSlc);
	if (nz),
		%allocate enough space
			slc.M = zeros(BrikInfo.DATASET_DIMENSIONS(DispManip(iplane).d(1)), BrikInfo.DATASET_DIMENSIONS(DispManip(iplane).d(2)), nz);
			bufslc = zeros(BrikInfo.DATASET_DIMENSIONS(DispManip(iplane).d(1)), BrikInfo.DATASET_DIMENSIONS(DispManip(iplane).d(2)));
			if (DispManip(iplane).orpermute),
				slcDisp.M = zeros(BrikInfo.DATASET_DIMENSIONS(DispManip(iplane).d(2)), BrikInfo.DATASET_DIMENSIONS(DispManip(iplane).d(1)), nz);
			else
				slcDisp.M = zeros(BrikInfo.DATASET_DIMENSIONS(DispManip(iplane).d(1)), BrikInfo.DATASET_DIMENSIONS(DispManip(iplane).d(2)), nz);
			end
			
		%first grab the slices
		switch DispManip(iplane).SliceDim,
			case 1, 				
				for (j=1:1:nz),
					slc.M(:,:,j) = Brik(Opt.iSlc(j)+1, :, :, indx);
				end
			case 2
				for (j=1:1:nz),
					slc.M(:,:,j) = Brik(:, Opt.iSlc(j)+1, :, indx);
				end
			case 3
				for (j=1:1:nz),
					slc.M(:,:,j) = Brik(:, :, Opt.iSlc(j)+1, indx);
				end
		end %switch DispManip(iplane).SliceDim,
		
		%now manipulate them to have them display properly
		for (j=1:1:nz),
			bufslc = slc.M(:,:,j); %this is inefficient, might want to dump this step
			if (DispManip(iplane).orpermute), bufslc = permute(bufslc,[2 1]);	end
			if (DispManip(iplane).udflip), bufslc = flipud(bufslc);	end
			if (DispManip(iplane).lrflip), bufslc = fliplr(bufslc);	end
			slcDisp.M(:,:,j) = bufslc;
		end
	end %if (nz)


err = 0;
return;

