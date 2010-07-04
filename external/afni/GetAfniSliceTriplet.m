function [err, slc] = GetAfniSliceTriplet (Brik, BrikInfo, DM, Opt)
%
%   [err, slc] = GetAfniSliceTriplet (Brik, BrikInfo, DM, Opt)
%
%Purpose:
%   Get the slices to display as they appear in the AFNI window
%   
%   
%Input Parameters:
%   Brik: is a vector or a matrix containing the brik values
%   BrikInfo is a structure containing Header info
%      (both are output by BrikLoad function)
%
%   DM is the output of function AFNI_SliceDispManip
%
%   Opt is a an options structure with the following fields
%    .iSlc is an 1x3 vector with the Axial, Sagittal and Coronal slice indices to show (a la AFNI, first slice is indexed 0)
%        for examle [45 12 3] retrieves Axial slice 45, Sagittal slice 12 and Coronal slice 3.
%        if you need just one slice, then pass -1 where you do not need any slices
%        for example [45 -1 3] retrieves the axial and coronal slices only
%   .index is the sub-brick index, default is 0
%
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%   slc is a structure [3x1] vector containing the slices requested in .iSlc
%			(where iSlc is -1, slc is empty)
%    .M   [N x M] the slice as present in the brick
%    .Mdisp [O x P] the slice ready for display
%   
%      
%Key Terms:
%   
%More Info :
%   Test_DispAFNISlice
%   GetSurfSliceTriplet
%
%
%     Author : Ziad Saad
%     Date : Tue Aug 22 18:36:04 PDT 2000
%     LBC/NIMH/ National Institutes of Health, Bethesda Maryland
%    last modified: Mon Aug 27 13:29:42 PDT 2001


%Define the function name for easy referencing
FuncName = 'GetAfniSliceTriplet';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;

%turn Brick into a matrix if it isn't already
if (ndims(Brik) < 3),
	Brik = reshape(Brik, BrikInfo.DATASET_DIMENSIONS(1), BrikInfo.DATASET_DIMENSIONS(2),...
						BrikInfo.DATASET_DIMENSIONS(3), BrikInfo.DATASET_RANK(2));
end

if (~isfield(Opt, 'index') | isempty(Opt.index)), 
	Opt.index = 0;
end

for (ip = 1:1:3),
	if (Opt.iSlc(ip) > 0),
		
		if (DM(ip).zflip), 
			%added -1 on Mon Aug 27 13:29:42 PDT 2001
			Opt.iSlc(ip) = BrikInfo.DATASET_DIMENSIONS(DM(ip).SliceDim) - Opt.iSlc(ip) -1;
		end

		%grab the slices
		switch DM(ip).SliceDim,
			case 1, 				
				slc(ip).M = permute(Brik(Opt.iSlc(ip)+1, :, : , Opt.index+1), [2 3 1]);
			case 2
				slc(ip).M = permute(Brik(:, Opt.iSlc(ip)+1, :, Opt.index+1), [1 3 2]); 
			case 3
				slc(ip).M = Brik(:, :, Opt.iSlc(ip)+1, Opt.index+1);
		end 
		%now manipulate them to have them display properly 
			if (DM(ip).orpermute), slc(ip).Mdisp = permute(slc(ip).M,[2 1]);	
				else slc(ip).Mdisp = slc(ip).M; end
			if (DM(ip).udflip), slc(ip).Mdisp = flipud(slc(ip).Mdisp);	end
			if (DM(ip).lrflip), slc(ip).Mdisp = fliplr(slc(ip).Mdisp);	end
	else
		slc(ip).M = [];
		slc(ip).Mdisp = [];
	end
end %plane

err = 0;
return;

