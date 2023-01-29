function [DispManip] = AFNI_SliceDispManip (B_Info)
%
%   [DispManip] = AFNI_SliceDispManip (B_Info)
%
%Purpose:
%
%   Determine the manipulations necessary to apply to a
%  slice for proper display a la AFNI
%
%Input Parameters:
%   B_Info the output of BrickInfo function
%
%
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%   DispManip is a 3x1 structure with the following structures
%      .plane : 2 character string specifying the plane ('Ax', 'Sa', Co')
%            this is a bit redundant since the three elements of
%             DispManip always refer to 'Ax', 'Sa', 'Co', respectively
%   .SliceDim : The dimension in Brik that the slices are stored in
%      .zflip : must flip in the z-dimension
%
% The following flags are used for manipulating the slice once extracted
% from the brik. Think of them as being applied to .orient matrix to produce
% .Disporient matrix which is what AFNI displays
%      .orpermute
%      .udflip
%      .lrflip
%
%      .orient : a 2x2 character matrix. The 1st and 2nd rows indicates
%          the orientation  of the voxels in the 1st  and 2nd dimension
%          of the slice as stored in the Brik, respectively.
%      .Disporient : this represents the slice orientations once the
%          manipulatons are performed on the slices (or .orient)
%
%Key Terms:
%
%More Info :
%   see GetAfniSliceTriplet
%   and GetSurfSliceTriplet
%
%
%     Author : Ziad Saad
%     Date : Wed Sep 20 15:53:29 PDT 2000
%     LBC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'AFNI_SliceDispManip';

%Debug Flag
DBG = 1;

for (i = 1:1:3),
	switch i,
		case 1,
			DispManip(i).plane = 'Ax';
			[DispManip(i).SliceDim,SliceDim_col] = find(B_Info.Orientation == 'I'); %The first slice should be the Inferior one
		case 2,
			DispManip(i).plane = 'Sa';
			[DispManip(i).SliceDim,SliceDim_col] = find(B_Info.Orientation == 'L'); %The first slice should be the left one
		case 3,
			DispManip(i).plane = 'Co';		
			[DispManip(i).SliceDim,SliceDim_col] = find(B_Info.Orientation == 'A'); %The first slice should be the Anterior one	
	end
	
	%determine if flipping in Z-dimension is necessary
	if (SliceDim_col == 2), DispManip(i).zflip = 1;
		else DispManip(i).zflip = 0;
	end
	
	switch DispManip(i).SliceDim, %determine the Brik dimensions that contain the slice
		case 1,	% slices are in the first dimension
			DispManip(i).d = [2 3];
		case 2, %slices are in the second dimension
			DispManip(i).d = [1 3];
		case 3,	
			DispManip(i).d = [1 2];
	end
	
	%determine orientation in slice plane
		DispManip(i).orient = B_Info.Orientation(DispManip(i).d,:);
	
	%determine what manipulatons are required to have the slices display as in AFNI
	DispManip(i).udflip = 0;
	DispManip(i).lrflip = 0;
	DispManip(i).orpermute = 0;
	switch i, %different planes
		case 1,
			[iA,jA] = find (DispManip(i).orient == 'A'); [iL,jL] = find (DispManip(i).orient == 'R');
			if (iA == 2),  %must permute orders, A--P must be along the columns),
				DispManip(i).orpermute = 1;		end
			if (jA == 2), %must flip up to down, the column dimention should go in A--P not P--A
				DispManip(i).udflip = 1; end
			if (jL == 2), %must go L--R (radiologic format)
				DispManip(i).lrflip = 1; end
		case 2,
			[iA,jA] = find (DispManip(i).orient == 'S'); [iL,jL] = find (DispManip(i).orient == 'A');
			if (iA == 2),  %must permute orders, S--I must be along the columns),
				DispManip(i).orpermute = 1;		end
			if (jA == 2), %must flip up to down, the column dimention should go in S--I not I--S
				DispManip(i).udflip = 1; end
			if (jL == 2), %must go A--S not S--A
				DispManip(i).lrflip = 1; end
		case 3,
			[iA,jA] = find (DispManip(i).orient == 'S'); [iL,jL] = find (DispManip(i).orient == 'R');
			if (iA == 2),  %must permute orders, S--I must be along the columns),
				DispManip(i).orpermute = 1;		end
			if (jA == 2), %must flip up to down, the column dimention should go in S--I not I--S
				DispManip(i).udflip = 1; end
			if (jL == 2), %must go L--R (radiologic format)
				DispManip(i).lrflip = 1; end	
	end

	%apply the manipulations to the slice orientation matrix, this is good for debugging
		DispManip(i).Disporient = DispManip(i).orient;
		if (DispManip(i).orpermute), DispManip(i).Disporient = permute(DispManip(i).orient,[2 1]); end
		if (DispManip(i).udflip), DispManip(i).Disporient = flipud(DispManip(i).Disporient); end
		if (DispManip(i).lrflip), DispManip(i).Disporient = fliplr(DispManip(i).Disporient); end

end



err = 0;
return;

