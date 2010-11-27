function [err,maplocation, mapsign, Mtrans] = AFNI_CoordChange (Orig, Trans, Morig)
%
%   [err, maplocation, mapsign, Mtrans] = AFNI_CoordChange (Orig, Trans, [Morig])
%
%Purpose:
%   Change the AFNI coordinate system between the 48 different coordinate systems allowed
%   Coordinates are as they appear on the top left corner of AFNI's controller
%   Coord Order plugin is used to change the coordinate system in AFNI
%   
%Input Parameters:
%   Orig : the 1x3 letter or number code for the original orientation, like RAI or [1 3 4]
%   Trans : the 1x3 letter or number code for the new orientation
%   Morig : an (optional) Nx3 matrix containing the XYZ coordinates of N points
%   
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%   maplocation: 1x3 vector containig the map from the old coordinate to the new coordinate system 
%         This specifies where each dimension in Mtrans is located in Morig
%   mapsign: 1x3 vector (of 1 or -1) containing the sign of the map from the old coordinate system to the new one 
%        This specifies if the dimension has a negative direction (see more Info for an example)
%   Mtrans: if Morig is specified, Mtrans is Morig in the new coordinate system
%   
%      
%Key Terms:
%   
%More Info :
%
%   maplocation and mapsign are used as such
%
%	for (i=1:1:3),
%		Mtrans(:,i) = mapsign(i).* Morig(:,maplocation(i));
%	end
%   
% example:  
%[err,maplocation, mapsign, Mtrans] = AFNI_CoordChange ('RAI', 'ALI', [1 2 3; 4 5 6])
%
%     Author : Ziad Saad
%     Date : Tue Sep 5 18:56:40 PDT 2000
%     LBC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'AFNI_CoordChange';

%Debug Flag
DBG = 1;

if (nargin == 2),
	Morig = [];
end

if (~isempty(Morig) & size(Morig,2) ~= 3),
	err = ErrEval(FuncName,'Err_Bad size for Morig');
	return;
end


%initailize return variables
err = 1;
Mtrans= [];
maplocation = [0 0 0];
mapsign = [0 0 0];

if (ischar(Orig)),
	[err, OrCode] = AFNI_OrientCode (Orig);
else
	OrCode = Orig;
end

if (ischar(Trans)),
	[err, TrCode] = AFNI_OrientCode (Trans);
else
	TrCode = Trans;
end

for (i=1:1:3),
	%look for the orientation
	itmp = find (OrCode == TrCode(i));
	if (~isempty(itmp)), %found, no need for flipping
		maplocation(i) = itmp;
		mapsign(i) = 1;
	else  %look for opposite orientation, need flipping
		if (rem(TrCode(i),2)), 
			shft = -1;
		else
			shft = 1;
		end
		itmp = find(OrCode == TrCode(i)+shft);
		if (isempty(itmp)),
			err = ErrEval(FuncName,'Err_Bad code duuude');
		end
		maplocation(i) = itmp;
		mapsign(i) = -1;
	end
end %for i

if (~isempty(Morig)),
	Mtrans = Morig;
	for (i=1:1:3),
		Mtrans(:,i) = mapsign(i).* Morig(:,maplocation(i));
	end
end

err = 0;
return;

