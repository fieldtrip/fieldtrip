function [err, Itrans] = AFNI_IndexChange (Info, Iorig, Direction, DispOrient)
%
%   [err, Itrans] = AFNI_IndexChange (Info, Iorig, Direction, [DispOrient])
%
%Purpose:
%   Change the AFNI Index system between the display's and AFNI's coordinate system
%
%
%Input Parameters:
%   Info: The data structure output from BrikInfo
%   Iorig : anNx3 matrix containing the Ix Iy Iz indices of N points
%   Direction : A string : either 'A2D' or 'D2A' meaning Afni to Display or vice versa
%   DispOrient : The orientation string (or vector) for AFNI's display. This parameter is optional
%    and the defualt is 'RAI", it could be 'LAI' if AFNI's using Left is Left option.
%
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%   Itrans: if Iorig is specified, Itrans is Iorig in the new coordinate system
%
%
%
%Key Terms:
%
%More Info :
%   Test_AFNI_IndexChange
%
%   see also AFNI_CoordChange
%
%
%     Author : Ziad Saad
%     Date : Fri Sep 8 12:21:01 PDT 2000
%     LBC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'AFNI_IndexChange';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;
Itrans= [];
maplocation = [0 0 0];
mapsign = [0 0 0];


if (nargin == 3),
	DispOrient = [0 3 4]; %that's RAI
end


%Assume we're going from Display to Afni ('D2A)
if (ischar(DispOrient)),
	[err, OrCode] = AFNI_OrientCode (DispOrient);
else
	OrCode = DispOrient;
end

	[err, TrCode] = AFNI_OrientCode (Info.Orientation(:,1)');

%check if that's what is required
if (strmatch(Direction, 'A2D')),
	%we're going from Afni to Display coordinates
		tmp = TrCode;
		TrCode = OrCode;
		OrCode = tmp;
else
	if (~strmatch(Direction, 'D2A')),
		err = ErrEval(FuncName,'Err_Bad Direction string');
		return;
	end
end	
	

%get maplocation and mapsign, automatically from AFNI_CoordChange (no need to rewrite the code here)
	[err,maplocation, mapsign] = AFNI_CoordChange (OrCode, TrCode);

Itrans = Iorig;
if (strmatch(Direction, 'A2D')),
		fprintf ('Doing A2D\n');
		%(1-mapsign(i))./2 is 0 when mapsign(i) = 1 and 1 when mapsign(i) = -1;  This way, the if condition for mapsign(i) can be done away with
		%I left the simple method for the second loop for clarity. i don't think there's much efficiency difference between the two.
	for (i=1:1:3),
		Itrans(:,i) = ( (1-mapsign(i))./2 .* (Info.DATASET_DIMENSIONS(maplocation(i)) - 1) ) + ( mapsign(i) .* Iorig(:,maplocation(i)) );
	end
else
		fprintf ('Doing D2A: \n');
	for (i=1:1:3),
		Itrans(:,i) = Iorig(:,maplocation(i));
		if (mapsign(i) < 0),
			Itrans(:,i) = Info.DATASET_DIMENSIONS((i)) -1 - Itrans(:,i);
		end
	end	
end


err = 0;
return;

