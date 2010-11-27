function [err,Eq] = Plane_Equation (Triplets, verbose)
%
%   [err,Eq] = Plane_Equation (Triplets, [verbose])
%
%Purpose:
%   Determine the equation of the plane passing through three points
%   
%   
%Input Parameters:
%   Triplets is an Nx1 vector of strucutres. Each structure defines a plane
%    .XYZ : is a 3x3 matrix containing the XYZ of each of the three points
%           Each row is a point. ie: XYZ = [0 0 0; 1 0 0; 1 1 1]; is for the
%           three points (0 0 0), (1 0 0) and (1 1 1).
%           If the three points are colinear, Eq = [0 0 0 0]
%
%   verbose (0/1), default is 1
%   
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%   
%   Eq is a Nx4 matrix containing the equation of the plane containing each
%       triplet in  Triplets. The plane passing by triplet i is speicifed in
%       Eq(i,:) the plane would be Eq(i,1)x + Eq(i,2)y + Eq(i,3)z + Eq(i,4) = 0
%      
%More Info :
%   
%   see also ShowPlane
%   try 
%       Triplets(1).XYZ = [0 0 0; 1 0 0; 1 1 0];
%       Triplets(2).XYZ = [0 0 0; 1 0 0; 1 1 1];
%       Triplets(3).XYZ = [0 5 0; 1 5 0; 1 1 1];
%
%       [err,Eq] = Plane_Equation (Triplets);
%       [err,PatchHandles] = ShowPlane (Eq); view(3)
%
%     Author : Ziad Saad
%     Date : Thu Oct 22 16:09:56 CDT 1998 


%Define the function name for easy referencing
FuncName = 'Plane_Equation';

%initailize return variables
err = 1;

if (nargin == 1), 	verbose = 1;	end

if (is_row(Triplets) == -1),	err = ErrEval(FuncName,'Err_Triplets must be an Nx1 vector');	return;	end

Triplets = Triplets(:);
Nplanes = size(Triplets,1);
Eq = zeros(Nplanes,4); %allocate

for (i=1:1:Nplanes),
	%Form the equation of the plane
	x1 = Triplets(i).XYZ(1,1); y1 = Triplets(i).XYZ(1,2); z1 = Triplets(i).XYZ(1,3);
	x2 = Triplets(i).XYZ(2,1); y2 = Triplets(i).XYZ(2,2); z2 = Triplets(i).XYZ(2,3);
	x3 = Triplets(i).XYZ(3,1); y3 = Triplets(i).XYZ(3,2); z3 = Triplets(i).XYZ(3,3);
	Eq(i,1) = y1.*(z2-z3) + y2.*(z3-z1) + y3.*(z1-z2);
	Eq(i,2) = z1.*(x2-x3) + z2.*(x3-x1) + z3.*(x1-x2);
	Eq(i,3) = x1.*(y2-y3) + x2.*(y3-y1) + x3.*(y1-y2);
	Eq(i,4) = -x1.*(y2.*z3 - y3.*z2) - x2.*(y3.*z1 - y1.*z3) - x3.*(y1.*z2 - y2.*z1);
	
	if (verbose),
		if (~rem (i,200)),
			fprintf (1,'[%g/%g]\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b',i,Nplanes);
		end

	end
end

err = 0;
return;

