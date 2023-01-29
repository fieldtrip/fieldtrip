function [err,PatchHandles] = ShowPlane (EqPlane,Opt)
%
%   [err,PatchHandles] = ShowPlane (EqPlane,Opt)
%
%Purpose:
%   This function displays a plane with equation ax + by + cz +d = 0
%   The plane is defined by one sqaure patch
%   The planes span the current axis limits of the figure.
%
%Input Parameters:
%   EqPlane is a Nx4 matrix containing the equation of the plane containing each
%       triplet in  Triplets. The plane passing by triplet i is speicifed in
%       EqPlane(i,:) the plane would be
%       EqPlane(i,1)x + EqPlane(i,2)y + EqPlane(i,3)z + EqPlane(i,4) = 0
%   Opt is the options structure
%     .Fig is a handle to the figure you want the planes displayed in,
%         default is the current figure.
%     .WriteIV if this string is not empty, the planes that are displayed
%         on the graph are written to an inventor format file
%     .units 'mm' or 'tesscon'. default is tesscon
%       if you specify mm, then the xyz coordinates are transformed to tesscon
%       before writing them out. (*319.7). The idea is to write all inventor
%       files in tesscon units. This option will only be used if WriteIV is not empty.
%     .OvrWrite (0/1) default is 0, flag for overwriting existing .iv file
%
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%
%   PatchHandles : the handle to the patches displayed on the figure
%
%More Info :
%   see Plane_Equation
%
%
%
%     Author : Ziad Saad
%     Date : Thu Oct 22 20:19:36 CDT 1998


%Define the function name for easy referencing
FuncName = 'ShowPlane';

%initailize return variables
err = 1;


%check on the size of input data
if (nargin == 1),	
	Opt.Fig = [];	
	Opt.OvrWrite = 0;
	Opt.WriteIV = '';
end

if (~isfield(Opt,'OvrWrite') | isempty(Opt.OvrWrite)), Opt.OvrWrite = 0; end

if (size(EqPlane,2) ~= 4),	err = ErrEval(FuncName,'Err_Bad size for EqPlane');	return;	end

Nplanes = size(EqPlane,1);
Nnodes = 4.* Nplanes;

%pop up a figure
if (~isfield(Opt,'Fig') | isempty(Opt.Fig)),
	Opt.Fig = gcf;
end

figure(Opt.Fig);

%get the axis roperties of the figure
Xlim = get(gca,'Xlim');
Ylim = get(gca,'ylim');
Zlim = get(gca,'Zlim');

inode = 0;
Node = zeros(Nnodes,3);
Pat = zeros(Nplanes,4);
ztmp = zeros(1,4);

for (i=1:1:Nplanes),
	%using the XY limits, find the corrspondign z values
	if (EqPlane(i,3) ~= 0),
		ztmp(1) = (-EqPlane(i,4) - EqPlane(i,1).*Xlim(1) - EqPlane(i,2).*Ylim(1)) ./ EqPlane(i,3);
		ztmp(2) = (-EqPlane(i,4) - EqPlane(i,1).*Xlim(2) - EqPlane(i,2).*Ylim(1)) ./ EqPlane(i,3);
		ztmp(3) = (-EqPlane(i,4) - EqPlane(i,1).*Xlim(2) - EqPlane(i,2).*Ylim(2)) ./ EqPlane(i,3);
		ztmp(4) = (-EqPlane(i,4) - EqPlane(i,1).*Xlim(1) - EqPlane(i,2).*Ylim(2)) ./ EqPlane(i,3);
		%form the four points on the plane
		Node(inode+1,:) = [Xlim(1) Ylim(1) ztmp(1)];
		Node(inode+2,:) = [Xlim(2) Ylim(1) ztmp(2)];
		Node(inode+3,:) = [Xlim(2) Ylim(2) ztmp(3)];
		Node(inode+4,:) = [Xlim(1) Ylim(2) ztmp(4)];
		
	elseif (EqPlane(i,2) ~= 0),
		ytmp(1) = (-EqPlane(i,4) - EqPlane(i,1).*Xlim(1) - EqPlane(i,3).*Zlim(1)) ./ EqPlane(i,2);
		ytmp(2) = (-EqPlane(i,4) - EqPlane(i,1).*Xlim(2) - EqPlane(i,3).*Zlim(1)) ./ EqPlane(i,2);
		ytmp(3) = (-EqPlane(i,4) - EqPlane(i,1).*Xlim(2) - EqPlane(i,3).*Zlim(2)) ./ EqPlane(i,2);
		ytmp(4) = (-EqPlane(i,4) - EqPlane(i,1).*Xlim(1) - EqPlane(i,3).*Zlim(2)) ./ EqPlane(i,2);
		%form the four points on the plane
		Node(inode+1,:) = [Xlim(1) ytmp(1) Zlim(1)];
		Node(inode+2,:) = [Xlim(2) ytmp(2) Zlim(1)];
		Node(inode+3,:) = [Xlim(2) ytmp(3) Zlim(2)];
		Node(inode+4,:) = [Xlim(1) ytmp(4) Zlim(2)];
	elseif (EqPlane(i,1) ~= 0),
		xtmp(1) = (-EqPlane(i,4) - EqPlane(i,2).*Ylim(1) - EqPlane(i,3).*Zlim(1)) ./ EqPlane(i,1);
		xtmp(2) = (-EqPlane(i,4) - EqPlane(i,2).*Ylim(2) - EqPlane(i,3).*Zlim(1)) ./ EqPlane(i,1);
		xtmp(3) = (-EqPlane(i,4) - EqPlane(i,2).*Ylim(2) - EqPlane(i,3).*Zlim(2)) ./ EqPlane(i,1);
		xtmp(4) = (-EqPlane(i,4) - EqPlane(i,2).*Ylim(1) - EqPlane(i,3).*Zlim(2)) ./ EqPlane(i,1);
		%form the four points on the plane
		Node(inode+1,:) = [xtmp(1) Ylim(1) Zlim(1)];
		Node(inode+2,:) = [xtmp(2) Ylim(2) Zlim(1)];
		Node(inode+3,:) = [xtmp(3) Ylim(2) Zlim(2)];
		Node(inode+4,:) = [xtmp(4) Ylim(1) Zlim(2)];
	end	

	%verify that all points are on plane for debugging only
	%sum(EqPlane.*[Node(inode+1,:) 1])
	%sum(EqPlane.*[Node(inode+2,:) 1])
	%sum(EqPlane.*[Node(inode+3,:) 1])
	%sum(EqPlane.*[Node(inode+4,:) 1])
	
	%form the faceset connection, for the patch
	Pat(i,:) = [inode+1 inode+2 inode+3 inode+4];

	inode = inode + 4;
	
end

%Now display those patches
vc = 0:1./Nnodes:(1-1./Nnodes);
tcolor = [vc' (1-vc)' vc'];

PatchHandles = patch('vertices',Node,'faces',Pat,...
          'FaceVertexCData',tcolor,'FaceColor','flat');

%write to file ?
if (isfield(Opt,'WriteIV') & ~isempty(Opt.WriteIV)),
	fprintf (1,'Saving patches to iv file %s ...\n',Opt.WriteIV);
   Opt.OptIV.BaseCol = tcolor;
   if (~isfield(Opt,'units') | isempty(Opt.units)),	Opt.units = 'tesscon'; end
	Opt.OptIV.units = Opt.units;
	Opt.OptIV.OvrWrite = Opt.OvrWrite;
	Opt.OptIV.verbose = 0;
   [err] = WriteInv21Surf(Opt.WriteIV,Node,Pat,Opt.OptIV);  %exclude redundant last node
end

err = 0;
return;

