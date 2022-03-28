function [err,M] = MakeColorMap (Fiducials,Ncols,Opt)
%
%   [err,M] = MakeColorMap (Fiducials,Ncols,Opt)
%
%Purpose:
%   returns the RGB colormap containing Ncols that vary linearily
%   from the first color in Fiducials to the last.
%
%Input Parameters:
%   Fiducials : Nx3 matix specifying the RGB values (0-255) or
%              (0-1) depending on the value of Opt.Range. Those
%              fiducial colours will be equally spaced on the map
%   Ncols : Total number of colours in the map
%           You are somewhat restricted in the total number of
%           colours you choose. You must choose a number that
%           allows you to have the same number of colours between
%           each successive fiducials. Do not worry, the function
%           will suggest a good number closest to the one you chose.
%				Such a mighty function !
%   Opt is a structure containing the following fields
%   	.Range (225 / 1) This specifies if RGB values are specified
%           in both Fiducials and M form 0-255 (integers)
%           or 0-1 floats.
%     .SkipLast (0/1) if set to 0, then the last color specified in
%          Fiducials is the last color in M. If set to 1, the last
%          color in M represents the color that would come right
%          before the last one in Fifucials. This last option is
%          usefull when you're crating cyclical color maps where
%          the last color in Fiduciasl is like the first.
%     .Showme (0/1) optional parameter to show a display of the map
%     .Write optional string. If supplied, M is written to the file
%       specified by .Write
%
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%     .verbose (0/1) optional verbose parameter, default is 0
%More Info :
%   example
%      Fiducials = [255 0 0; 0 255 0; 0 0 255];
%      Opt.Range = 255; Opt.SkipLast = 1; Opt.Write = '';
%      [err,M] = MakeColorMap (Fiducials,6,Opt)
% gives M =
%   255     0     0
%   128   128     0
%     0   255     0
%     0   128   128
%     0     0   255
%   128     0   128
%
%If you set Opt.SkipLast = 0, you'll get
%     [err,M] = MakeColorMap (Fiducials,7,Opt);
%M =
%
%   255     0     0
%   128   128     0
%     0   255     0
%     0   128   128
%     0     0   255
%   128     0   128
%   255     0     0
%
%   see also rgbdectohex, and ScaleToMap
%            ROIcmap, ShowCmap, readXcol
%
%     Author : Ziad Saad
%     Date : Wed Apr 08 12:51:29 CDT 1998


%Define the function name for easy referencing
FuncName = 'MakeColorMap';

%initailize return variables
err = 1;
M = [];

if (~isfield(Opt,'Range') | ~isfield(Opt,'SkipLast')),
	fprintf(2, 'Error %s: Range or SkipLast fields missing.\n', FuncName);
   err = 1; return;
end

if (~isfield(Opt,'verbose')),
	Opt.verbose = 0;
end

if (~isfield(Opt,'Write')),
	Opt.Write = '';
else
	if (filexist(Opt.Write)),
		err = ErrEval(FuncName,'Err_FileExist');	return
	end
end

if (~isfield(Opt,'Showme')),
	Opt.Showme = 1;
end

if (size(Fiducials,2) ~= 3),
	err = ErrEval(FuncName,'Err_BadOptSize');	return
end

if (size(Fiducials,1) > Ncols),
	err = ErrEval(FuncName,'Err_More fiducials than colours needed.');	return
end

%Check for some weird input
tmp =  find (Fiducials < 0);
if (~isempty(tmp)),
	err= ErrEval(FuncName,'Err_No negative values allowed in Fiducials.');	return
end

tmp =  find (Fiducials > Opt.Range);
if (~isempty(tmp)),
	serr = sprintf('Err_Values in Fiducials should not exceed Opt.Range (%g).',Opt.Range);
	err= ErrEval(FuncName,serr);	return
end


%initialize
M = -ones(Ncols,3);
Nfid = size (Fiducials,1); %number of fiducial colours

if (~Opt.SkipLast),
	Ninter = Ncols - Nfid; %total number of intermediate colours
else
	Ninter = Ncols - (Nfid -1);
end

Ngap = Nfid - 1;	%total number of gaps to fill

%You must have an equal number of colours in each gap
Npergap = Ninter ./ Ngap;

if (Npergap ~= round(Npergap)),
	if (Opt.SkipLast) Ncolsgood = round(Npergap) .* Ngap + Nfid -1;	
		else	Ncolsgood = round(Npergap) .* Ngap + Nfid;	end
	serr = sprintf('Err_The choice of Ncols does not work with the number\nof fiducial colours.\nTry Ncols = %g ',Ncolsgood);
	M = [];
   err = ErrEval(FuncName,serr);	return;
end

%Start forming M
	cnt = 0;
	im = 1;
	for (i=1:1:Ngap),
		if (Fiducials(i,1)~=Fiducials(i+1,1)),
			M1 = linspace(Fiducials(i,1),Fiducials(i+1,1),Npergap+2)';
		else
			M1 = Fiducials(i,1) .* ones(Npergap+2,1);	
		end
		
		if (Fiducials(i,2)~=Fiducials(i+1,2)),
			M2 = linspace(Fiducials(i,2),Fiducials(i+1,2),Npergap+2)';
		else
			M2 = Fiducials(i,2) .* ones(Npergap+2,1);	
		end
		
		if (Fiducials(i,3)~=Fiducials(i+1,3)),
			M3 = linspace(Fiducials(i,3),Fiducials(i+1,3),Npergap+2)';
		else
			M3 = Fiducials(i,3) .* ones(Npergap+2,1);	
		end
		
		
		im2 = im+Npergap+1;
		if (i<Ngap | ~Opt.SkipLast),
			M(im:im2,:) = [M1 M2 M3];
		else
			M(im:im2-1,:) = [M1(1:Npergap+1) M2(1:Npergap+1) M3(1:Npergap+1)];
		end
		im = im2;
	
	end	


%make sure format for output is good
if (Opt.Range == 255),	%needs to be integer output
	M = round(M);
end

if (Opt.Showme),
	figure;
	if (Opt.Range == 255),
		Mrgb = M ./ 255;
	else
		Mrgb = M;
	end
	colormap (Mrgb);
	subplot 211;
	image ([1:1:length(Mrgb(:,1))]);
	subplot 212;
	pie (ones(1,length(Mrgb(:,1))));
end

if (~isempty(Opt.Write)),
	wryte2(M,3,Opt.Write,'on');
end

err = 0;
return;

