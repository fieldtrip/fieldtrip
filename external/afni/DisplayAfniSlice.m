function [err] = DisplayAfniSlice (Slc, Opt)
%
%   [err] = DisplayAfniSlice (Slc, Opt)
%
%Purpose:
%
%    OBSOLETE,  see DispAFNISlice
%    the function is left here for backward compatibility
%
%   Display slices a la AFNI.
%
%
%Input Parameters:
%   Slc is a structure obtained from the output of function GetAfniSlice
%
%   Opt is a structure wth the following fields
%      .fhandle: figure handle, default is current figure
%      .subplot: (2x1 integer vector) number of images per figure
%      .colrange: (2x1 vector) the percentile of values to map to the color
%             maps. A la Afni, [0 0]  for min-max, [2 98] to match afni's window
%      .plane string containing the plane displayed (same as in GetAfniSlice).
%      .Info is the structure returned from BrikInfo
%      .iSlc is the vector containing the slice indices used by GetAfniSlice
%
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%
%
%
%Key Terms:
%
%More Info :
%   GetAfniSlice
%   Test_GetAfniSlice
%
%  DispAFNISlice, quite better.
%
%     Author : Ziad Saad
%     Date : Mon Aug 28 11:47:32 PDT 2000
%     LBC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'DisplayAfniSlice';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;


%if (~isfield(Opt,'') | isempty(Opt.)),	Opt. = ; end
if (~isfield(Opt,'fhandle') | isempty(Opt.fhandle)),	Opt.fhandle = gcf; end
if (~isfield(Opt,'colrange') | isempty(Opt.colrange)),	Opt.colrange = [0 0]; end
if (~isfield(Opt,'subplot') | isempty(Opt.subplot)),	Opt.subplot = [1 1]; end
if (~isfield(Opt,'plane') | isempty(Opt.plane)),	err = ErrEval(FuncName, 'Err_must specify plane'); return; end
if (~isfield(Opt,'Info') | isempty(Opt.Info)),	err = ErrEval(FuncName, 'Err_must specify BRIK Info'); return; end
if (~isfield(Opt,'iSlc') | isempty(Opt.iSlc)),	err = ErrEval(FuncName, 'Err_must specify slice numbers'); return; end

switch Opt.plane,
	case 'Ax',
		planename = 'Axial';
	case 'Sa',
		planename = 'Sagittal';
	case 'Co',
		planename = 'Coronal';
end

%number of slices
	szslc = size(Slc.M);
	prod(Opt.subplot);
	
%make sure enough subplots are there
if (prod(Opt.subplot) < szslc(3)),
	err= ErrEval(FuncName,'Err_Not enough subplots to display all slices');
	return;
end

% perform required scaling;	
%I'd like to normalize to cover the color map
%from 0 to the number of colors in the map
	map = colormap;
	nmap = max(size(map));
	lb = 0;
	nM = prod(szslc);
	
	if (Opt.colrange(1) | Opt.colrange(2)), %not all points are to be scaled
		M = reshape(Slc.M,[nM 1]);
		[Msort, iMsort] = sort(M);
		lbval = Msort(round(Opt.colrange(1).*nM./100));
		ubval = Msort(round(Opt.colrange(2).*nM./100));
		iToScale = find (Msort >= lbval & Msort <= ubval);
		
		Mmin = min ( Msort(iToScale) );
		Mmax = max ( Msort(iToScale) );
		if (Mmin == Mmax),
			Slc.M = ones(szslc).*nmap;
		else
			M(iMsort(1:min(iToScale)-1)) = lb;
			M(iMsort(iToScale)) = (((M(iMsort(iToScale)) - Mmin) ./ (Mmax - Mmin)) .* (nmap-lb)) + lb;
			M(iMsort(max(iToScale)+1:nM)) = nmap;
			Slc.M = reshape(M,szslc);
		end
	else
		Mmin = min (Slc.M(:));
		Mmax = max (Slc.M(:));
		if (Mmin == Mmax),
			Slc.M = ones(szslc).*nmap;
		else
			Slc.M = (((Slc.M - Mmin) ./ (Mmax - Mmin)) .* (nmap-lb)) + lb;
		end
	end

	
cf = figure(Opt.fhandle); clf
stit = sprintf ('%s--- %s', Opt.Info.RootName, planename);
set (cf,'Name',stit);

nslc = size(Slc.M,3);
for (i=1:1:nslc),
	subplot (Opt.subplot(1), Opt.subplot(2), i);
	image (Slc.M(:,:,i));
	axis equal;
	stit = sprintf('Slice %g', Opt.iSlc(i));
	title (stit);
	%axis ([0 szslc(1) 0 szslc(2)]);
end

err = 0;
return;

