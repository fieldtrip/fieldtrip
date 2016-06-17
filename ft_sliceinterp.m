function [outim] = ft_sliceinterp(cfg, ininterp)

% FT_SLICEINTERP plots a 2D-montage of source reconstruction and anatomical MRI
% after these have been interpolated onto the same grid.
%
% Use as
%   ft_sliceinterp(cfg, interp)
%      or
%   [rgbimage] = ft_sliceinterp(cfg, interp), rgbimage is the monatage image
%
% where interp is the output of sourceinterpolate and cfg is a structure
% with any of the following fields:
%
% cfg.funparameter  string with the functional parameter of interest (default = 'source')
% cfg.maskparameter parameter used as opacity mask (default = 'none')
% cfg.clipmin       value or 'auto' (clipping of source data)
% cfg.clipmax       value or 'auto' (clipping of source data)
% cfg.clipsym       'yes' or 'no' (default) symmetrical clipping
% cfg.colormap      colormap for source overlay (default is jet(128))
% cfg.colmin        source value mapped to the lowest color (default = 'auto')
% cfg.colmax        source value mapped to the highest color (default = 'auto')
% cfg.maskclipmin   value or 'auto' (clipping of mask data)
% cfg.maskclipmax   value or 'auto' (clipping of mask data)
% cfg.maskclipsym   'yes' or 'no' (default) symmetrical clipping
% cfg.maskmap       opacitymap for source overlay
%                   (default is linspace(0,1,128))
% cfg.maskcolmin    mask value mapped to the lowest opacity, i.e.
%                   completely transparent (default ='auto')
% cfg.maskcolmin    mask value mapped to the highest opacity, i.e.
%                   non-transparent (default = 'auto')
% cfg.alpha         value between 0 and 1 or 'adaptive' (default)
% cfg.nslices       integer value, default is 20
% cfg.dim           integer value, default is 3 (dimension to slice)
% cfg.spacemin      'auto' (default) or integer (first slice position)
% cfg.spacemax      'auto' (default) or integer (last slice position)
% cfg.resample      integer value, default is 1 (for resolution reduction)
% cfg.rotate        number of ccw 90 deg slice rotations (default = 0)
% cfg.title         optional title (default is '')
% cfg.whitebg       'yes' or 'no' (default = 'yes')
% cfg.flipdim       flip data along the sliced dimension, 'yes' or 'no'
%                   (default = 'no')
% cfg.marker        [Nx3] array defining N marker positions to display
% cfg.markersize    radius of markers (default = 5);
% cfg.markercolor   [1x3] marker color in RGB (default = [1 1 1], i.e. white)
% cfg.interactive   'yes' or 'no' (default), interactive coordinates
%                   and source values
%
% if cfg.alpha is set to 'adaptive' the opacity of the source overlay
% linearly follows the source value: maxima are opaque and minima are
% transparent.
%
% if cfg.spacemin and/or cfg.spacemax are set to 'auto' the sliced
% space is automatically restricted to the evaluated source-space
%
% if cfg.colmin and/or cfg.colmax are set to 'auto' the colormap is mapped
% to source values the following way: if source values are either all
% positive or all negative the colormap is mapped to from
% min(source) to max(source). If source values are negative and positive
% the colormap is symmetrical mapped around 0 from -max(abs(source)) to
% +max(abs(source)).
%
% If cfg.maskparameter specifies a parameter to be used as an opacity mask
% cfg.alpha is not used. Instead the mask values are maped to an opacitymap
% that can be specified using cfg.maskmap. The mapping onto that
% opacitymap is controlled as for the functional data using the
% corresponding clipping and min/max options.
%
% if cfg.whitebg is set to 'yes' the function estimates the head volume and
% displays a white background outside the head, which can save a lot of black
% printer toner.
%
% if cfg.interactive is set to 'yes' a button will be displayed for
% interactive data evaluation and coordinate reading. After clicking the
% button named 'coords' you can click on any position in the slice montage.
% After clicking these coordinates and their source value are displayed in
% a text box below the button. The coordinates correspond to indeces in the
% input data array:
%
%   f = interp.source(coord_1,coord_2,coord_3)
%
% The coordinates are not affected by any transformations used for displaying
% the data such as cfg.dim, cfg.rotate,cfg.flipdim or cfg.resample.
%
% See also FT_SOURCEANALYSIS, FT_VOLUMERESLICE

% Copyright (C) 2004, Markus Siegel, markus.siegel@fcdonders.kun.nl
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar ininterp
ft_preamble provenance ininterp
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
ininterp = ft_checkdata(ininterp, 'datatype', 'volume', 'feedback', 'yes');

% set the defaults
if ~isfield(cfg, 'clipmin');      cfg.clipmin = 'auto';        end
if ~isfield(cfg, 'clipmax');      cfg.clipmax = 'auto';        end
if ~isfield(cfg, 'clipsym');      cfg.clipsym = 'no';          end
if ~isfield(cfg, 'alpha');        cfg.alpha   = 'adaptive';    end
if ~isfield(cfg, 'nslices');      cfg.nslices = 20;            end
if ~isfield(cfg, 'dim');          cfg.dim = 3;                 end
if ~isfield(cfg, 'colormap');     cfg.colormap = jet(128);     end
if ~isfield(cfg, 'spacemin');     cfg.spacemin = 'auto';       end
if ~isfield(cfg, 'spacemax');     cfg.spacemax = 'auto';       end
if ~isfield(cfg, 'colmin');       cfg.colmin = 'auto';         end
if ~isfield(cfg, 'colmax');       cfg.colmax = 'auto';         end
if ~isfield(cfg, 'resample');     cfg.resample = 1;            end
if ~isfield(cfg, 'rotate');       cfg.rotate = 0;              end
if ~isfield(cfg, 'title');        cfg.title = '';              end
if ~isfield(cfg, 'whitebg');      cfg.whitebg = 'no';          end
if ~isfield(cfg, 'flipdim');      cfg.flipdim = 'no';          end
if ~isfield(cfg, 'marker');       cfg.marker = [];             end
if ~isfield(cfg, 'markersize');   cfg.markersize = 5;          end
if ~isfield(cfg, 'markercolor');  cfg.markercolor = [1,1,1];   end
if ~isfield(cfg, 'interactive');  cfg.interactive = 'no';      end
if ~isfield(cfg, 'maskclipmin');  cfg.maskclipmin = 'auto';           end
if ~isfield(cfg, 'maskclipmax');  cfg.maskclipmax = 'auto';           end
if ~isfield(cfg, 'maskclipsym');  cfg.maskclipsym = 'no';             end
if ~isfield(cfg, 'maskmap');      cfg.maskmap = linspace(0,1,128);    end
if ~isfield(cfg, 'maskcolmin');   cfg.maskcolmin = 'auto';            end
if ~isfield(cfg, 'maskcolmax');   cfg.maskcolmax = 'auto';            end
if ~isfield(cfg, 'maskparameter');cfg.maskparameter = [];             end

% perform some checks on the configuration for backward compatibility
if ~isfield(cfg, 'funparameter') && isfield(ininterp, 'source')
  % if present, the default behavior should be to use this field for plotting
  cfg.funparameter = 'source';
end

% make the selection of functional and mask data consistent with the data
cfg.funparameter  = parameterselection(cfg.funparameter, ininterp);
cfg.maskparameter = parameterselection(cfg.maskparameter, ininterp);
% only a single parameter should be selected
try, cfg.funparameter  = cfg.funparameter{1}; end
try, cfg.maskparameter = cfg.maskparameter{1}; end

% check anatomical data
if isfield(ininterp,'anatomy');
  interp.anatomy = reshape(ininterp.anatomy, ininterp.dim);
else
  error('no anatomical data supplied');
end

% check functional data
if ~isempty(cfg.funparameter)
  interp.source = double(reshape(getsubfield(ininterp, cfg.funparameter), ininterp.dim));
else
  error('no functional data supplied');
end

% check mask data
if ~isempty(cfg.maskparameter)
  interp.mask = double(reshape(getsubfield(ininterp,cfg.maskparameter), ininterp.dim));
  maskdat = 1;
else
  fprintf('no opacity mask data supplied\n');
  interp.mask = [];
  maskdat = 0;
end

% only work with the copy of the relevant parameters in "interp"
clear ininterp;

% convert anatomy data type and optimize contrast
if isa(interp.anatomy, 'uint8') || isa(interp.anatomy, 'uint16')
  fprintf('converting anatomy to floating point values...');
  interp.anatomy = double(interp.anatomy);
  fprintf('done\n');
end
fprintf('optimizing contrast of anatomical data ...');
minana = min(interp.anatomy(:));
maxana = max(interp.anatomy(:));
interp.anatomy = (interp.anatomy-minana)./(maxana-minana);
fprintf('done\n');

% store original data if 'interactive' mode
if strcmp(cfg.interactive,'yes')
  data.source = interp.source;
end

% place markers
marker = zeros(size(interp.anatomy));
if ~isempty(cfg.marker)
  fprintf('placing markers ...');
  [x,y,z] = ndgrid([1:size(interp.anatomy,1)],[1:size(interp.anatomy,2)],[1:size(interp.anatomy,3)]);
  for imarker = 1:size(cfg.marker,1)
    marker(find(sqrt((x-cfg.marker(iarker,1)).^2 + (y-cfg.marker(imarker,2)).^2 + (z-cfg.marker(imarker,3)).^2)<=cfg.markersize)) = 1;
  end
  fprintf('done\n');
end

% shift dimensions
fprintf('sorting dimensions...');
interp.anatomy = shiftdim(interp.anatomy,cfg.dim-1);
interp.source = shiftdim(interp.source,cfg.dim-1);
interp.mask = shiftdim(interp.mask,cfg.dim-1);
marker = shiftdim(marker,cfg.dim-1);
fprintf('done\n');

% flip dimensions
if strcmp(cfg.flipdim,'yes')
  fprintf('flipping dimensions...');
  interp.anatomy = flipdim(interp.anatomy,1);
  interp.source = flipdim(interp.source,1);
  interp.mask = flipdim(interp.mask,1);
  marker = flipdim(marker,1);
  fprintf('done\n');
end

% set slice space
if ischar(cfg.spacemin)
  fprintf('setting first slice position...');
  spacemin = min(find(~isnan(max(max(interp.source,[],3),[],2))));
  fprintf('%d...done\n',spacemin);
else
  spacemin = cfg.spacemin;
end

if ischar(cfg.spacemax)
  fprintf('setting last slice position...');
  spacemax = max(find(~isnan(max(max(interp.source,[],3),[],2))));
  fprintf('%d...done\n',spacemax);
else
  spacemax = cfg.spacemax;
end

% clip funtional data
if ~ischar(cfg.clipmin)
  fprintf('clipping functional minimum...');
  switch cfg.clipsym
  case 'no'
    interp.source(find(interp.source<cfg.clipmin)) = nan;
  case 'yes'
    interp.source(find(abs(interp.source)<cfg.clipmin)) = nan;
  end
  fprintf('done\n');
end
if ~ischar(cfg.clipmax)
  fprintf('clipping functional maximum...');
  switch cfg.clipsym
  case 'no'
    interp.source(find(interp.source>cfg.clipmax)) = nan;
  case 'yes'
    interp.source(find(abs(interp.source)>cfg.clipmax)) = nan;
  end
  fprintf('done\n');
end

% clip mask data
if maskdat
  if ~ischar(cfg.maskclipmin)
    fprintf('clipping mask minimum...');
    switch cfg.maskclipsym
    case 'no'
      interp.mask(find(interp.mask<cfg.maskclipmin)) = nan;
    case 'yes'
      interp.mask(find(abs(interp.mask)<cfg.maskclipmin)) = nan;
    end
    fprintf('done\n');
  end
  if ~ischar(cfg.maskclipmax)
    fprintf('clipping mask maximum...');
    switch cfg.maskclipsym
    case 'no'
      interp.mask(find(interp.mask>cfg.maskclipmax)) = nan;
    case 'yes'
      interp.mask(find(abs(interp.mask)>cfg.maskclipmax)) = nan;
    end
    fprintf('done\n');
  end
end

% scale functional data
fprintf('scaling functional data...');
fmin = min(interp.source(:));
fmax = max(interp.source(:));
if ~ischar(cfg.colmin)
  fcolmin = cfg.colmin;
else
  if sign(fmin)==sign(fmax)
    fcolmin = fmin;
  else
    fcolmin = -max(abs([fmin,fmax]));
  end
end
if ~ischar(cfg.colmax)
  fcolmax = cfg.colmax;
else
  if sign(fmin)==sign(fmax)
    fcolmax = fmax;
  else
    fcolmax = max(abs([fmin,fmax]));
  end
end
interp.source = (interp.source-fcolmin)./(fcolmax-fcolmin);
if ~ischar(cfg.colmax)
  interp.source(find(interp.source>1)) = 1;
end
if ~ischar(cfg.colmin)
  interp.source(find(interp.source<0)) = 0;
end
fprintf('done\n');

% scale mask data
if maskdat
  fprintf('scaling mask data...');
  fmin = min(interp.mask(:));
  fmax = max(interp.mask(:));
  if ~ischar(cfg.maskcolmin)
    mcolmin = cfg.maskcolmin;
  else
    if sign(fmin)==sign(fmax)
      mcolmin = fmin;
    else
      mcolmin = -max(abs([fmin,fmax]));
    end
  end
  if ~ischar(cfg.maskcolmax)
    mcolmax = cfg.maskcolmax;
  else
    if sign(fmin)==sign(fmax)
      mcolmax = fmax;
    else
      mcolmax = max(abs([fmin,fmax]));
    end
  end
  interp.mask = (interp.mask-mcolmin)./(mcolmax-mcolmin);
  if ~ischar(cfg.maskcolmax)
    interp.mask(find(interp.mask>1)) = 1;
  end
  if ~ischar(cfg.maskcolmin)
    interp.mask(find(interp.mask<0)) = 0;
  end
  fprintf('done\n');
end

% merge anatomy, functional data and mask
fprintf('constructing overlay...');
if ischar(cfg.colormap)
  % replace string by colormap using standard MATLAB function
  cfg.colormap = colormap(cfg.colormap);
end
cmap = cfg.colormap;
cmaplength = size(cmap,1);
maskmap = cfg.maskmap(:);
maskmaplength = size(maskmap,1);
indslice = round(linspace(spacemin,spacemax,cfg.nslices));
nvox1 = length(1:cfg.resample:size(interp.anatomy,2));
nvox2 = length(1:cfg.resample:size(interp.anatomy,3));
if mod(cfg.rotate,2)
  dummy = nvox1;
  nvox1 = nvox2;
  nvox2 = dummy;
end
out = zeros(nvox1,nvox2,3,cfg.nslices);
for islice = 1:cfg.nslices
  sel1 = 1:cfg.resample:size(interp.anatomy,2);
  sel2 = 1:cfg.resample:size(interp.anatomy,3);
  
  dummy1 = reshape(interp.anatomy(indslice(islice),sel1,sel2), [numel(sel1) numel(sel2)]);
  dummy2 = reshape(interp.source(indslice(islice),sel1,sel2),  [numel(sel1) numel(sel2)]);
  indmarker = find(reshape(marker(indslice(islice),sel1,sel2), [numel(sel1) numel(sel2)]));
  indsource = find(~isnan(dummy2));
  if maskdat
    dummymask = reshape(interp.mask(indslice(islice),sel1,sel2), [numel(sel1) numel(sel2)]);
    indsource = find(~isnan(dummy2) & ~isnan(dummymask));
  end
  for icol = 1:3
    dummy3 = dummy1;
    if not(maskdat)
      if ~ischar(cfg.alpha)
        try
          dummy3(indsource) = ...
            (1-cfg.alpha) * dummy3(indsource) + ...
            cfg.alpha * cmap(round(dummy2(indsource)*(cmaplength-1))+1,icol);
        end
      else
        try
          dummy3(indsource) = ...
            (1-dummy2(indsource)) .* dummy3(indsource) + ...
            dummy2(indsource) .* cmap(round(dummy2(indsource)*(cmaplength-1))+1,icol);
        end
      end
    else
      dummy3(indsource) = ...
        (1-maskmap(round(dummymask(indsource)*(maskmaplength-1))+1)).* ...
        dummy3(indsource) + ...
        maskmap(round(dummymask(indsource)*(maskmaplength-1))+1) .* ...
        cmap(round(dummy2(indsource)*(cmaplength-1))+1,icol);
    end
    dummy3(indmarker) = cfg.markercolor(icol);
    out(:,:,icol,islice) = rot90(dummy3,cfg.rotate);
  end
  if strcmp(cfg.whitebg,'yes')
    bgmask = zeros(nvox1,nvox2);
    bgmask(find(conv2(mean(out(:,:,:,islice),3),ones(round((nvox1+nvox2)/8))/(round((nvox1+nvox2)/8).^2),'same')<0.1)) = 1;
    for icol = 1:3
      out(:,:,icol,islice) = bgmask.*ones(nvox1,nvox2) + (1-bgmask).* out(:,:,icol,islice);
    end
  end
end
fprintf('done\n');

clf;
fprintf('plotting...');
axes('position',[0.9 0.3 0.02 0.4]);
image(permute(cmap,[1 3 2]));
set(gca,'YAxisLocation','right');
set(gca,'XTick',[]);
set(gca,'YDir','normal');
set(gca,'YTick',linspace(1,cmaplength,5));
set(gca,'YTickLabel',linspace(fcolmin,fcolmax,5));
set(gca,'Box','on');
axes('position',[0.01 0.01 0.88 0.90]);
[h,nrows,ncols]=slicemon(out);
xlim=get(gca,'XLim');
ylim=get(gca,'YLim');
text(diff(xlim)/2,-diff(ylim)/100,cfg.title,'HorizontalAlignment','center','Interpreter','none');
drawnow;
fprintf('done\n');
if nargout > 0
  outim=get(h,'CData');
end

if strcmp(cfg.interactive,'yes')
  data.sin = size(interp.source);
  data.nrows = nrows;
  data.ncols = ncols;
  data.out = out;
  data.indslice = indslice;
  data.cfg = cfg;
  data.hfig = gcf;
  uicontrol('Units','norm', 'Position', [0.9 0.2 0.08 0.05], 'Style','pushbutton', 'String','coords',...
    'Callback',@getcoords,'FontSize',7);
  data.hcoords = uicontrol('Units','norm', 'Position', [0.9 0.05 0.08 0.13], 'Style','text', 'String','','HorizontalAlign','left','FontSize',7);
  guidata(data.hfig,data);
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble history ininterp
ft_postamble provenance


% ---------------- subfunctions ----------------

function getcoords(h,eventdata,handles,varargin)
data = guidata(gcf);
[xi,yi] = ginput(1);

co(2,1) = round(mod(yi,size(data.out,1)));
co(3,1) = round(mod(xi,size(data.out,2)));
switch mod(data.cfg.rotate,4)
case 1,
  t1 = co(2);
  co(2) = co(3);
  co(3) = data.sin(3)-t1;
case 2,
  co(2) = data.sin(2)-co(2);
  co(3) = data.sin(3)-co(3);
case 3,
  t1 = co(3);
  co(3) = co(2);
  co(2) = data.sin(2)-t1;
end

try
  co(1) = data.indslice(fix(xi/size(data.out,2)) + fix(yi/size(data.out,1))*data.ncols + 1);
catch
  co(1) = NaN;
end

if strcmp(data.cfg.flipdim, 'yes')
  co(1) = data.sin(1) - co(1) + 1;
end
co = co(:);

co(2:3) = round(co(2:3)*data.cfg.resample);
for ishift = 1:data.cfg.dim-1
 co = [co(3);co(1);co(2)];
end
set(data.hcoords,'String',sprintf('1: %d\n2: %d\n3: %d\nf: %0.4f',co(1),co(2),co(3),data.source(co(1),co(2),co(3))));

function [h,nrows,ncols] = slicemon(a) % display the montage w/o image_toolbox
siz = [size(a,1) size(a,2) size(a,4)];
nn = sqrt(prod(siz))/siz(2);
mm = siz(3)/nn;
if (ceil(nn)-nn) < (ceil(mm)-mm),
  nn = ceil(nn); mm = ceil(siz(3)/nn);
else
  mm = ceil(mm); nn = ceil(siz(3)/mm);
end
b = a(1,1);
b(1,1) = 0;
b = repmat(b, [mm*siz(1), nn*siz(2), size(a,3), 1]);
rows = 1:siz(1); cols = 1:siz(2);
for i=0:mm-1,
  for j=0:nn-1,
    k = j+i*nn+1;
    if k<=siz(3),
      b(rows+i*siz(1),cols+j*siz(2),:) = a(:,:,:,k);
    end
  end
end
hh = image(b);
axis image;
box off;
set(gca,'XTick',[],'YTick',[],'Visible','off');
if nargout > 0
  h = hh;
  nrows = mm;
  ncols = nn;
end
