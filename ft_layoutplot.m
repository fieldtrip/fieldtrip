function [cfg] = ft_layoutplot(cfg, data)

% FT_LAYOUTPLOT makes a figure with the 2-D layout of the channel positions
% for topoplotting and the individual channel axes (i.e. width and height
% of the subfigures) for multiplotting. A correct 2-D layout is a
% prerequisite  for plotting the topographical distribution of the
% potential or field distribution, or for plotting timecourses in a
% topographical arrangement.
%
% This function uses the same configuration options as prepare_layout and
% as the topoplotting and multiplotting functions. The difference is that
% this function plots the layout without any data, which facilitates
% the validation of your 2-D layout.
%
% Use as
%   ft_layoutplot(cfg, data)
%
% There are several ways in which a 2-D layout can be made: it can be read
% directly from a *.lay file, it can be created based on 3-D electrode or
% gradiometer positions in the configuration or in the data, or it can be
% created based on the specification of an electrode of gradiometer file.
%
% You can specify either one of the following configuration options
%   cfg.layout      = filename containg the layout
%   cfg.rotate      = number, rotation around the z-axis in degrees (default = [], which means automatic)
%   cfg.projection  = string, 2D projection method can be 'stereographic', 'ortographic', 'polar', 'gnomic' or 'inverse' (default = 'orthographic')
%   cfg.elec        = structure with electrode definition
%   cfg.grad        = structure with gradiometer definition
%   cfg.elecfile    = filename containing electrode definition
%   cfg.gradfile    = filename containing gradiometer definition
%   cfg.output      = filename to which the layout will be written (default = [])
%   cfg.montage     = 'no' or a montage structure (default = 'no')
%   cfg.image       = filename, use an image to construct a layout (e.g. usefull for ECoG grids)
%   cfg.visible     = string, 'yes' or 'no' whether figure will be visible (default = 'yes')
%   cfg.box         = string, 'yes' or 'no' whether box should be plotted around electrode (default = 'yes')
%   cfg.mask        = string, 'yes' or 'no' whether the mask should be plotted (default = 'yes')
%
% Alternatively the layout can be constructed from either
%   data.elec     structure with electrode positions
%   data.grad     structure with gradiometer definition
%
% Alternatively, you can specify
%   cfg.layout = 'ordered'
% which will give you a 2-D ordered layout. Note that this is only suited
% for multiplotting and not for topoplotting.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure.
%
% See also FT_PREPARE_LAYOUT, FT_TOPOPLOTER, FT_TOPOPLOTTFR, FT_MULTIPLOTER, FT_MULTIPLOTTFR

% Undocumented options
%   cfg.montage

% Copyright (C) 2006-2008, Robert Oostenveld
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
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% the data can be passed as input argument or can be read from disk
hasdata = exist('data', 'var');

if hasdata
  % check if the input data is valid for this function
  data = ft_checkdata(data);
end

% set the defaults
cfg.visible = ft_getopt(cfg, 'visible', 'yes');
cfg.box     = ft_getopt(cfg, 'box', 'yes');
cfg.mask    = ft_getopt(cfg, 'mask', 'yes');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract/generate layout information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lay = [];

% try to use the layout structure from input if specified
if isfield(cfg, 'layout')
  % brief check to determine if cfg.layout is a valid layout (lay) structre
  if isstruct(cfg.layout)
    if all(isfield(cfg.layout, {'pos';'width';'height';'label'}))
      lay = cfg.layout;
    end
  end
end

% otherwise create the layout structure
if isempty(lay)
  if hasdata
    lay = ft_prepare_layout(cfg, data);
  else
    lay = ft_prepare_layout(cfg);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot all details pertaining to the layout in one figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if istrue(cfg.visible)
  f = figure('visible', 'on');
else
  f = figure('visible', 'off');
end

% set the figure window title
funcname = mfilename();
if hasdata
  if isfield(cfg, 'inputfile') && ~isempty(cfg.inputfile)
    dataname = cfg.inputfile;
  else
    dataname = inputname(2);
  end
else
  dataname = [];
end

if ~isempty(dataname)
  set(gcf, 'Name', sprintf('%d: %s: %s', double(gcf), funcname, dataname));
else
  set(gcf, 'Name', sprintf('%d: %s', double(gcf), funcname));
end
set(gcf, 'NumberTitle', 'off');

if isfield(cfg, 'image') && ~isempty(cfg.image)
  % start with the background image
  fprintf('reading background image from %s\n', cfg.image);
  img = imread(cfg.image);
  img = flipdim(img, 1); % in combination with "axis xy"
  
  bw = 1;
  
  if bw
    % convert to greyscale image
    img = mean(img, 3);
    imagesc(img);
    colormap gray
  else
    % plot as RGB image
    image(img);
  end
  axis equal
  axis off
  axis xy
end

ft_plot_lay(lay, 'point', true, 'box', istrue(cfg.box), 'label', true, 'mask', istrue(cfg.mask), 'outline', true);

% the following code can be used to verify a bipolar montage, given the
% layout of the monopolar channels
if isfield(cfg, 'montage') && ~isempty(cfg.montage)
  fprintf('plotting an arrow for each of the bipolar electrode pairs\n');
  % the arrow begins at the +1 electrode
  % the arrow ends   at the -1 electrode
  for i=1:length(cfg.montage.labelnew)
    begindx = find(cfg.montage.tra(i,:)==+1);
    endindx = find(cfg.montage.tra(i,:)==-1);
    if ~numel(begindx)==1 || ~numel(endindx)==1
      % the re-referenced channel does not seem to be a bipolar pair
      continue
    end
    % find the position of the begin and end of the arrow
    beglab = cfg.montage.labelold{begindx};
    endlab = cfg.montage.labelold{endindx};
    begindx = find(strcmp(lay.label, beglab)); % the index in the layout
    endindx = find(strcmp(lay.label, endlab)); % the index in the layout
    if ~numel(begindx)==1 || ~numel(endindx)==1
      % one of the channels in the bipolar pair does not seem to be in the layout
      continue
    end
    
    begpos = lay.pos(begindx,:);
    endpos = lay.pos(endindx,:);
    arrow(begpos, endpos, 'Length', 5)
    
  end % for all re-referenced channels
end % if montage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
ft_postamble previous data
