function [linecolor, linestyle, linewidth] = lineattributes_common(cfg, varargin)

% LINEATTRIBUTES_COMMON implements consistent line attributes for multiple channels/conditions
% 
% This function is used by 
%   ft_databrowser
%   ft_multiplotER
%   ft_singleplotER
% 
% It is not yet used by
%   ft_connectivityplot
% 
% Use as
% 
% [linecolor, linestyle, linewidth] = lineattributes_common(cfg, varargin)
% 
% The input varargin are the data object(s) which are the input of the caller function.
% The output consists of:
% 
%   linecolor = N x 3 x M matrix with rgb-values for N channels and M data arguments
%   linestyle = N x M cell array with linestyle for N channels and M data arguments
%   linewidth = N x M matrix with linewidth for N channels and M data arguments  
% 
% The configuration can have the following parameters:
%   cfg.colorgroups = char or numeric vector determining whether the
%                       different values are going to be distributed across channels
%                       ('sequential'), or across data arguments ('condition'). Other
%                       possibilities are 'allblacks', 'chantype', 'labelcharI', where
%                       I is a scalar number indicating the I'th character of the label
%                       based on which the grouping is done
%   cfg.stylegroups = char or numeric vector, same possibilities as above, save for 'allblacks'
%   cfg.widthgroups = char or numeric vector, same possibilities as above, save for 'allblacks'
%   cfg.linecolor   = char, Nx3 matrix, or Nx3xM matrix
%   cfg.linestyle   = char, or cell-array
%   cfg.linewidth   = scalar, or NxM matrix
%
% If cfg.linecolor is a char, it should either be a sequence of characters that can be translated into
% and rgb value (i.e. any of 'rbgcmykw'), or it can be 'spatial', in which case a color will be assigned
% based on the layout.color field. Typically, this will be a color that is based on the x/y/z position of 
% the corresponding sensor.

% Copyright (C) 2022, Donders Centre for Cognitive Neuroimaging
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

cfg.linecolor      = ft_getopt(cfg, 'linecolor',   []);
cfg.linestyle      = ft_getopt(cfg, 'linestyle',   '-');
cfg.linewidth      = ft_getopt(cfg, 'linewidth',   0.5);
cfg.colorgroups    = ft_getopt(cfg, 'colorgroups', 'condition'); 
cfg.stylegroups    = ft_getopt(cfg, 'stylegroups', 'condition');
cfg.widthgroups    = ft_getopt(cfg, 'widthgroups', 'condition');

% it is assumed here that all datasets have the same channels
label = varargin{1}.label;
Nchan = length(label);
Ndata = length(varargin);

% not all options can be combined
if isfield(cfg, 'maskstyle') && strcmp(cfg.maskstyle, 'difference') && ~strcmp(cfg.colorgroups, 'condition')
  ft_error('cfg.colorgroups=''%s'' is not supported', cfg.colorgroups);
end

% not all options can be combined
if Ndata>1 && ~strcmp(cfg.colorgroups, 'condition')
  ft_error('cfg.colorgroups=''%s'' is not supported', cfg.colorgroups);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linecolor
if isnumeric(cfg.linecolor) && size(cfg.linecolor,2)==3
  % overrule the default colors, because it seems to be Nx3 numeric
elseif ischar(cfg.linecolor) && all(ismember(cfg.linecolor, 'rgbcmykw'))
  % this is handled below
elseif isequal(cfg.linecolor, 'spatial')
  if isfield(cfg, 'layout') && isfield(cfg.layout, 'color')
    cfg.linecolor = cfg.layout.color;
  else
    ft_error('this does not work yet'); % FIXME in principle the color can be derived from the 3D chanpos if present in the data
  end
elseif isempty(cfg.linecolor)
  cfg.linecolor = [
    0    0    1
    0.75 0    0
    0    1    0
    0.44 0.19 0.63
    0    0.13 0.38
    0.5  0.5  0.5
    1    0.75 0
    1    0    0
    0.89 0.42 0.04
    0.85 0.59 0.58
    0.57 0.82 0.31
    0    0.69 0.94
    1    0    0.4
    0    0.69 0.31
    0    0.44 0.75];
end

if ischar(cfg.linecolor)
  % convert to rgb
  cfg.linecolor = char2rgb(cfg.linecolor);
end

% ensure it is a Nx3xM matrix, the right size will be dealt with below
if ndims(cfg.linecolor)==3
  % if cfg.linecolor is already 3D
  linecolor = repmat(cfg.linecolor, [ceil(Nchan/size(cfg.linecolor,1)) 1 ceil(Ndata/size(cfg.linecolor,3))]);
  linecolor = linecolor(1:Nchan, :, 1:Ndata);
  
else 
  % the operation depends on the cfg.colorgroups specs
  if isnumeric(cfg.colorgroups)
    % this is a per channel color definition, if more data arguments are defined, use the same color for the other inputs
    
    % channel groups are defined by the user
    if Nchan ~= length(cfg.colorgroups)
      ft_error('length(cfg.colorgroups) should correspond to the number of channels')
    end
    linecolor = repmat(cfg.linecolor(cfg.colorgroups(:),:), [1 1 Ndata]);
    
  elseif strcmp(cfg.colorgroups, 'condition')
    % no grouping of channels, each condition has a single color
    % ensure that the number of colors matches the number of conditions
    linecolor = repmat(cfg.linecolor, ceil(Ndata/size(cfg.linecolor,1)), 1);
    
    % make 3D    
    linecolor = repmat(shiftdim(linecolor(1:Ndata, :)', -1), [Nchan 1 1]);
    
  elseif strcmp(cfg.colorgroups, 'sequential')
    % each channel has its color, same across conditions
    linecolor = repmat(cfg.linecolor, ceil(Nchan/size(cfg.linecolor,1)), 1);
    
    % make 3D
    linecolor = repmat(linecolor(1:Nchan,:), [1 1 Ndata]);
    
  elseif strcmp(cfg.colorgroups, 'allblack')
    
    linecolor = zeros(Nchan, 3, Ndata);
    
  elseif strcmp(cfg.colorgroups, 'chantype')
    % channel groups are defined by the chantype
    type = ft_chantype(varargin{1});
    [tmp1, tmp2, colorgroups] = unique(type);
    fprintf('%3d colorgroups were identified\n',length(tmp1));
    linecolor = cfg.linecolor(colorgroups(:),:);
  
    % make 3D
    linecolor = repmat(linecolor, [1 1 Ndata]);
    
  elseif startsWith(cfg.colorgroups, 'labelchar')
    % channel groups are defined by the Nth letter of the channel label
    labelchar_num = sscanf(cfg.colorgroups, 'labelchar%d');
    vec_letters = num2str(zeros(Nchan,1));
    for iChan = 1:Nchan
      vec_letters(iChan) = label{iChan}(labelchar_num);
    end
    [tmp1, tmp2, colorgroups] = unique(vec_letters);
    fprintf('%3d colorgroups were identified\n', length(tmp1))
    linecolor = cfg.linecolor(colorgroups(:),:);
    
    % make 3D
    linecolor = repmat(linecolor, [1 1 Ndata]);

  else
    ft_error('unsupported value for cfg.colorgroups');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linestyle
if isempty(cfg.linestyle)
  cfg.linestyle = {'-' '--' ':' '-.'};
end

if ischar(cfg.linestyle)
  % convert to cell
  cfg.linestyle = {cfg.linestyle};
end

% ensure it is a Nx3xM matrix, the right size will be dealt with below
if numel(cfg.linestyle)==1
  % if cfg.linestyle is the same for all channels and conditions
  linestyle = repmat(cfg.linestyle, [ceil(Nchan/size(cfg.linestyle,1)) ceil(Ndata/size(cfg.linestyle,3))]);
  linestyle = linestyle(1:Nchan, 1:Ndata);
  
else 
  % the operation depends on the cfg.colorgroups specs
  if isnumeric(cfg.stylegroups)
    % this is a per channel style definition, if more data arguments are defined, use the same style for the other inputs
    
    % channel groups are defined by the user
    if Nchan ~= length(cfg.stylegroups)
      ft_error('length(cfg.stylegroups) should correspond to the number of channels')
    end
    linestyle = repmat(cfg.linestyle(cfg.stylegroups(:),:), [1 Ndata]);
    
  elseif strcmp(cfg.stylegroups, 'condition')
    % no grouping of channels, each condition has a single style
    % ensure that the number of styles matches the number of conditions
    linestyle = repmat(cfg.linestyle(:)', 1, ceil(Ndata/size(cfg.linestyle,1)));
    
    % make 2D    
    linestyle = repmat(linestyle(:, 1:Ndata), [Nchan 1]);
    
  elseif strcmp(cfg.stylegroups, 'sequential')
    % each channel has its style, same across conditions
    linestyle = repmat(cfg.linestyle(:), ceil(Nchan/size(cfg.linestyle,1)), 1);
    
    % make 2D
    linestyle = repmat(linestyle(1:Nchan), [1 Ndata]);
      
  elseif strcmp(cfg.stylegroups, 'chantype')
    % channel groups are defined by the chantype
    type = ft_chantype(varargin{1});
    [tmp1, tmp2, stylegroups] = unique(type);
    fprintf('%3d stylegroups were identified\n',length(tmp1));
    linestyle = cfg.linestyle(stylegroups(:),:);
  
    % make 3D
    linestyle = repmat(linestyle, [1 1 Ndata]);
    
  elseif startsWith(cfg.stylegroups, 'labelchar')
    % channel groups are defined by the Nth letter of the channel label
    labelchar_num = sscanf(cfg.stylegroups, 'labelchar%d');
    vec_letters = num2str(zeros(Nchan,1));
    for iChan = 1:Nchan
      vec_letters(iChan) = label{iChan}(labelchar_num);
    end
    [tmp1, tmp2, stylegroups] = unique(vec_letters);
    fprintf('%3d stylegroups were identified\n', length(tmp1))
    linestyle = cfg.linestyle(stylegroups(:),:);
    
    % make 3D
    linestyle = repmat(linestyle, [1 1 Ndata]);

  else
    ft_error('unsupported value for cfg.stylegroups');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linewidth
if isempty(cfg.linewidth)
  cfg.linewidth = 0.5;
end

% ensure it is a Nx3xM matrix, the right size will be dealt with below
if numel(cfg.linewidth)==1
  % if cfg.linewidth is the same for all channels and conditions
  linewidth = repmat(cfg.linewidth, [ceil(Nchan/size(cfg.linewidth,1)) ceil(Ndata/size(cfg.linewidth,3))]);
  linewidth = linewidth(1:Nchan, 1:Ndata);
  
else 
  % the operation depends on the cfg.colorgroups specs
  if isnumeric(cfg.widthgroups)
    % this is a per channel width definition, if more data arguments are defined, use the same width for the other inputs
    
    % channel groups are defined by the user
    if Nchan ~= length(cfg.widthgroups)
      ft_error('length(cfg.widthgroups) should correspond to the number of channels')
    end
    linewidth = repmat(cfg.linewidth(cfg.widthgroups(:),:), [1 Ndata]);
    
  elseif strcmp(cfg.widthgroups, 'condition')
    % no grouping of channels, each condition has a single width
    % ensure that the number of widths matches the number of conditions
    linewidth = repmat(cfg.linewidth(:)', 1, ceil(Ndata/size(cfg.linewidth,1)));
    
    % make 2D    
    linewidth = repmat(linewidth(:, 1:Ndata), [Nchan 1]);
    
  elseif strcmp(cfg.widthgroups, 'sequential')
    % each channel has its width, same across conditions
    linewidth = repmat(cfg.linewidth(:), ceil(Nchan/size(cfg.linewidth,1)), 1);
    
    % make 2D
    linewidth = repmat(linewidth(1:Nchan), [1 Ndata]);
      
  elseif strcmp(cfg.widthgroups, 'chantype')
    % channel groups are defined by the chantype
    type = ft_chantype(varargin{1});
    [tmp1, tmp2, widthgroups] = unique(type);
    fprintf('%3d widthgroups were identified\n',length(tmp1));
    linewidth = cfg.linewidth(widthgroups(:),:);
  
    % make 3D
    linewidth = repmat(linewidth, [1 1 Ndata]);
    
  elseif startsWith(cfg.widthgroups, 'labelchar')
    % channel groups are defined by the Nth letter of the channel label
    labelchar_num = sscanf(cfg.widthgroups, 'labelchar%d');
    vec_letters = num2str(zeros(Nchan,1));
    for iChan = 1:Nchan
      vec_letters(iChan) = label{iChan}(labelchar_num);
    end
    [tmp1, tmp2, widthgroups] = unique(vec_letters);
    fprintf('%3d widthgroups were identified\n', length(tmp1))
    linewidth = cfg.linewidth(widthgroups(:),:);
    
    % make 3D
    linewidth = repmat(linewidth, [1 1 Ndata]);

  else
    ft_error('unsupported value for cfg.widthgroups');
  end
end
