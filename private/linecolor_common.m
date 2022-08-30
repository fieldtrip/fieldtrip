function linecolor = linecolor_common(cfg, varargin)

% LINECOLOR_COMMON implements consistent color handling for multiple channels/conditions
%
% This function is used by 
%   ft_databrowser
%   ft_multiplotER
%   ft_singleplotER
%
% It is not yet used by
%   ft_connectivityplot
%
% The goal is to create a N x 3 x M matrix of rgb-values, where
% where N = Nchan, and M = Ncond, depending on the specification
% of cfg.colorgroups

cfg.linecolor   = ft_getopt(cfg, 'linecolor', []);
cfg.colorgroups = ft_getopt(cfg, 'colorgroups', 'condition'); 

if isempty(cfg.linecolor)
  cfg.linecolor = [0.75 0 0; 
    0 0 1; 
    0 1 0; 
    0.44 0.19 0.63; 
    0 0.13 0.38;
    0.5 0.5 0.5;
    1 0.75 0;
    1 0 0;
    0.89 0.42 0.04;
    0.85 0.59 0.58;
    0.57 0.82 0.31;
    0 0.69 0.94;
    1 0 0.4;
    0 0.69 0.31;
    0 0.44 0.75];
end


Ndata = length(varargin);

% it is assumed here that all datasets have the same channels
label = varargin{1}.label;
Nchan = length(label);

% not all options can be combined
if isfield(cfg, 'maskstyle') && strcmp(cfg.maskstyle, 'difference') && ~strcmp(cfg.colorgroups, 'condition')
  ft_error('cfg.colorgroups=''%s'' is not supported', cfg.colorgroups);
end

% not all options can be combined
if Ndata>1 && ~strcmp(cfg.colorgroups, 'condition')
  ft_error('cfg.colorgroups=''%s'' is not supported', cfg.colorgroups);
end

if Ndata==1
  Ncolor = Nchan;
else
  Ncolor = Ndata;
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
