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
  % ensure it is a column vector of the right length
  cfg.linecolor = repmat(cfg.linecolor(:), ceil(Ncolor/length(cfg.linecolor)), 1);
  cfg.linecolor = cfg.linecolor(1:Ncolor);
elseif isnumeric(cfg.linecolor)
  % ensure it is a Nx3 matrix of the right length
  cfg.linecolor = repmat(cfg.linecolor, ceil(Ncolor/size(cfg.linecolor,1)), 1);
  cfg.linecolor = cfg.linecolor(1:Ncolor,:);
end

if isnumeric(cfg.colorgroups)
  % channel groups are defined by the user
  if Nchan ~= length(cfg.colorgroups)
    ft_error('length(cfg.colorgroups) should correspond to the number of channels')
  end
  linecolor = cfg.linecolor(cfg.colorgroups(:),:);
elseif strcmp(cfg.colorgroups, 'condition')
  % no grouping of channels, each condition has a single color
  linecolor = cfg.linecolor;
elseif strcmp(cfg.colorgroups, 'sequential')
  % no grouping of channels
  linecolor = cfg.linecolor;
elseif strcmp(cfg.colorgroups, 'allblack')
  % all channels/components are show in black
  linecolor = zeros(Ncolor,3);
elseif strcmp(cfg.colorgroups, 'chantype')
  % channel groups are defined by the chantype
  type = ft_chantype(varargin{1});
  [tmp1, tmp2, cfg.colorgroups] = unique(type);
  fprintf('%3d colorgroups were identified\n',length(tmp1))
  linecolor = cfg.linecolor(cfg.colorgroups(:),:);
elseif startsWith(cfg.colorgroups, 'labelchar')
  % channel groups are defined by the Nth letter of the channel label
  labelchar_num = str2double(cfg.colorgroups(10));
  vec_letters = num2str(zeros(Nchan,1));
  for iChan = 1:Nchan
    vec_letters(iChan) = label{iChan}(labelchar_num);
  end
  [tmp1, tmp2, cfg.colorgroups] = unique(vec_letters);
  fprintf('%3d colorgroups were identified\n', length(tmp1))
  linecolor = cfg.linecolor(cfg.colorgroups(:),:);
else
  ft_error('do not understand cfg.colorgroups')
end
