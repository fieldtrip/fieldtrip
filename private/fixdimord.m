function [data] = fixdimord(data);

% FIXDIMORD ensures consistency between the dimord string and the axes
% that describe the data dimensions. The main purpose of this function
% is to ensure backward compatibility of all functions with data that has
% been processed by older FieldTrip versions
%
% Use as
%   [data] = fixdimord(data)
% This will modify the data.dimord field to ensure consistency.
% The name of the axis is the same as the name of the dimord, i.e. if
% dimord='freq_time', then data.freq and data.time should be present.
%
% The default dimensions in the data are described by
%  'time'
%  'freq'
%  'chan'
%  'refchan'
%  'rpt'
%  'subj'
%  'chancmb'
%  'rpttap'

% Copyright (C) 2009, Robert Oostenveld, Jan-Mathijs Schoffelen
%
% $Log: fixdimord.m,v $
% Revision 1.4  2009/10/01 12:23:10  jansch
% added {'pos'} as placeholder
%
% Revision 1.3  2009/08/03 15:53:32  ingnie
% fixed bug introduced in last revision, thanks to Esther
%
% Revision 1.2  2009/08/03 15:07:51  ingnie
% cut off dimord from source and volume data, dimord is not implemented in source and volume data yet, but some functions do add it which causes unexpected behavior
%
% Revision 1.1  2009/07/02 08:04:36  roboos
% moved fixdimord and fixinside from private to public
%
% Revision 1.8  2006/04/12 09:11:23  roboos
% added ; to the end of a line
%
% Revision 1.7  2006/04/10 12:13:52  roboos
% improved documentation, added refchan as dimension
%
% Revision 1.6  2006/03/10 12:34:38  roboos
% ensure that label cell-array is a column
%
% Revision 1.5  2006/02/28 11:22:11  roboos
% changed chancmb into chan
%
% Revision 1.4  2006/02/23 10:28:17  roboos
% changed dimord strings for consistency, changed toi and foi into time and freq, added fixdimord where neccessary
%
% Revision 1.3  2006/02/09 10:06:01  roboos
% do not automatically add all unknown axes, that is too messy and currently not yet desired
%
% Revision 1.2  2006/02/01 12:26:04  roboos
% made all uses of dimord consistent with the common definition of data dimensions, see the fixdimord() function
%
% Revision 1.1  2006/02/01 11:34:32  roboos
% new implementation
%

if strcmp('volume', datatype(data)) || strcmp('source', datatype(data));
  if isfield(data, 'dimord')
    % data should not have a dimord (is not implemented yet, but some
    % functions add a dimord to these data which leads to unexpected behavior)
    warning(sprintf('unexpected dimord "%s", dimord is removed from data', data.dimord));
    data = rmfield(data, 'dimord');
    return
  else
    %is okay
    return
  end
end

if ~isfield(data, 'dimord')
  if ~isfield(data, 'trial') || ~iscell(data.trial) || ...
      ~isfield(data, 'time')  || ~iscell(data.time)  || ...
      ~isfield(data, 'label') || ~iscell(data.label)
    error('The data does not contain a dimord, but it also does not resemble raw data');
  elseif isfield(data, 'topo')
    % the data resembles a component decomposition
    data.dimord = 'chan_comp';
  else
    % the data does not contain a dimord, but it resembles raw data -> that's ok
    return
  end
end

dimtok = tokenize(data.dimord, '_');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(dimtok)
  switch dimtok{i}
    case {'tim' 'time' 'toi' 'latency'}
      dimtok{i} = 'time';

    case {'frq' 'freq' 'foi' 'frequency'}
      dimtok{i} = 'freq';

    case {'sgn' 'label' 'chan'}
      dimtok{i} = 'chan';

    case {'rpt' 'trial'}
      dimtok{i} = 'rpt';

    case {'subj' 'subject'}
      dimtok{i} = 'subj';

    case {'comp'}
      % don't change, it is ok

    case {'sgncmb' 'labelcmb' 'chancmb'}
      dimtok{i} = 'chan';

    case {'rpttap'}
      % this is a 2-D field, coding trials and tapers along the same dimension
      % don't change, it is ok

    case {'refchan'}
      % don't change, it is ok

    case {'vox' 'repl' 'wcond'}
      % these are used in some fieldtrip functions, but are not considered standard
      warning(sprintf('unexpected dimord "%s"', data.dimord));

    case {'pos'}
      % this will be the future default for simple sources

    otherwise
      error(sprintf('unexpected dimord "%s"', data.dimord));

  end % switch dimtok
end % for length dimtok

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(data, 'tim'),         data.time      = data.tim         ; data = rmfield(data, 'tim')        ; end
if isfield(data, 'toi'),         data.time      = data.toi         ; data = rmfield(data, 'toi')        ; end
if isfield(data, 'latency'),     data.time      = data.latency     ; data = rmfield(data, 'latency')    ; end
if isfield(data, 'frq'),         data.freq      = data.frq         ; data = rmfield(data, 'frq')        ; end
if isfield(data, 'foi'),         data.freq      = data.foi         ; data = rmfield(data, 'foi')        ; end
if isfield(data, 'frequency'),   data.freq      = data.frequency   ; data = rmfield(data, 'frequency')  ; end
if isfield(data, 'sgn'),         data.label     = data.sgn         ; data = rmfield(data, 'sgn')        ; end
if isfield(data, 'chan'),        data.label     = data.chan        ; data = rmfield(data, 'chan')       ; end
% if isfield(data, 'trial'),         data.rpt     = data.trial         ; data = rmfield(data, 'trial')        ; end  % DO NOT CONVERT -> this is an exception
if isfield(data, 'subject'),     data.subj      = data.subject     ; data = rmfield(data, 'subject')    ; end
if isfield(data, 'sgncmb'),      data.labelcmb  = data.sgncmb      ; data = rmfield(data, 'sgncmb')     ; end
if isfield(data, 'chancmb'),     data.labelcmb  = data.chancmb     ; data = rmfield(data, 'chancmb')    ; end

% ensure that it is a column
if isfield(data, 'label')
  data.label = data.label(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if isfield(data, 'trial')
%   mat = data.trial;
% elseif isfield(data, 'individual')
%   mat = data.individual;
% elseif isfield(data, 'avg')
%   mat = data.avg;
% elseif isfield(data, 'crsspctrm')
%   mat = data.crsspctrm;
% elseif isfield(data, 'powspctrm')
%   mat = data.powspctrm;
% elseif isfield(data, 'fourierspctrm')
%   mat = data.fourierspctrm;
% end
%
% add the descriptive axis for each dimension
% for i=1:length(dimtok)
%   if isfield(data, dimtok{i})
%     % the dimension is already described with its own axis
%     % data = setfield(data, dimtok{i}, getfield(data, dimtok{i}));
%   else
%     % add an axis to the output data
%     data = setfield(data, dimtok{i}, 1:size(mat,i));
%   end
% end

% undo the tokenization
data.dimord = dimtok{1};
for i=2:length(dimtok)
  data.dimord = [data.dimord '_' dimtok{i}];
end

