function data = ft_math(cfg, varargin)

% FT_MATH performs mathematical operations on FieldTrip data structures,
% such as addition, subtraction, division, etc.
%
% Use as
%   data = ft_examplefunction(cfg, data1, data2, ...)
% with one or multiple FieldTrip data structures as input and where cfg is a
% configuration structure that should contain
%
%  cfg.operation  = string, can be 'add', 'subtract', 'divide', 'multiply'
%  cfg.parameter  = string, input data field on which the operation is performed
%
% If you specify only a single input data structure, the configuration should contain
%   cfg.value     = scalar value to be used in the operation
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_DATATYPE

% Copyright (C) 2012-2013, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial part deals with parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

revision = '$Id$';

ft_defaults                   % this ensures that the path is correct and that the ft_defaults global variable is available
ft_preamble init              % this will show the function help if nargin==0 and return an error
ft_preamble provenance        % this records the time and memory usage at teh beginning of the function
ft_preamble trackconfig       % this converts the cfg structure in a config object, which tracks the cfg options that are being used
ft_preamble debug
ft_preamble loadvar varargin  % this reads the input data in case the user specified the cfg.inputfile option

type = ft_datatype(varargin{1});
for i=1:length(varargin)
  % check that all data types are equal, and update old data structures
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', type);
end

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', {'operation', 'parameter'});

if length(varargin)>1
  % perform the operation on two or multiple input data structures
  cfg = ft_checkconfig(cfg, 'forbidden', {'value'});
else
  % the operation involves the data structure and the specified value
  % or the operation is a transformation such as log10
end

% this function only works for the upcoming (not yet standard) source representation without sub-structures
if ft_datatype(varargin{1}, 'source')
  % update the old-style beamformer source reconstruction
  for i=1:length(varargin)
    varargin{i} = ft_datatype_source(varargin{i}, 'version', 'upcoming');
  end
  if isfield(cfg, 'parameter') && length(cfg.parameter)>4 && strcmp(cfg.parameter(1:4), 'avg.')
    cfg.parameter = cfg.parameter(5:end); % remove the 'avg.' part
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~issubfield(varargin{1}, cfg.parameter)
  error('the requested parameter is not present in the data');
end

% ensure that the data in all inputs has the same channels, time-axis, etc.
tmpcfg = [];
tmpcfg.parameter = cfg.parameter;
[varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
% restore the provenance information
[cfg, varargin{:}] = rollback_provenance(cfg, varargin{:});

cfg.parameter = tmpcfg.parameter;

if isfield(varargin{1}, [cfg.parameter 'dimord'])
  dimord = varargin{1}.([cfg.parameter 'dimord']);
elseif isfield(varargin{1}, 'dimord')
  dimord = varargin{1}.dimord;
else
  error('the dimord of the requested parameter is unknown');
end

dimtok = tokenize(dimord, '_');

% this determines which descriptive fields will get copied over
haschan = any(strcmp(dimtok, 'chan'));
hasfreq = any(strcmp(dimtok, 'freq'));
hastime = any(strcmp(dimtok, 'time'));
haspos  = any(strcmp(dimtok, 'pos'));

% construct the output data structure
data = [];
if haschan
  data.label = varargin{1}.label;
end
if hasfreq
  data.freq = varargin{1}.freq;
end
if hastime
  data.time = varargin{1}.time;
end
if haspos
  data.pos = varargin{1}.pos;
end

fprintf('selecting %s from the first input argument\n', cfg.parameter);
tmp = getsubfield(varargin{1}, cfg.parameter);

if length(varargin)==1
  switch cfg.operation
    case 'add'
      fprintf('adding %f to the %s\n', cfg.value, cfg.parameter);
      tmp = tmp + cfg.value;
      
    case 'subtract'
      fprintf('subtracting %f from the %s\n', cfg.value, cfg.parameter);
      tmp = tmp - cfg.value;
      
    case 'multiply'
      fprintf('multiplying %s with %f\n', cfg.parameter, cfg.value);
      tmp = tmp .* cfg.value;
      
    case 'divide'
      fprintf('dividing %s by %f\n', cfg.parameter, cfg.value);
      tmp = tmp ./ cfg.value;
      
    case 'log10'
      fprintf('taking the log10 of %s\n', cfg.parameter);
      tmp = log10(tmp);
      
    otherwise
      error('unsupported operation "%s"', cfg.operation);
  end % switch
  
  
else
  
  switch cfg.operation
    case 'add'
      for i=2:length(varargin)
        fprintf('adding the %s input argument\n', nth(i));
        tmp = tmp + varargin{i}.(cfg.parameter);
      end
      
    case 'multiply'
      for i=2:length(varargin)
        fprintf('multiplying with the %s input argument\n', nth(i));
        tmp = tmp .* varargin{i}.(cfg.parameter);
      end
      
    case 'subtract'
      if length(varargin)>2
        error('the operation "%s" requires exactly 2 input arguments', cfg.operation);
      end
      fprintf('subtracting the 2nd input argument from the 1st\n');
      tmp = tmp - varargin{2}.(cfg.parameter);
            
    case 'divide'
      if length(varargin)>2
        error('the operation "%s" requires exactly 2 input arguments', cfg.operation);
      end
      fprintf('dividing the 1st input argument by the 2nd\n');
      tmp = tmp ./ varargin{2}.(cfg.parameter);
      
    otherwise
      error('unsupported operation "%s"', cfg.operation);
  end % switch
end % one or multiple input data structures

% store the result of the operation in the output structure
data = setsubfield(data, cfg.parameter, tmp);
data.dimord = dimord;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ft_postamble debug
ft_postamble trackconfig        % this converts the config object back into a struct and can report on the unused fields
ft_postamble provenance         % this records the time and memory at the end of the function, prints them on screen and adds this information together with the function name and matlab version etc. to the output cfg
ft_postamble previous varargin  % this copies the datain.cfg structure into the cfg.previous field. You can also use it for multiple inputs, or for "varargin"
ft_postamble history data       % this adds the local cfg structure to the output data structure, i.e. dataout.cfg = cfg
ft_postamble savevar data       % this saves the output data structure to disk in case the user specified the cfg.outputfile option

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = nth(n)
if rem(n,10)==1 && rem(n,100)~=11
  s = sprintf('%dst', n);
elseif rem(n,10)==2 && rem(n,100)~=12
  s = sprintf('%dnd', n);
elseif rem(n,10)==3 && rem(n,100)~=13
  s = sprintf('%drd', n);
else
  s = sprintf('%dth', n);
end
