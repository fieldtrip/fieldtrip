function [data] = ft_math(cfg, varargin)

% FT_MATH performs mathematical operations on FieldTrip data structures,
% such as addition, subtraction, division, etc.
%
% Use as
%   data = ft_math(cfg, data1, data2, ...)
% with one or multiple FieldTrip data structures as the input and the configuration
% structure cfg in which you specify the mathematical operation that is to be
% executed on the desired parameter from the data
%   cfg.parameter = string, field from the input data on which the operation is
%                   performed, e.g. 'pow' or 'avg'
%   cfg.operation = string, for example '(x1-x2)/(x1+x2)' or 'x1/6'
%
% In the specification of the mathematical operation, x1 is the parameter obtained
% from the first input data structure, x2 from the second, etc.
%
% Rather than specifying the operation as a string that is evaluated, you can also
% specify it as a single operation. The advantage is that it is computed faster.
%    cfg.operation = string, can be 'add', 'subtract', 'divide', 'multiply', 'log10', 'abs'
% If you specify only a single input data structure and the operation is 'add',
% 'subtract', 'divide' or 'multiply', the configuration should also contain:
%   cfg.scalar    = scalar value to be used in the operation
%
% The operation 'add' is implemented as follows
%   y = x1 + x2 + ....
% if you specify multiple input arguments, or as
%   y = x1 + s
% if you specify one input argument and a scalar value.
%
% The operation 'subtract' is implemented as follows
%   y = x1 - x2 - ....
% if you specify multiple input arguments, or as
%   y = x1 - s
% if you specify one input argument and a scalar value.
%
% The operation 'divide' is implemented as follows
%   y = x1 ./ x2
% if you specify two input arguments, or as
%   y = x1 / s
% if you specify one input argument and a scalar value.
%
% The operation 'multiply' is implemented as follows
%   y = x1 .* x2
% if you specify two input arguments, or as
%   y = x1 * s
% if you specify one input argument and a scalar value.
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

% Undocumented options:
%   cfg.matrix = rather than using a scalar, a matrix can be specified. In
%                this case, the dimensionality of cfg.matrix should be equal
%                to the dimensionality of data.(cfg.parameter). If used in
%                combination with cfg.operation, the operation should
%                involve element-wise combination of the data and the
%                matrix.

% Copyright (C) 2012-2015, Robert Oostenveld
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial part deals with parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar varargin
ft_preamble provenance varargin
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

type = ft_datatype(varargin{1});
for i=1:length(varargin)
  % check if the input data is valid for this function, that all data types are equal and update old data structures
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', type);
end

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', {'operation', 'parameter'});
cfg = ft_checkconfig(cfg, 'renamed', {'value', 'scalar'});
cfg = ft_checkconfig(cfg, 'renamedval', {'funparameter', 'avg.pow', 'pow'});
cfg = ft_checkconfig(cfg, 'renamedval', {'funparameter', 'avg.coh', 'coh'});
cfg = ft_checkconfig(cfg, 'renamedval', {'funparameter', 'avg.mom', 'mom'});

if ~iscell(cfg.parameter)
  cfg.parameter = {cfg.parameter};
end

if ft_datatype(varargin{1}, 'raw+comp')
    if length(varargin)>1
        ft_error('ft_math does not support more than one input argument if the input data is of type "raw" or "comp"')
    end
end

% this function only works for the upcoming (not yet standard) source representation without sub-structures
if ft_datatype(varargin{1}, 'source')
  % update the old-style beamformer source reconstruction
  for i=1:length(varargin)
    varargin{i} = ft_datatype_source(varargin{i}, 'version', 'upcoming');
  end
  for p = 1:length(cfg.parameter)
    if strncmp(cfg.parameter{p}, 'avg.', 4)
      cfg.parameter{p} = cfg.parameter{p}(5:end); % remove the 'avg.' part
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for p=1:length(cfg.parameter)
  if ~issubfield(varargin{1}, cfg.parameter{p})
    ft_error('the requested parameter is not present in the data');
  end
end

% ensure that the data in all inputs has the same channels, time-axis, etc.
tmpcfg = [];
tmpcfg.parameter = cfg.parameter;
[varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
% restore the provenance information
[cfg, varargin{:}] = rollback_provenance(cfg, varargin{:});
% restore the user-specified parameter option
cfg.parameter = tmpcfg.parameter;

for p = 1:length(cfg.parameter)
  dimordtmp{p} = getdimord(varargin{1}, cfg.parameter{p});
  if p>1 && ~strcmp(dimordtmp{1}, dimordtmp{p})
    ft_error('the dimord of multiple parameters must be the same');
  end
end
clear dimordtmp

% construct the output data structure; make sure descriptive fields will get copied over
% some ugly things need to be done in order to get the correct xxxdimord
% fields in the output
fn  = fieldnames(varargin{1});
dimordfields = fn(~cellfun(@isempty, strfind(fn, 'dimord')))';
if numel(dimordfields)==1 && strcmp(dimordfields{1},'dimord')
    % this is OK and counts for most data structures
else
    % this is in the case of one or more xxxdimord fields, in which case
    % only the requested parameters' xxxdimord fields should be returned in
    % the output
    ok = false(1,numel(dimordfields));
    for p = 1:length(cfg.parameter)
        ok(p) = any(~cellfun(@isempty, strfind(dimordfields, cfg.parameter{p})));
    end
    dimordfields = dimordfields(ok);
end
data = keepfields(varargin{1}, [dimordfields {'label', 'labelcmb', 'freq', 'time', 'pos', 'dim', 'transform'}]);

for p = 1:length(cfg.parameter)
  fprintf('selecting %s from the first input argument\n', cfg.parameter{p});
  % create the local variables x1, x2, ...
  for i=1:length(varargin)
    assign_var(sprintf('x%i', i), getsubfield(varargin{i}, cfg.parameter{p}));
  end

  % create the local variables s and m
  s = ft_getopt(cfg, 'scalar');
  m = ft_getopt(cfg, 'matrix');

  % check the dimensionality of m against the input data
  if ~isempty(m)
    for i=1:length(varargin)
      ok = isequal(size(getsubfield(varargin{i}, cfg.parameter{p})),size(m));
      if ~ok, break; end
    end
    if ~ok
      ft_error('the dimensions of cfg.matrix do not allow for element-wise operations');
    end
  end

  % only one of these can be defined at the moment (i.e. not allowing for
  % operations such as (x1+m)^s for now
  if ~isempty(m) && ~isempty(s)
    ft_error('you can either specify a cfg.matrix or a cfg.scalar, not both');
  end

  % touch it to keep track of it in the output cfg
  if ~isempty(s), cfg.scalar; end
  if ~isempty(m), cfg.matrix; end

  % replace s with m, so that the code below is more transparent
  if ~isempty(m)
    s = m; clear m;
  end

  if length(varargin)==1
    switch cfg.operation
      case 'add'
        if isscalar(s)
          fprintf('adding %f to the %s\n', s, cfg.parameter{p});
        else
          fprintf('adding the contents of cfg.matrix to the %s\n', cfg.parameter{p});
        end
        if iscell(x1)
          y = cellplus(x1, s);
        else
          y = x1 + s;
        end

      case 'subtract'
        if isscalar(s)
          fprintf('subtracting %f from the %s\n', s, cfg.parameter{p});
        else
          fprintf('subtracting the contents of cfg.matrix from the %s\n', cfg.parameter{p});
        end
        if iscell(x1)
          y = cellminus(x1, s);
        else
          y = x1 - s;
        end

      case 'multiply'
        if isscalar(s)
          fprintf('multiplying %s with %f\n', cfg.parameter{p}, s);
        else
          fprintf('multiplying %s with the content of cfg.matrix\n', cfg.parameter{p});
        end
        fprintf('multiplying %s with %f\n', cfg.parameter{p}, s);
        if iscell(x1)
          y = celltimes(x1, s);
        else
          y = x1 .* s;
        end

      case 'divide'
        if isscalar(s)
          fprintf('dividing %s by %f\n', cfg.parameter{p}, s);
        else
          fprintf('dividing %s by the content of cfg.matrix\n', cfg.parameter{p});
        end
        if iscell(x1)
          y = cellrdivide(x1, s);
        else
          y = x1 ./ s;
        end

      case 'log10'
        fprintf('taking the log10 of %s\n', cfg.parameter{p});
        if iscell(x1)
          y = celllog10(x1);
        else
          y = log10(x1);
        end

      case 'abs'
        fprintf('taking the abs of %s\n', cfg.parameter{p});
        if iscell(x1)
          y = cellabs(x1);
        else
          y = abs(x1);
        end

      otherwise
        % assume that the operation is descibed as a string, e.g. x1^s
        % where x1 is the first argument and s is obtained from cfg.scalar

        arginstr = sprintf('x%i,', 1:length(varargin));
        arginstr = arginstr(1:end-1); % remove the trailing ','
        eval(sprintf('operation = @(%s) %s;', arginstr, cfg.operation));

        if ~iscell(varargin{1}.(cfg.parameter{p}))
          % gather x1, x2, ... into a cell-array
          arginval = eval(sprintf('{%s}', arginstr));
          eval(sprintf('operation = @(%s) %s;', arginstr, cfg.operation));
          if numel(s)<=1
            y = arrayfun(operation, arginval{:});
          elseif size(s)==size(arginval{1})
            y = feval(operation, arginval{:});
          end
        else
          y = cell(size(x1));
          % do the same thing, but now for each element of the cell-array
          for i=1:numel(y)
            for j=1:length(varargin)
              % rather than working with x1 and x2, we need to work on its elements
              % xx1 is one element of the x1 cell-array
              assign_var(sprintf('xx%d', j), eval(sprintf('x%d{%d}', j, i)))
            end

            % gather xx1, xx2, ... into a cell-array
            arginstr = sprintf('xx%i,', 1:length(varargin));
            arginstr = arginstr(1:end-1); % remove the trailing ','
            arginval = eval(sprintf('{%s}', arginstr));
            if numel(s)<=1
              y{i} = arrayfun(operation, arginval{:});
            else
              y{i} = feval(operation, arginval{:});
            end
          end % for each element
        end % iscell or not

    end % switch


  else

    switch cfg.operation
      case 'add'
        for i=2:length(varargin)
          fprintf('adding the %s input argument\n', nth(i));
          if iscell(x1)
            y = cellplus(x1, varargin{i}.(cfg.parameter{p}));
          else
            y = x1 + varargin{i}.(cfg.parameter{p});
          end
        end

      case 'multiply'
        for i=2:length(varargin)
          fprintf('multiplying with the %s input argument\n', nth(i));
          if iscell(x1)
            y = celltimes(x1, varargin{i}.(cfg.parameter{p}));
          else
            y = x1 .* varargin{i}.(cfg.parameter{p});
          end
        end

      case 'subtract'
        if length(varargin)>2
          ft_error('the operation "%s" requires exactly 2 input arguments', cfg.operation);
        end
        fprintf('subtracting the 2nd input argument from the 1st\n');
        if iscell(x1)
          y = cellminus(x1, varargin{2}.(cfg.parameter{p}));
        else
          y = x1 - varargin{2}.(cfg.parameter{p});
        end

      case 'divide'
        if length(varargin)>2
          ft_error('the operation "%s" requires exactly 2 input arguments', cfg.operation);
        end
        fprintf('dividing the 1st input argument by the 2nd\n');
        if iscell(x1)
          y = cellrdivide(x1, varargin{2}.(cfg.parameter{p}));
        else
          y = x1 ./ varargin{2}.(cfg.parameter{p});
        end

      case 'log10'
        if length(varargin)>2
          ft_error('the operation "%s" requires exactly 2 input arguments', cfg.operation);
        end
        fprintf('taking the log difference between the 2nd input argument and the 1st\n');
        y = log10(x1 ./ varargin{2}.(cfg.parameter{p}));

      otherwise
        % assume that the operation is descibed as a string, e.g. (x1-x2)/(x1+x2)

        % ensure that all input arguments are being used
        for i=1:length(varargin)
          assert(~isempty(regexp(cfg.operation, sprintf('x%i', i), 'once')), 'not all input arguments are assigned in the operation')
        end

        arginstr = sprintf('x%i,', 1:length(varargin));
        arginstr = arginstr(1:end-1); % remove the trailing ','
        eval(sprintf('operation = @(%s) %s;', arginstr, cfg.operation));

        if ~iscell(varargin{1}.(cfg.parameter{p}))
          % gather x1, x2, ... into a cell-array
          arginval = eval(sprintf('{%s}', arginstr));
          eval(sprintf('operation = @(%s) %s;', arginstr, cfg.operation));
          if numel(s)<=1
            y = arrayfun(operation, arginval{:});
          else
            y = feval(operation, arginval{:});
          end
        else
          y = cell(size(x1));
          % do the same thing, but now for each element of the cell-array
          for i=1:numel(y)
            for j=1:length(varargin)
              % rather than working with x1 and x2, we need to work on its elements
              % xx1 is one element of the x1 cell-array
              assign_var(sprintf('xx%d', j), eval(sprintf('x%d{%d}', j, i)))
            end

            % gather xx1, xx2, ... into a cell-array
            arginstr = sprintf('xx%i,', 1:length(varargin));
            arginstr = arginstr(1:end-1); % remove the trailing ','
            arginval = eval(sprintf('{%s}', arginstr));
            if numel(s)<=1
              y{i} = arrayfun(operation, arginval{:});
            else
              y{i} = feval(operation, arginval{:});
            end
          end % for each element
        end % iscell or not

    end % switch
  end % one or multiple input data structures

  % store the result of the operation in the output structure
  data = setsubfield(data, cfg.parameter{p}, y);
end % p over length(cfg.parameter)

% certain fields should remain in the output, but only if they are identical in all inputs
keepfield = {'grad', 'elec', 'opto', 'inside', 'trialinfo', 'sampleinfo', 'tri'};
for j=1:numel(keepfield)
  if isfield(varargin{1}, keepfield{j})
    tmp  = varargin{1}.(keepfield{j});
    keep = true;
  else
    keep = false;
  end
  for i=1:numel(varargin)
    if ~isfield(varargin{i}, keepfield{j}) || ~isequal(varargin{i}.(keepfield{j}), tmp)
      keep = false;
      break
    end
  end
  if keep
    data.(keepfield{j}) = tmp;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   varargin
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data


function assign_var(var, val)
% Note: using an anonymous function as follows does not work in Octave:
%
% **    assign_var = @(var, val) assignin('caller', var, val);
%
% Also using the name 'assign' does not seem to work, hence 'assign_var'

   assignin('caller', var, val);


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS for doing math on each element of a cell-array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = cellplus(x, y)
if ~iscell(y)
  y = repmat({y}, size(x));
end
z = cellfun(@plus, x, y, 'UniformOutput', false);

function z = cellminus(x, y)
if ~iscell(y)
  y = repmat({y}, size(x));
end
z = cellfun(@minus, x, y, 'UniformOutput', false);

function z = celltimes(x, y)
if ~iscell(y)
  y = repmat({y}, size(x));
end
z = cellfun(@times, x, y, 'UniformOutput', false);

function z = cellrdivide(x, y)
if ~iscell(y)
  y = repmat({y}, size(x));
end
z = cellfun(@rdivide, x, y, 'UniformOutput', false);

function z = celllog10(x)
z = cellfun(@log10, x, 'UniformOutput', false);

function z = cellabs(x)
z = cellfun(@abs, x, 'UniformOutput', false);
