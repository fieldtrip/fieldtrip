function y = config(x, varargin)

% CONFIG constructs an object that can be used to track used and unused
% options that are represented in the cfg structure in the FieldTrip
% toolbox.
%
% y = CONFIG( ) creates an empty config object.
% y = CONFIG(x) creates a config object from the
%         existing config object x, or from the existing
%         structure x.
% y = CONFIG('field1', '..', 'field2', '..' ) creates a
%         config object with the specified fields set to
%         the specified values.
% The assignment and reference counters are initialized to zero.
%
% A config object is a structure-like object, but with the
% additional feature that the number of assignments and references
% to each field are being counted.
%
% The config object can be used throughout your FieldTrip function
% as if it were a plain MATLAB structure. At the end of your FieldTrip
% function you can report the fields that were
%   - not used (references=0)
%   - locally assigned as compared to assigned by the end user (assignments>0)
% Subsequently you are able to the fields that are locally assigned but not
% being used (e.g. irrelevant default values).
%
% As an example, construct a config object with two fields
%   A = config('a', 1, 'b', 2);
%   A.a = 1;
%   A.a = 1;
%   A.a = 1;
%   access(A)   % this shows the values, the number of assignments and references
%   disp(A.b);
%   disp(A.b);
%   access(A)   % this shows the values, the number of assignments and references
%
% A config object can also contain structures
%   A = config('a', 1, 'b', 2, 'c', struct('a', 1));
%
% A config object can even contain nested config objects
%   A = config('a', 1, 'b', 2, 'c', config('a', 1))
%
% The nesting of config objects can be arbitrary deep
%   a = 1
%   b = config('a', a)
%   c = config('b', b)
%   d = config('c', c)

% Copyright (C) 2012-2015, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

if nargin==1
  if isa(x, 'config')
    % formally no conversion is needed, but a copy will be made in which the counters are reset
    y = deepcopy(x);
    key = fieldnames(x);
    for i=1:length(key)
      setzero(y.assign   .(key{i}));      % reset the counter to zero
      setzero(y.reference.(key{i}));      % reset the counter to zero
      setzero(y.original .(key{i}));      % first set to zero and then increment with one,
      increment(y.original.(key{i}));   % since all fields were present in the original
    end
  elseif isa(x, 'struct')
    % convert the input structure into a config object
    key = fieldnames(x);
    for j=1:numel(x)
      val = {};
      for i=1:length(key)
        try
          val{i} = x(j).(key{i});
        catch
          val{i} = [];
        end
        % use recursion to let some other part of the code handle the remainder
        if isa(val{i}, 'struct')
          val{i} = config(val{i});
        end
      end
      tmp           = struct();
      tmp.value     = struct();
      tmp.assign    = struct();
      tmp.reference = struct();
      tmp.original  = struct();
      tmp.hidden    = struct(); % this contains hidden fields which are not tracked
      for i=1:length(key)
        tmp.value.(key{i})     = val{i};
        tmp.assign.(key{i})    = deepcopy(0); % ensure that a unique scalar is created for each counter
        tmp.reference.(key{i}) = deepcopy(0); % ensure that a unique scalar is created for each counter
        tmp.original.(key{i})  = deepcopy(1); % ensure that a unique scalar is created for each counter
      end
      y(j) = class(tmp,'config');
    end
    if numel(x)
      y = reshape(y, size(x));
    else
      y = config;
    end
  else
    error('Unsupported input class ''%s'' for constructing a config object', class(x));
  end

elseif nargin>1
  if mod(nargin,2)
    error('Incorrect number of input arguments (should be key-value pairs)')
  end
  varargin = {x varargin{:}};
  key = varargin(1:2:end);
  val = varargin(2:2:end);

  % When having y.assign and y.reference point to the same scalars, there is
  % a side effect in the increment function that reveals that the scalars
  % representing the different counters all point to the same physical
  % memory address. Therefore I have to ensure that there is a unique
  % scalar for each individual counter.
  assign    = {};
  reference = {};
  original  = {};
  for i=1:length(key)
    assign   {end+1} = key{i};
    reference{end+1} = key{i};
    original {end+1} = key{i};
    assign   {end+1} = deepcopy(0);  % ensure that a unique scalar is created for each counter
    reference{end+1} = deepcopy(0);  % ensure that a unique scalar is created for each counter
    original {end+1} = deepcopy(1);  % ensure that a unique scalar is created for each counter
  end

  for i=1:length(val)
    if isa(val{i}, 'struct')
      % use recursion to convert sub-structures into sub-configs
      val{i} = config(val{i});
    end
  end

  y.value     = struct(varargin{:});
  y.assign    = struct(assign{:});
  y.reference = struct(reference{:});
  y.original  = struct(original{:});
  y.hidden    = struct();  % this contains hidden fields which are not tracked
  y = class(y,'config');

else
  % create an empty config object
  y           = struct();
  y.value     = struct();
  y.assign    = struct();
  y.reference = struct();
  y.original  = struct();
  y.hidden    = struct(); % this contains hidden fields which are not tracked
  y = class(y,'config');

end

