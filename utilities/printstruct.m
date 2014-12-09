function str = printstruct(name, val)

% PRINTSTRUCT converts a MATLAB structure into a multi-line string that can be
% interpreted by MATLAB, resulting in the original structure.
%
% Use as
%   str = printstruct(val)
% or
%   str = printstruct(name, val)
% where "val" is any MATLAB variable, e.g. a scalar, vector, matrix, structure, or
% cell-array. If you pass the name of the variable, the output is a piece of MATLAB code
% that you can execute, i.e. an ASCII serialized representation of the variable.
%
% Example
%   a.field1 = 1;
%   a.field2 = 2;
%   s = printstruct(a)
%
%   b = rand(3);
%   s = printstruct(b)
%
%   s = printstruct('c', randn(10)>0.5)

% Copyright (C) 2006-2013, Robert Oostenveld
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

if nargin==1
  val  = name;
  name = inputname(1);
end

if isa(val, 'config')
  % this is fieldtrip specific: the @config object resembles a structure but tracks the
  % access to each field.  In this case it is to be treated as a normal structure.
  val = struct(val);
end

% note here that because we don't know the final size of the string,
% iteratively appending is actually faster than creating a cell array and
% subsequently doing a cat(2, strings{:}) (also sprintf() is slow)
str = '';

% note further that in the string concatenations I use the numerical value
% of a newline (\n), which is 10

if numel(val) == 0
  if iscell(val)
    str = [name ' = {};' 10];
  else
    str = [name ' = [];' 10];
  end
elseif isstruct(val)
  if numel(val)>1
    str = cell(size(val));
    for i=1:numel(val)
      str{i} = printstruct(sprintf('%s(%d)', name, i), val(i));
    end
    str = cat(2, str{:});
    return
  else
    % print it as a named structure
    fn = fieldnames(val);
    for i=1:length(fn)
      fv = val.(fn{i});
      switch class(fv)
        case 'char'
          line = printstr([name '.' fn{i}], fv);
        case {'single' 'double' 'int8' 'int16' 'int32' 'int64' 'uint8' 'uint16' 'uint32' 'uint64' 'logical'}
          if ismatrix(fv)
            line = [name '.' fn{i} ' = ' printmat(fv) ';' 10];
          else
            line = '''FIXME: printing multidimensional arrays is not supported''';
          end
        case 'cell'
          line = printcell([name '.' fn{i}], fv);
        case 'struct'
          line = [printstruct([name '.' fn{i}], fv) 10];
        case 'function_handle'
          line = printstr([name '.' fn{i}], func2str(fv));
        otherwise
          error('unsupported');
      end
      str  = [str line];
    end
  end
elseif ~isstruct(val)
  % print it as a named variable
  switch class(val)
    case 'char'
      str = printstr(name, val);
    case {'double' 'single' 'int8' 'int16' 'int32' 'int64' 'uint8' 'uint16' 'uint32' 'uint64' 'logical'}
      str = [name ' = ' printmat(val)];
    case 'cell'
      str = printcell(name, val);
    otherwise
      error('unsupported');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = printcell(name, val)
siz = size(val);
if isempty(val)
  str = sprintf('%s = {};\n', name);
  return;
end
if all(size(val)==1)
  str = sprintf('%s = { %s };\n', name, printval(val{1}));
else
  str = sprintf('%s = {\n', name);
  for i=1:siz(1)
    dum = '';
    for j=1:(siz(2)-1)
      dum = [dum ' ' printval(val{i,j}) ',']; % add the element with a comma
    end
    dum = [dum ' ' printval(val{i,siz(2)})]; % add the last one without comma
    
    str = [str dum 10];
  end
  str = sprintf('%s};\n', str);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = printstr(name, val)
siz = size(val);
if siz(1)>1
  str = sprintf('%s = \n', name);
  for i=1:siz(1)
    str = [str sprintf('  %s\n', printval(val(i,:)))];
  end
elseif siz(1)==1
  str = sprintf('%s = %s;\n', name, printval(val));
else
  str = sprintf('%s = '''';\n', name);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = printmat(val)
if numel(val) == 0
  str = '[]';
elseif numel(val) == 1
  % an integer will never get trailing decimals when using %g
  str = sprintf('%g', val);
elseif ismatrix(val)
  if isa(val, 'double')
    str = mat2str(val);
  else
    % add class information for non-double numeric matrices
    str = mat2str(val, 'class');
  end
  str = strrep(str, ';', [';' 10]);
else
  warning('multidimensional arrays are not supported');
  str = '''FIXME: printing multidimensional arrays is not supported''';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = printval(val)

switch class(val)
  case 'char'
    str = ['''' val ''''];
    
  case {'single' 'double' 'int8' 'int16' 'int32' 'int64' 'uint8' 'uint16' 'uint32' 'uint64' 'logical'}
    str = printmat(val);
    
  case 'function_handle'
    str = ['@' func2str(val)];
    
  case 'struct'
    warning('cannot print structure at this level');
    str = '''FIXME: printing structures at this level is not supported''';
    
  otherwise
    warning('cannot print unknown object at this level');
    str = '''FIXME: printing unknown objects is not supported''';
end
