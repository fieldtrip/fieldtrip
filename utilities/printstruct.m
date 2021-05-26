function str = printstruct(name, val, varargin)

% PRINTSTRUCT converts a MATLAB structure into a multiple-line string that, when
% evaluated by MATLAB, results in the original structure. It also works for most
% other standard MATLAB classes, such as numbers, vectors, matrices, and cell-arrays.
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
%
% See also DISP

% Copyright (C) 2006-2018, Robert Oostenveld
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
  val  = name;
  name = inputname(1);
end

if isa(val, 'config')
  % This is FieldTrip specific: the @config object resembles a structure but tracks
  % the access to each field. In this case it is to be treated as a normal structure.
  val = struct(val);
end

% get the options and set defaults
transposed    = ft_getopt(varargin, 'transposed', false);     % print a column as a transposed row, or a row as a transposed column
linebreaks    = ft_getopt(varargin, 'linebreaks', true);      % print with or without linebreaks
lastnewline   = ft_getopt(varargin, 'lastnewline', false);    % end with a newline

% when called recursively, it should keep the last newline intact
varargin = ft_setopt(varargin, 'lastnewline', true);

% Note that because we don't know the final size of the string, iteratively appending
% is actually faster than creating a cell-array and subsequently doing a cat(2,
% strings{:}). Note also that sprintf() is slow.
str = '';

% Note further that in the string concatenations I use the numerical value
% of a newline (\n), which is 10.

if numel(val) == 0
  if iscell(val)
    str = [name ' = {};' 10];
  elseif isstruct(val)
    str = [name ' = [];' 10];
  elseif isnumeric(val)
    str = [name ' = [];' 10];
  end
  
elseif isstruct(val)
  if numel(val)>1
    % print it as a struct-array
    str = cell(size(val));
    for i=1:numel(val)
      str{i} = printstruct(sprintf('%s(%d)', name, i), val(i), varargin{:});
    end
    str = cat(2, str{:});
  elseif numel(val)==1
    % print it as a named structure
    fn = fieldnames(val);
    for i=1:length(fn)
      fv = val.(fn{i});
      switch class(fv)
        case {'char' 'string'}
          line = printstr([name '.' fn{i}], fv, linebreaks, transposed);
        case {'single' 'double' 'int8' 'int16' 'int32' 'int64' 'uint8' 'uint16' 'uint32' 'uint64' 'logical'}
          if ismatrix(fv)
            line = printmat([name '.' fn{i}], fv, linebreaks, transposed);
          else
            line = '''ERROR: multidimensional arrays are not supported''';
          end
        case 'cell'
          line = printcell([name '.' fn{i}], fv, linebreaks, transposed);
        case 'struct'
          line = [printstruct([name '.' fn{i}], fv, varargin{:}) 10];
        case 'function_handle'
          line = printstr([name '.' fn{i}], func2str(fv), linebreaks, transposed);
        otherwise
          ft_error('unsupported');
      end
      str = [str line newline];
    end % for fieldnames
    if isempty(fn)
      % print it as an empty structure
      str = [name ' = struct();'];
    end
  end
  
elseif ~isstruct(val)
  % print it as a named variable
  switch class(val)
    case {'char' 'string'}
      str = printstr(name, val, linebreaks, transposed);
    case {'double' 'single' 'int8' 'int16' 'int32' 'int64' 'uint8' 'uint16' 'uint32' 'uint64' 'logical'}
      str = printmat(name, val, linebreaks, transposed);
    case 'cell'
      str = printcell(name, val, linebreaks, transposed);
    otherwise
      ft_error('unsupported');
  end
end

% replace double newlines with a single one
str = strrep(str, [newline newline], newline);

if isempty(name)
  % remove the ' = ' part from the string
  str = str(4:end);
end

if str(end)==newline && ~lastnewline
  % remove the last newline
  str = str(1:end-1);
end

if ~nargout
  disp(str)
  clear str
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = printstr(name, val, linebreaks, transposed)
siz = size(val);
if siz(1)>1
  str = sprintf('%s = \n', name);
  for i=1:siz(1)
    str = [str sprintf('  %s\n', printval(val(i,:), linebreaks, transposed))];
  end
elseif siz(1)==1
  str = sprintf('%s = %s;\n', name, printval(val, linebreaks, transposed));
else
  str = sprintf('%s = '''';\n', name);
end
if ~linebreaks
  str(str==10) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = printcell(name, val, linebreaks, transposed)
siz = size(val);
if isempty(val)
  str = sprintf('%s = {};\n', name);
elseif numel(siz) == 2 && siz(2) == 1 && transposed
  % print column vector as (non-conjugate) transposed row
  str = printcell(name, transpose(val), linebreaks, false);
  str = sprintf('%s'';\n', str(1:end-2));
elseif numel(siz) == 2 && siz(1) == 1 && transposed
  % print row vector as (non-conjugate) transposed column
  str = printcell(name, transpose(val), linebreaks, false);
  str = sprintf('%s'';\n', str(1:end-2));
else
  str = sprintf('%s = {\n', name);
  for i=1:siz(1)
    for j=1:(siz(2)-1)
      dum = printval(val{i,j}, linebreaks, transposed);
      if endsWith(dum, ';')
        dum = dum(1:end-1);
      elseif endsWith(dum, [';' newline])
        dum = dum(1:end-2);
      end
      str = [str ' ' dum ',']; % add the element with a comma
    end
    dum = printval(val{i,siz(2)}, linebreaks, transposed);
    if endsWith(dum, ';')
      dum = dum(1:end-1);
    elseif endsWith(dum, [';' newline])
      dum = dum(1:end-2);
    end
    str = [str ' ' dum 10];  % add the last element of each row with a linebreak
  end
  if siz(1)==1
    str(str==10) = [];
  end
  str = sprintf('%s };\n', str);
end
if ~linebreaks
  str(str==10) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = printmat(name, val, linebreaks, transposed)
siz = size(val);
if numel(val) == 0
  str = sprintf('%s = [];\n', name);
elseif numel(val) == 1
  % an integer will never get trailing decimals when using %g
  str = sprintf('%s = %g;', name, val);
elseif numel(siz) == 2 && siz(2) == 1 && transposed
  % print column vector as (non-conjugate) transposed row
  str = printmat(name, transpose(val), linebreaks, false);
  str = sprintf('%s''', str);
elseif numel(siz) == 2 && siz(1) == 1 && transposed
  % print row vector as (non-conjugate) transposed column
  str = printmat(name, transpose(val), linebreaks, false);
  str = sprintf('%s''', str);
elseif ismatrix(val)
  if isa(val, 'double')
    str = sprintf('%s = %s;', name, mat2str(val));
  else
    % add class information for non-double numeric matrices
    str = sprintf('%s = %s;', name, mat2str(val, 'class'));
  end
  % add a newline after each semicolon
  str = strrep(str, ';', [';' newline]);
else
  ft_warning('multidimensional arrays are not supported');
  str = '''ERROR: multidimensional arrays are not supported''';
end
if ~linebreaks
  str(str==10) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = printval(val, linebreaks, transposed)

switch class(val)
  case 'char'
    str = ['''' val ''''];
    
  case 'string'
    % these use " in the declaration rather than '
    str = ['"' char(val) '";'];
    
  case {'single' 'double' 'int8' 'int16' 'int32' 'int64' 'uint8' 'uint16' 'uint32' 'uint64' 'logical'}
    str = printmat('', val, linebreaks, transposed);
    str = str(4:end);
    if endsWith(str, ';')
      str = str(1:end-1);
    end
    
  case 'function_handle'
    str = ['@' func2str(val)];
    
  case 'struct'
    % print it as an anonymous structure
    str = 'struct(';
    fn = fieldnames(val);
    for i=1:numel(fn)
      str = [str '''' fn{i} '''' ', ' printval(val.(fn{i}), linebreaks, transposed)];
    end
    str = [str ')'];
    
  otherwise
    ft_warning('cannot print unknown object at this level');
    str = '''ERROR: cannot print unknown object at this level''';
end
if ~linebreaks
  str(str==10) = [];
end
