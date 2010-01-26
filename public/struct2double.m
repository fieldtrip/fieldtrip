function [x] = struct2double(x, maxdepth);

% STRUCT2DOUBLE converts all single precision numeric data in a structure 
% into double precision. It will also convert plain matrices and
% cell-arrays.
%
% Use as
%    x = struct2double(x);
%
% Starting from Matlab 7.0, you can use single precision data in your
% computations, i.e. you do not have to convert back to double precision.
%
% Matlab version 6.5 and older only support single precision for storing
% data in memory or on disk, but do not allow computations on single
% precision data. Therefore you should converted your data from single to
% double precision after reading from file.
%
% See also STRUCT2SINGLE

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: struct2double.m,v $
% Revision 1.5  2009/02/25 09:25:55  roboos
% also deal with intXX
%
% Revision 1.4  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.3  2005/07/18 14:35:04  roboos
% added support for structure arrays (caused error)
%
% Revision 1.2  2005/06/30 13:56:49  roboos
% corrected the function name inside the file (thanks to Christian)
%
% Revision 1.1  2005/06/29 09:50:31  roboos
% new implementation, intended to reduce the size of Matlab files on disk
%

fieldtripdefs

if nargin<2
  maxdepth = inf;
end

% convert the data, work recursively through the complete structure
x = convert(x, 0, maxdepth);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this subfunction does the actual work
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a] = convert(a, depth, maxdepth);

if depth>maxdepth
  error('recursive depth exceeded');
end

switch class(a)
  case 'struct'
    % process all fields of the structure recursively
    fna = fieldnames(a);
    % process all elements of the array
    for j=1:length(a(:))
      % warning, this is a recursive call to traverse nested structures
      for i=1:length(fna)
        fn = fna{i};
        ra = getfield(a(j), fn);
        ra = convert(ra, depth+1, maxdepth);
        a(j) = setfield(a(j), fn, ra);
      end
    end 

  case 'cell'
    % process all elements of the cell-array recursively
    % warning, this is a recursive call to traverse nested structures
    for i=1:length(a(:))
      a{i} = convert(a{i}, depth+1, maxdepth);
    end

  case {'single' 'int32' 'uint32' 'int16' 'uint16'}
    % convert the values to double precision
    a = double(a);

  case 'double'
    % keep as it is

  otherwise
    warning('not converting class %s', class(a))
    % do nothing
end
