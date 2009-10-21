function [val] = keyval(key, varargin)

% KEYVAL returns the value that corresponds to the requested key in a
% key-value pair list of variable input arguments
%
% Use as
%   [val] = keyval(key, varargin)
%
% See also VARARGIN

% Copyright (C) 2005-2007, Robert Oostenveld
%
% $Log: keyval.m,v $
% Revision 1.5  2009/08/04 11:58:28  roboos
% perform a case-insensitive string comparison for keys
%
% Revision 1.4  2009/07/14 16:11:02  roboos
% speed up the input checks
%
% Revision 1.3  2009/05/14 19:24:02  roboos
% removed ; at end of function declaration
%
% Revision 1.2  2009/01/06 09:05:26  roboos
% added additional check on optional input arguments: the 1st, 3rd, etc. should be strings (keys)
%
% Revision 1.1  2008/11/13 09:55:36  roboos
% moved from fieldtrip/private, fileio or from roboos/misc to new location at fieldtrip/public
%
% Revision 1.2  2007/07/18 12:43:53  roboos
% test for an even number of optional input arguments
%
% Revision 1.1  2005/11/04 10:24:46  roboos
% new implementation
%

if length(varargin)==1 && iscell(varargin{1})
  varargin = varargin{1};
end

if mod(length(varargin),2)
  error('optional input arguments should come in key-value pairs, i.e. there should be an even number');
end

% the 1st, 3rd, etc. contain the keys, the 2nd, 4th, etc. contain the values
keys = varargin(1:2:end);
vals = varargin(2:2:end);

if ~all(cellfun(@ischar, keys))
  error('optional input arguments should come in key-value pairs, the optional input argument %d is invalid (should be a string)', i);
end

hit = find(strcmpi(key, keys));
if isempty(hit)
  % the requested key was not found
  val = [];
elseif length(hit)==1  
  % the requested key was found
  val = vals{hit};
else
  error('multiple input arguments with the same name');
end

