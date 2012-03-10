function str = fixname(str)

% FIXNAME changes all inappropriate characters in a sting into '_' 
% such that it can be used as a filename or as a structure field name.
%
% Use as
%   str = fixname(str)
%
% See also DEBLANK

% Copyright (C) 2012, Robert Oostenveld

str = lower(str);
str(str=='-') = '_'; % fix dashes
str(str==' ') = '_'; % fix spaces
str(str=='/') = '_'; % fix forward slashes
str(str=='\') = '_'; % fix backward slashes
str(str=='!') = '_';
str(str=='@') = '_';
str(str=='#') = '_';
str(str=='$') = '_';
str(str=='%') = '_';
str(str=='^') = '_';
str(str=='&') = '_';
str(str=='*') = '_';
str(str=='(') = '_';
str(str==')') = '_';
str(str=='{') = '_';
str(str=='}') = '_';
str(str=='[') = '_';
str(str==']') = '_';
str(str=='<') = '_';
str(str=='>') = '_';
str(str=='?') = '_';
