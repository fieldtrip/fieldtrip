function unit = estimate_units(size)

% ESTIMATE_UNITS tries to determine the units of a geometrical object by
% looking at its size and by relating this to the size of the human
% brain.
%
% Use as
%   unit = estimate_units(size)
%
% This function will return one of the following strings
%   'm'
%   'dm'
%   'cm'
%   'mm'

% Copyright (C) 2009, Robert Oostenveld
%
% $Log: estimate_units.m,v $
% Revision 1.1  2009/03/26 14:57:04  roboos
% moved the unit estimation to a seperate function
%

% do some magic based on the size
unit = {'m', 'dm', 'cm', 'mm'};
indx = round(log10(size)+2-0.2);

if indx>length(unit)
  indx = length(unit);
  warning('assuming that the units are "%s"', unit{indx});
end

if indx<1
  indx = 1;
  warning('assuming that the units are "%s"', unit{indx});
end

unit = unit{indx};
