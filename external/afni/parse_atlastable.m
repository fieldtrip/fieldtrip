function T = parse_atlastable(atlastable)

% function to convert a string (read from an afni .HEAD file) with
% atlastable information into a struct

% assume that the string is well-behaved
ap = regexp(atlastable, '<ATLAS_POINT');
st = regexp(atlastable, 'STRUCT');
v  = regexp(atlastable, 'VAL');
o  = regexp(atlastable, 'OKEY');
g  = regexp(atlastable, 'GYoAR');

assert(numel(ap)==numel(st));
assert(numel(ap)==numel(v));
assert(numel(ap)==numel(o));
assert(numel(ap)==numel(g));

ap(end+1) = numel(atlastable);

T = [];
for k = 1:numel(st)
  T(k).struct = getcontent(atlastable(st(k):(v(k)-1)));
  T(k).val    = str2double(getcontent(atlastable(v(k):(o(k)-1))));
  T(k).okey   = str2double(getcontent(atlastable(o(k):(g(k)-1))));
  T(k).gyoar  = str2double(getcontent(atlastable(g(k):(ap(k+1)-1))));
end

function out = getcontent(substring)

sel = regexp(substring, '"');
out = substring((sel(1)+1):(sel(2)-1));
