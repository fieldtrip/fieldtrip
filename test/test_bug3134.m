function test_bug3134

% WALLTIME 00:10:00
% MEM 1gb


list = {
  'isvector'
  'ismatrix'
  'iscolumn'
  'isrow'
  'isreal'
};

% These should not generally be present. They were part of the system identification toolbox from 2007a up to 2011b, but were removed in 2012a.
%  'isrealvec'
%  'isrealmat'

for i=1:numel(list)
  assert(~isempty(which(list{i})), sprintf('%s is missing\n', list{i}));
end
