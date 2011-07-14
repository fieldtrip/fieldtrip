function setglobal(var)

% SETGLOBAL assigns global variables from a structure
%
% Use as
%   var = getglobal;
%   setglobal(var);
%
% See also GETGLOBAL

if isempty(var)
  return
end

list = fieldnames(var);

for i=1:length(list)
  eval(sprintf('global %s', list{i}));
  eval(sprintf('%s = var.%s;', list{i}, list{i}));
end
