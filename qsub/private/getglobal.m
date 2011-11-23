function var = getglobal

% GETGLOBAL gets all global variables and puts them in a structure
%
% Use as
%   var = getglobal;
%   setglobal(var);
%
% See also SETGLOBAL

list = whos('global');

var = [];
for i=1:length(list)
  eval(sprintf('global %s', list(i).name));
  var.(list(i).name) = eval(list(i).name);
end
