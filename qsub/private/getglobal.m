function var = getglobal

% GETGLOBAL gets all global variables and puts them in a structure
%
% Use as
%   var = getglobal;
%   setglobal(var);
%
% See also SETGLOBAL

persistent previous_list previous_varout

list = whos('global');

if isequal(list, previous_list)
    var = previous_varout;
    return;
end
    
var = [];
for i=1:length(list)
  eval(sprintf('global %s', list(i).name));
  var.(list(i).name) = eval(list(i).name);
end

previous_list = list;
previous_varout = var;
