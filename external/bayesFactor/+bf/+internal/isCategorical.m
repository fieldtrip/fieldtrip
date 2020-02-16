function value = isCategorical(tbl,name)  
% Determine which columns in a table contain categorical variables.
% 
% tbl  - The table
% name - Optional. Specify a specific column or cell array of columns to
%           check. 
% OUTPUT
% value - logical matching the columns in tbl or the name cell-array
isCatColumn = @(x) (isa(x,'categorical') || iscellstr(x) || isstring(x) || ischar(x) || islogical(x));
value = logical(varfun(isCatColumn,tbl,'OutputFormat','Uniform'));
if nargin>1
    stay = ismember(tbl.Properties.VariableNames,name);
    value= value(stay);
end