function keyval = struct2keyval(s)
% STRUCT2KEYVAL converts a structure to a key-value pair cell array, e.g.
% the following struct:
%
% s = 
%   struct with fields:
%     a: 1
%     b: 2
%
% will be converted into the following:
%
% >> struct2keyval(s)
% ans =
%   1×4 cell array
%     {'a'}    {[1]}    {'b'}    {[2]}

s = struct(s);
keyval = [fieldnames(s).'; struct2cell(s).'];
keyval = keyval(:).';

end