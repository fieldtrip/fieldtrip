function [input, output] = funargname(fname)

% FUNARGNAME returns the input and output arguments of the function
% by parsing the m-file
%
% Use as
%    [input, output] = funargname(fname)
% where the input and output function arguments will be returned
% as cell-arrays containing strings.

fid = fopen(which(fname), 'r');
line = '';
while isempty(regexp(line, '^function', 'once')) && ~feof(fid)
  line = fgetl(fid);
end
if isempty(regexp(line, '^function', 'once'))
  error('could not locate function declaration in %s', fname);
end

% the declaration of the function should look like this
% function [out1, out2] = funname(in1, in2); % optional comment

% clean up the line
line(line=='[') = ' ';
line(line==']') = ' ';
line(line=='(') = ' ';
line(line==')') = ' ';
line(line==',') = ' ';
line(line==';') = ' ';
if any(line=='%')
  sel = find(line=='%', 'first');
  line(sel:end) = [];
end
line = strtrim(line);
% the line will now look like
% function out1 out2 funname in1 in2

% cut it into pieces and return the relevant ones
line = tokenize(line, ' ', true);
f1 = find(strcmp(line, '=')); % the index of =
f2 = f1+1;                    % the index of the function name
output = line(2:(f1-1));
input  = line((f2+1):end);
