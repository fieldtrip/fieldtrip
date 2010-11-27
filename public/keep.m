function keep(varargin)

% KEEP clears all the variables in the workspace except the ones you specify
% after the keep command. It works just like clear, but only for variables
%
% KEEP is a modified version of the KEEP M-file written by 
% Xiaoning (David) Yang (1998). It allows the use of the wildcard * as in the
% clear command
%
% Examples:
%   keep set1data set2data set3data
%   keep *data
%
% Martin Barugel (mbarugel@utdt.edu)

%       Keep all
if isempty(varargin) | sum(strcmp('*',varargin))>=1
        return
end

%       See what are in caller workspace
wh = evalin('caller','who');

del=' ';
for i=1:length(wh)
    del=[del wh{i} '  '];
end
        
%       Check workspace variables
if isempty(wh)
        error('  There is nothing to keep!')
end

%       Create string of variables to be cleared
for i=1:length(varargin)
    name=varargin{i};
    name=strrep(name,'*','\w*');
    [s,f]=regexp(del,['\s' name '\s']);
    for j=1:length(s)
        del(s(j):f(j))=char(32*ones(1,f(j)-s(j)+1));
    end
end

%       Clear them
if length(del)==sum(isspace(del))
    return
else
    evalin('caller',['clear ' del]);
end
