function y=spcharout(x)

%SPCHAROUT replaces special character codes by their values
%
%Syntax y=spcharout(x)
%
%Description
%   x is a character or 2D cell array of characters
%   y is the corresponding version with special character codes 
%     '%ascii;' replaced by their values
%
%See also: spcharin
%
% Jonas Almeida, almeidaj@musc.edu 8 Oct 2002, MAT4NAT Tbox

if iscell(x)
    [n,m]=size(x);
    for i=1:n
        for j=1:m
            y(i,j)={spcharout(x{i,j})};
        end
    end   
elseif ischar(x)
    y=eval(['[''',regexprep(x,'\#(\d+);',''',char($1),'''),''']']);
else
    w=whos('x')
    error(['string expected, ',w.class',' found instead'])
end