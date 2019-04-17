function y=spcharin(x)

%SPCHARIN replaces special characters by their codes
%
% Syntax y=spcharin(x)
%
% Description
%   x is a character or 2D cell array of characters
%   y is the corresponding version with reserved (special) characters 
%   replaced by '%ascii;' codes
%
% See also: SPCHAROUT
%
% Jonas Almeida, almeidaj@musc.edu, 8 Oct 2002, MAT4NAT Tbox

if iscell(x)
    [n,m]=size(x);
    for i=1:n
        for j=1:m
            y(i,j)={spcharin(x{i,j})};
        end
    end   
elseif ischar(x)
    x=['AAA',x,'AAA']; % add polyA header and tail
    % replace delimiters first (% and ;)
    x=strrep(x,';','*59;');
    x=strrep(x,'#','#35;');
    x=strrep(x,'*59;','#59;');
    ascii=x*1; % Find special characters
    sp=find(~(((x>47)&(x<58))|((x>96)&(x<123))|((x>64)&(x<91))|(x==59)|(x==35)));
    % Replace them by ascii code delimited by % and ;
    for i=length(sp):-1:1
        x=[x(1:(sp(i)-1)),'#',num2str(ascii(sp(i))),';',x((sp(i)+1):end)];
    end
    y=x(4:end-3); % note plyA head and tail are removed
    %y=x;
else
    w=whos('x')
    error(['string expected, ',w.class',' found instead'])
end