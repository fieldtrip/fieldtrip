function s=unsplit_string(c, sep)
% 'zips' a cell with strings and a seperator symbol; the inverse of
% SPLIT_STRING
%
% S=UNSPLIT_STRING(C,SEP) takes a cell with strings C and a seperator
% string SEP, and returns the string [C{1} SEP C{2} SEP ... SEP C{end}].
%
% It is assumed that C and SEP are strings with exactly one row.
%
% NNO Jan 2010

if ~iscell(c)
    error('Expected cell input');
end

n=numel(c);
if n==0
    s='';
    return
end
    
r=cell(1,2*n);
for k=1:n
    r{k*2-1}=sep;
    r{k*2}=c{k};
end

s=[r{2:end}];




