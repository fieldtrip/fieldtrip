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

if ~iscellstr(c)
    error('Expected cell string input');
end

n=numel(c);
if n==0
    s='';
    return
end

% insert separator string
parts=cell(1,2*n-1);
for k=1:n
    parts{k*2-1}=c{k};
    if k<n
        parts{k*2}=sep;
    end
end

% concatenate
s=sprintf('%s',parts{:});
