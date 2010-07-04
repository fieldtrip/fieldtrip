function split = split_string(in, sep, varargin)
% splits a string based on a seperator string
%
% P=SPLIT_STRING(S,D) splits a string S based on delimiter D, and
% returns a cell P with the different parts of string S. Both S and D
% should have only one row.
%
% If P contains N occurences of delimeter D, then P contains N strings so
% that strcat(P{1}, D, P{2}, D, ...,D, P{n}) == S.
%
% Extended syntax (Jan 2010):
% 
% Pi=SPLIT_STRING(S,D,I) returns the I-th element from S, that is 
% P{i} if used with the syntax above
%
% S=SPLIT_STRING(S,D1,I1,D2,I2,...) yields the same as
% SPLIT_STRING(SPLIT_STRING(S,D1,I1),D2,I2,...).
%
% Examples: 
%   split_string('a,b;c,d;e,f',';')            returns {'a,b','c,d','e,f'} 
%   split_string('a,b;c,d;e,f',';',2)          returns 'c,d'
%   split_string('a,b;c,d;e,f',';',2,',')      returns {'c','d'}
%
% NNO Nov 2009

if ~ischar(in) || ~ischar(sep)
    error('Input should be strings')
end

if size(in,1) ~= 1 || size(sep,1) ~= 1
    error('Input should have only one row')
end

idxs=findstr(in,sep);
startidxs=[1 idxs+length(sep)];
endidxs=[idxs-1 length(in)];

count=length(idxs)+1;
split=cell(1,count);
for i=1:length(idxs)+1
    split{i}=in(startidxs(i):endidxs(i));
end;

% process further arguments, if given
if numel(varargin)>0
    in=split{varargin{1}};
    if numel(varargin)>1
        % recursion
        split=split_string(in,varargin{2:end});
    else
        split=in;
    end
end