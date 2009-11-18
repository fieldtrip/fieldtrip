function y = cellsum(x,dim)
% CELLSUM sums a cell array over dimension dim where the elements of the
% cell array must be of equal form
%
% y = cellsum(x,dim)
%
% Copyright (C) 2007, Marcel van Gerven
%
% $Log: cellsum.m,v $
% Revision 1.1.1.1  2008/02/27 14:42:54  roboos
% Marcel van Gerven, version 27 Feb 2008
%
  
szx = size(x);

% we create an index
index = cell(1,length(szx));
for i=1:length(index)
    index{i} = ':';
end

% construct a new cell structure
szx(dim) = 1; y = cell(szx);

% fill new cell with zeros
szc = size(x{1});
for j=1:numel(y)
   y{j} = zeros(szc); 
end

for j=1:size(x,dim)

    % construct a reduced cell
    index{dim} = j; rc = x(index{:});
    
    % add reduced cell to new cell
    for k=1:numel(y)
        y{k} = y{k} + rc{k};
    end
    
end


