function M=dimindex(A,dim,idx);
% function M=dimindex(A,dim,idx);
% 
% The purpose of the function is shown by the following example:
% If A is a multidimensional matrix,
% dimindex(A,4,23); is the same as A(:,:,:,23,:,:,...);
% The advantage is that the dimention is selected by a scalar, not by
% the place between the brackets.
% A(2,4,3); will then be written as: dimindex(A,[1,2,3],[2,4,3]);
% A(4,:,[5:10]); will then be written as: dimindex(A,[1,3],{4,[5:10]});
%
% See also the function DIMASSIGN

% Copyright (C) 2005, Geerten Kramer
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

if(~iscell(idx))
	if(~any(size(dim)==1)||~any(size(idx)==1)||ndims(dim)>2||ndims(idx)>2||...
		length(dim)~=length(idx))
		error('dim and idx must be both scalars or both vectors of the same size');
	end;
	dummi=[];
	for(i=1:length(idx))
		dummi{i}=idx(i);
	end;
	idx=dummi;
	clear dummi;
end;
if(~any(size(dim)==1)||~any(size(idx)==1)||ndims(dim)>2||ndims(idx)>2||...
	length(dim)~=length(idx))
	error('dim and idx must be both scalars or both must have the same length');
end;


if(length(dim)>ndims(A)||any(dim>ndims(A)))
	error('dim, or one of its contents are larger than the number of dimentions in A');
end;
if(~isequal(unique(dim),sort(dim)))
	error('dim must be unique, every dimention can be addressed only once');
end;

N=ndims(A);
for(i=1:N)
	ref=find(dim==i);
	if(isempty(ref))
		C{i}=':';
	else
		C{i}=idx{ref};
	end;
end;
M=A(C{:});
	
