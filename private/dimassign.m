function M=dimassign(A,dim,idx,B);
% function M=dimassign(A,dim,idx,B);
% 
% The purpose of the function is shown by the following example:
% If A and B are multidimensional matrixes,
% A=dimassign(A,4,23,B); is the same as A(:,:,:,23,:,:,...)=B;
% The difference is that the dimention is selected by a scalar, not by
% the place between the brackets.
% A(2,4,3)=B; will then be written as: A=dimassign(A,[1,2,3],[2,4,3],B);
% In this last case B, of cource, must be a scalar.
% A([1,2],:,3)=B; can be written as: A=dimassign(A,[1,3],{[1,2],3},B);
% Of cource, again, the dimensions of B must fit! 
% (size(B)==size(A([1,2],:,3) in this particular case)
%
% See also the function DIMINDEX

% Copyright (C) 2005, Geerten Kramer
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

if(~iscell(idx))
	if(~any(size(dim)==1)||~any(size(idx)==1)||ndims(dim)>2||ndims(idx)>2||...
		length(dim)~=length(idx))
		error('dim and idx must be both scalars oor both must have the same length');
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

if(~isequal(unique(dim),sort(dim)))
	error('dim must be unique, every dimention can be addressed only once');
end;

Na=ndims(A);
for(i=1:max([max(dim),Na]))
	ref=find(dim==i);
	if(isempty(ref))
		C{i}=':';
	else
		C{i}=idx{ref};
	end;
end;
M=A;		
try
	M(C{:})=B;
catch
	error('Subscripted assignment dimension mismatch.');
end;
