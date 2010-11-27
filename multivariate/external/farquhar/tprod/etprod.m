function [C,Atp,Btp]=etprod(Cidx,A,Aidx,B,Bidx)
% tprod wrapper to make calls more similar to Einstein Summation Convention
%
% [C,Atp,Btp]=etprod(Cidx,A,Aidx,B,Bidx);
% Wrapper function for tprod to map between Einstein summation
% convetion (ESC) and tprod's numeric calling convention e.g.
% 1) Matrix Matrix product: 
%   ESC: C_ij = A_ik B_kj <=> C = etprod('ij',A,'ik',B,'kj'); 
% 2) Vector outer product
%   ESC: C_ij = A_i B_j   <=> C = etprod('ij',A,'i',B,'j');   % A,B col vec
%                             C = etprod('ij',A,' i',B,' j'); % A,B row vec
% N.B. use spaces ' ' to indicate empty/ignored *singlenton* dimensions
% 3) Matrix vector product
%   ESC: C_i  = A_ik B_k  <=> C = etprod('i',A,'ik',B,'k');
% 4) Spatial Filtering
%   ESC: FX_fte = A_cte B_cf <=> C = etprod('fte',A,'cte',B,'cf') 
% OR:
%   C = etprod({'feat','time','epoch'},A,{'ch','time','epoch'},B,{'ch','feat'})
%
% Inputs:
%  Cidx -- the list of dimension labels for the output
%  A    -- [n-d] array of the A values
%  Aidx -- [ndims(A) x 1] (array, string, or cell array of strings) 
%          list of dimension labels for A array
%  B    -- [m-d] array of the B values
%  Bidx -- [ndims(B) x 1] (array, string or cell array of strings)
%          list of dimension labels for B array
% Outputs:
%  C    -- [p-d] array of output values.  Dimension labels are as in Cidx
%  Atp  -- [ndims(A) x 1] A's dimspec as used in the core tprod call
%  Btp  -- [ndims(B) x 1] B's dimspec as used in the core tprod call
%
% See Also: tprod, tprod_testcases
%
% Copyright 2006-     by Jason D.R. Farquhar (jdrf@zepler.org)
% Permission is granted for anyone to copy, use, or modify this
% software and accompanying documents for any uncommercial
% purposes, provided this copyright notice is retained, and note is
% made of any changes that have been made. This software and
% documents are distributed without any warranty, express or
% implied
if ( iscell(Aidx)~=iscell(Bidx) || iscell(Cidx)~=iscell(Aidx) )
 error('Aidx,Bidx and Cidx cannot be of different types, all cells or arrays');
end

Atp = zeros(size(Aidx));
Btp = zeros(size(Bidx));
% Map inner product dimensions, to unique *negative* index
for i=1:numel(Aidx)
   if ( iscell(Aidx) ) Bmatch = strcmp(Aidx{i}, Bidx);
   else                Bmatch =       (Aidx(i)==Bidx);
   end
   if ( any(Bmatch) ) Btp(Bmatch)=-i; Atp(i)=-i; end;
end
% Spaces/empty values in the input become 0's, i.e. ignored dimensions
if ( iscell(Aidx) ) 
   Btp(strcmp(' ',Bidx))=0;Btp(strcmp('',Bidx))=0;
   Atp(strcmp(' ',Aidx))=0;Atp(strcmp('',Aidx))=0;
else   
   Btp(' '==Bidx)=0;
   Atp(' '==Aidx)=0; 
end

% Map to output position numbers, to correct *positive* index
for i=1:numel(Cidx);
   if ( iscell(Aidx) ) 
      Atp(strcmp(Cidx{i}, Aidx))=i;
      Btp(strcmp(Cidx{i}, Bidx))=i;
   else
      Atp(Cidx(i)==Aidx)=i;
      Btp(Cidx(i)==Bidx)=i;
   end
end

% now do the tprod call.
global LOG; if ( isempty(LOG) ) LOG=0; end; % N.B. set LOG to valid fd to log
if ( LOG>0 ) 
   fprintf(LOG,'tprod(%s,[%s], %s,[%s])\n',mxPrint(A),sprintf('%d ',Atp),mxPrint(B),sprintf('%d ',Btp));
end
C=tprod(A,Atp,B,Btp,'n');
return;

function [str]=mxPrint(mx)
sz=size(mx);
if ( isa(mx,'double') ) str='d'; else  str='s'; end;
if ( isreal(mx)) str=[str 'r']; else str=[str 'c']; end;
str=[str ' [' sprintf('%dx',sz(1:end-1)) sprintf('%d',sz(end)) '] '];
return;

%----------------------------------------------------------------------------
function testCase();
A=randn(10,10); B=randn(10,10);
C2 = tprod(A,[1 -2],B,[-2 2]);
C  = etprod('ij',A,'ik',B,'kj'); mad(C2,C)
C  = etprod({'i' 'j'},A,{'i' 'k'},B,{'k' 'j'});
A=randn(100,100);B=randn(100,100,4);
C3 = tprod(A,[-1 -2],B,[-1 -2 3]);
C3 = tprod(B,[-1 -2 1],A,[-1 -2]);
C3 = tprod(A,[-1 -2],B,[-1 -2 1]);
C  = etprod('3',A,'12',B,'123');
C  = etprod([3],A,[1 2],B,[1 2 3]);
C  = etprod({'3'},A,{'1' '2'},B,{'1' '2' '3'})
