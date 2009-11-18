function [C,Atp,Btp]=etprod(Cidx,A,Aidx,B,Bidx)
% Wrapper for tprod to make calls more similar to Einstein Summation
% Convention
%
% [C,Atp,Btp]=etprod(Cidx,A,Aidx,B,Bidx);
% Wrapper function for tprod to map between Einstein summation
% convetion (ESC) and tprods internal calling convention e.g.
% 1) Matrix Matrix product: 
%   ESC: C_ij = A_ik B_kj <=> C = etprod('ij',A,'ik',B,'kj'); 
% 2) Vector outer product
%   ESC: C_ij = A_i B_j   <=> C = etprod('ij',A,'i',B,'j');   % A,B col vec
%                             C = etprod('ij',A,' i',B,' j'); % A,B row vec
% N.B. use spaces ' ' to indicate empty/ignored dimensions
% 3) Matrix vector product
%   ESC: C_i  = A_ik B_k  <=> C = etprod('i',A,'ik',B,'k');
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

Atp = zeros(size(Aidx));
Btp = zeros(size(Bidx));

% Map inner product dimensions, to unique *negative* index
for i=1:numel(Aidx)
  Bmatch=(Bidx==Aidx(i));
  if ( any(Bmatch) ) Btp(Bmatch)=-i; Atp(i)=-i; end;
end
% Spaces in the input become 0's, i.e. ignored dimensions
Btp(Bidx==' ')=0;
Atp(Aidx==' ')=0; 

% Map to output position numbers, to correct *positive* index
for i=1:numel(Cidx);
   Atp(Aidx==Cidx(i))=i;
   Btp(Bidx==Cidx(i))=i;
end

% now do the tprod call.
%fprintf('C=tprod(A,[%s],B,[%s])\n',sprintf('%d ',Atp),sprintf('%d ',Btp));
C=tprod(A,Atp,B,Btp,'n');
return;

%----------------------------------------------------------------------------
function testCase();
A=randn(10,10); B=randn(10,10);
C = etprod('ij',A,'ik',B,'kj');