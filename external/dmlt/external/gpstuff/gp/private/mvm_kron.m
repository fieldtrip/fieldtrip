function x = mvm_kron(pj, n, Ka, Kb, H, B_m, jitterSigma2,b)
%MVM_KRON Fast matrix vector multiplication using Kronecker product
%         for Logistic-Gaussian Process density model in two
%         dimensional cases
%
%  See also
%    GPLA_ND_E, LGPDENS, DEMO_LGPDENS
%
% Copyright (c) 2012 Jaakko Riihimäki

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

pjsq=sqrt(pj);
Rb=sqrt(n)*(-pj*(pj'*(b./pjsq))+pjsq.*b);

n1=size(Ka,1);
n2=size(Kb,1);

q=Kb*reshape(Rb,n2,n1)*Ka;
q=q(:);

q=q+jitterSigma2*Rb;

if ~isempty(H)
  q=q+H'*(B_m*(H*Rb));
end
x=b+sqrt(n)*( pjsq.*q - (pj*(pj'*q))./pjsq );

