function [unmixing] = csp(dat1, dat2)

% CSP Common spatial pattern decomposition
%
% Use as
%   [unmixing] = csp(dat1, dat2)
%
% This implements Ramoser, H., Gerking, M., and Pfurtscheller, G. "Optimal
% spatial filtering of single trial EEG during imagined hand movement."
% IEEE Trans. Rehab. Eng 8 (2000), 446, 441.

% The initial version was coded by James Ethridge and William Weaver.
% See http://www.mathworks.com/matlabcentral/fileexchange/22915-common-spatial-patterns
% Some cleanups by Robert Oostenveld, 2012

R1 = dat1*dat1';
R1 = R1/trace(R1);
R2 = dat2*dat2';
R2 = R2/trace(R2);

% Ramoser equation (2)
Rsum = R1+R2;

% Find Eigenvalues and Eigenvectors of RC
% Sort eigenvalues in descending order
[EVecsum,EValsum] = eig(Rsum);
[EValsum,ind] = sort(diag(EValsum),'descend');
EVecsum = EVecsum(:,ind);

% Find Whitening Transformation Matrix - Ramoser Equation (3)
W = sqrt(pinv(diag(EValsum))) * EVecsum';

% Whiten Data Using Whiting Transform - Ramoser Equation (4)
S1 = W * R1 * W';
S2 = W * R2 * W';

% Ramoser equation (5)
% [U1,Psi1] = eig(S1);
% [U2,Psi2] = eig(S2);

%generalized eigenvectors/values
[B,D] = eig(S1,S2);

% Simultanous diagonalization
% Should be equivalent to [B,D]=eig(S1);

% verify algorithim
%disp('test1:Psi1+Psi2=I')
%Psi1+Psi2

% sort ascending by default
%[Psi1,ind] = sort(diag(Psi1)); U1 = U1(:,ind);
%[Psi2,ind] = sort(diag(Psi2)); U2 = U2(:,ind);
[D,ind]=sort(diag(D));
B=B(:,ind);

% Resulting Projection Matrix-these are the spatial filter coefficients
unmixing = B'*W;
