function [err, Qd, dfx, dmat2] = QRDecom(dmat, cmat)
%
%   [err,] = QRDecom.m ()
%
%Purpose:
%
%
%
%Input Parameters:
%
%
%
%Output Parameters:
%   err : 0 No Problem
%       : 1  Problems
%
%
%
%Key Terms:
%
%More Info :
%
%
%
%
%     Author : Gang Chen
%     Date : Tue Mar 23 16:09:23 EST 2004
%     SSCC/NIMH/ National Institutes of Health, Bethesda MD 20892


%Define the function name for easy referencing
FuncName = 'QRDecom.m';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;

% Tolerance for computing rank from diag(R) after QR decomposition
%[nrows,ncols] = size(dmat);

% Find the null space of the constraints matrix
[Qc,Rc,Ec] = qr(cmat');
pc = Rrank(Rc);
Qc0 = Qc(:,pc+1:end);

% Do qr decomposition on design matrix projected to null space
Dproj = dmat*Qc0;
[Qd,Rd,Ed] = qr(Dproj,0);
dfx = Rrank(Rd);
Qd = Qd(:,1:dfx);
%Rd = Rd(1:dfx,1:dfx);

% Return reduced design matrix if requested
if nargout>3
   dmat2 = Qd' * dmat;
end

err = 0;
return;
