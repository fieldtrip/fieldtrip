function [mu,sd,COV,xc,M,R2]=decovm(XCN,NN)
% decompose extended covariance matrix into mean (mu), 
% standard deviation, the (pure) Covariance (COV), 
% correlation (xc) matrix and the correlation coefficients R2.
% NaN's are condsidered as missing values. 
% [mu,sd,COV,xc,N,R2]=decovm(ECM[,NN])
%
% ECM 	is the extended covariance matrix
% NN	is the number of elements, each estimate (in ECM) is based on 
%
% see also: MDBC, ECOVM, R2

%	Version 0.22
%	13.10.2002
%	Copyright (c) 1999-2002 by  Alois Schloegl
%	a.schloegl@ieee.org	

% V0.10	20.10.1998
% V0.11	10.11.1999    Made Octave compatible
% V0.12	26.04.2000    ecovm.m included
% V0.13	02.05.2000    correlation coefficient R2 included
% V0.20 27.03.2001    missing values (i.e. NaN) considered

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the  License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

[r,c]=size(XCN);
if r~=c,
        fprintf(2,'Warning DECOVM: input argument is not a square matrix\n');
        XCN = ecovm(XCN);
        c = c + 1
else
	M = XCN(1,1);
	if nargin<2,
                XCN = XCN/(XCN(1,1));
        else %if nargin==2
                XCN = XCN./(NN);
        end;

        if any(isnan(XCN(:))),
	        warning('DECOVM: Extended Covariance Matrix should not contain NaN''s');
	end;
        if det(XCN)<0,
                warning('DECOVM: Extended Covariance Matrix must be non-negative definite');
        end;
end;

mu  = XCN(1,2:c);
COV = XCN(2:c,2:c) - mu'*mu;
sd  = sqrt(diag(COV))';        
xc  = COV./(sd'*sd);
R2  = xc.*xc;

return;
        
M=XCN(1,1);
mu=XCN(2:N,1)/XCN(1,1);
COV=(XCN(2:N,2:N)/XCN(1,1)-XCN(2:N,1)*XCN(1,2:N)/XCN(1,1)^2);
sd=sqrt(diag(COV));
xc=COV./(sd*sd');
R2=xc.*xc;

% function [ECM] = ecovm(signal);
% Generates extended Covariance matrix, 
% ECM= [l signal]'*[l signal]; % l is a matching column of 1's
% ECM is additive, i.e. it can be applied to subsequent blocks and summed up afterwards
% [ECM1] = ecovm(s1);
% [ECM2] = ecovm(s1);
% [ECM]  = ecovm([s1;s2]);
% ECM1+ECM2==ECM;
%
% SS=sum(signal); ECM=[[size(signal,1),SS];[SS',signal'*signal]];
