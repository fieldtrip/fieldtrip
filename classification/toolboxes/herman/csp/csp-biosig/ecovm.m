function [XCN,NN,mu,sd,COV,corrcoef,M,R2] = ecovm(X,Y);
% ECOVM produces an extended Covariance matrix, 
% function [ECM] = ecovm(X);
% ECM= [1 X]'*[1 X]; % l is a matching column of 1's
% ECM is additive, i.e. it can be applied to subsequent blocks and summed up afterwards
%
% function [ECM] = ecovm(X,Y);
% ECM= [1 X]'*[1 Y]; % l is a matching column of 1's
%
% if   [ECM_1,NN1] = ecovm(s_1);
% and  [ECM_n,NN2] = ecovm(s_n);
% with [ECM,NN]   = ecovm([s_1+...+s_n]);
% then ECM_1 + ... + ECM_n == ECM
% and  NN1 + ... + NNn == NN
%
% see also: DECOVM.M, R2.M

%	$Revision$
%	$Id$
%	Copyright (c) 1999-2002, 2004 by Alois Schloegl
%	a.schloegl@ieee.org	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

% V0.10	20.10.1998
% V0.11	10.11.1999    Made Octave compatible
% V0.12	26.04.2000    ecovm.m included
% V0.13	02.05.2000    correlation coefficient R2 included
% V0.20 27.03.2001    missing values (i.e. NaN) considered
% V0.21 12.04.2001    extended to Y 
% V0.22 08.06.2001    Y is used as class information 
% V0.23 06.01.2002    bug fixed	

% see also: NN, SUM, MEAN, STD, VAR, SUMSQ, MEANSQ, RMS, COV, 
%  CORRCOEF, R2: might be incorrect in case of missing values, 
%         especially in case of a correlation with missing value

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

if nargin==1,
        [r1,c1]=size(X);
        [r2,c2]=size(X);
        if (c1>r1),
                fprintf(2,'Warning ECOVM: Covariance is ill-defined, because of too few observations (rows).\n');
        end;
        tmp=[ones(r1,1),~isnan(X)];
        NN = tmp'*tmp;
        
        X(isnan(X)) = 0; % skip NaN's
        tmp = [ones(r1,1),X];
        XCN = tmp'*tmp;
elseif nargin==2
        [h,x]=histo(Y,1);
        for k=1:length(x),
                tmp = [ones(h(k),1),X(Y==x(k),:)];
                XCN{k} = tmp'*tmp;
        end;
elseif 0
        [r1,c1]=size(X);
        [r2,c2]=size(Y);
        if r1~=r2,
                fprintf(2,'Error ECOVM: X and Y must have the same number of observations (rows).\n');
                return;
        end;
        if (c1>r1) | (c2>r2),
                fprintf(2,'Warning ECOVM: Covariance is ill-defined, because of too few observations (rows).\n');
        end;
                
        NN = [ones(r1,1),~isnan(X)]'*[ones(r2,1),~isnan(Y)];
        X(isnan(X)) = 0; % skip NaN's
        Y(isnan(Y)) = 0; % skip NaN's
        XCN = [ones(r1,1),X]'*[ones(r2,1),Y];
else
        
end;

if nargout==2.
        return;
elseif nargout<2
        XCN=XCN./max(NN-[zeros(1,c2+1);zeros(c1,1),ones(c1,c2)],0);
else %if nargout>2
        XCN = XCN./max(NN-[zeros(1,c2+1);zeros(c1,1),ones(c1,c2)],0);
        %mu = XCN(2:c1+1,1);
        mu  = XCN(1,2:c2+1);          %nonsense bei (X,Y)
        COV = XCN(2:c1+1,2:c2+1) - XCN(2:c1+1,1)*XCN(1,2:c2+1); % OK
        if nargin<2
                sd  = sqrt(diag(COV))';       
                corrcoef = COV./(sd'*sd);           
        else
                corrcoef = COV./(std(X,1)'*std(Y,1));
        end;
        R2  = corrcoef .*corrcoef ;             
end;
return;

if r>1
        SS  = sum(X); 
        XCN = [[size(X,1),SS]; [SS',X'*X]];
else
        XCN = [[1,X]; [X',X'*X]];
end;

if nargout>1
        %[mu,sd,COV,xc,M]=decovm(XCN)
        N=size(XCN,2);
        M=XCN(1,1);
        mu=XCN(2:N,1)/XCN(1,1);
        COV=(XCN(2:N,2:N)/XCN(1,1)-XCN(2:N,1)*XCN(1,2:N)/XCN(1,1)^2);
        sd=sqrt(diag(COV));
        xc=COV./(sd*sd');
        R2=xc.*xc;
end;
