function [CC,NN] = covm(X,Y,Mode);
% COVM generates covariance matrix
% X and Y can contain missing values encoded with NaN.
% NaN's are skipped, NaN do not result in a NaN output. 
% The output gives NaN only if there are insufficient input data
%
% COVM(X,Mode);
%      calculates the (auto-)correlation matrix of X
% COVM(X,Y,Mode);
%      calculates the crosscorrelation between X and Y
%
% Mode = 'M' minimum or standard mode [default]
% 	C = X'*X; or X'*Y correlation matrix
%
% Mode = 'E' extended mode
% 	C = [1 X]'*[1 X]; % l is a matching column of 1's
% 	C is additive, i.e. it can be applied to subsequent blocks and summed up afterwards
% 	the mean (or sum) is stored on the 1st row and column of C
%
% Mode = 'D' or 'D0' detrended mode
%	the mean of X (and Y) is removed. If combined with extended mode (Mode='DE'), 
% 	the mean (or sum) is stored in the 1st row and column of C. 
% 	The default scaling is factor (N-1). 
% Mode = 'D1' is the same as 'D' but uses N for scaling. 
%
% C = covm(...); 
% 	C is the scaled by N in Mode M and by (N-1) in mode D.
% [C,N] = covm(...);
%	C is not scaled, provides the scaling factor N  
%	C./N gives the scaled version. 

%	$Revision: 1.18 $
%	$Id: covm.m,v 1.18 2003/09/15 14:07:49 schloegl Exp $
%	CopyLeft (C) 2000-2002 by Alois Schloegl <a.schloegl@ieee.org>	

%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


if nargin<3,
        if nargin==2,
		if isnumeric(Y),
                        Mode='M';
                else
		        Mode=Y;
                        Y=[];
                end;
        elseif nargin==1,
                Mode = 'M';
                Y = [];
        elseif nargin==0,
                error('Missing argument(s)');
        end;
end;        

Mode = upper(Mode);

[r1,c1]=size(X);
if ~isempty(Y)
        [r2,c2]=size(Y);
        if r1~=r2,
                error('X and Y must have the same number of observations (rows).');
                return;
        end;
else
        [r2,c2]=size(X);
end;

if (c1>r1) | (c2>r2),
        warning('Covariance is ill-defined, because of too less observations (rows)');
end;

if ~isempty(Y),
        if (~any(Mode=='D') & ~any(Mode=='E')), % if Mode == M
        	NN = real(~isnan(X)')*real(~isnan(Y));
	        X(isnan(X)) = 0; % skip NaN's
	        Y(isnan(Y)) = 0; % skip NaN's
        	CC = X'*Y;
        else  % if any(Mode=='D') | any(Mode=='E'), 
	        [S1,N1] = sumskipnan(X,1);
                [S2,N2] = sumskipnan(Y,1);
                
                NN = real(~isnan(X)')*real(~isnan(Y));
        
	        if any(Mode=='D'), % detrending mode
        		X  = X - ones(r1,1)*(S1./N1);
                        Y  = Y - ones(r1,1)*(S2./N2);
                        if any(Mode=='1'),  %  'D1'
                                NN = NN;
                        else   %  'D0'       
                                NN = max(NN-1,0);
                        end;
                end;
                
                X(isnan(X)) = 0; % skip NaN's
        	Y(isnan(Y)) = 0; % skip NaN's
                CC = X'*Y;
                
                if any(Mode=='E'), % extended mode
                        NN = [r1, N2; N1', NN];
                        CC = [r1, S2; S1', CC];
                end;
	end;
        
else        
        if (~any(Mode=='D') & ~any(Mode=='E')), % if Mode == M
        	tmp = real(~isnan(X));
                NN  = tmp'*tmp; 
                X(isnan(X)) = 0; % skip NaN's
	        CC = X'*X;
        else  % if any(Mode=='D') | any(Mode=='E'), 
	        [S,N] = sumskipnan(X,1);
        	tmp = real(~isnan(X));
                NN  = tmp'*tmp; 
                if any(Mode=='D'), % detrending mode
	                X  = X - ones(r1,1)*(S./N);
                        if any(Mode=='1'),  %  'D1'
                                NN = NN;
                        else  %  'D0'      
                                NN = max(NN-1,0);
                        end;
                end;
                
                X(isnan(X)) = 0; % skip NaN's
                CC = X'*X;
                
                if any(Mode=='E'), % extended mode
                        NN = [r1, N; N', NN];
                        CC = [r1, S; S', CC];
                end;
	end
end;

if nargout<2
        CC = CC./NN; % unbiased
end;
return; 
