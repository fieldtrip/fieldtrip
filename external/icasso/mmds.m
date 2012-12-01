function [X,s,U,V]=mmds(D)
%
%function [X,s,U,V]=mmds(D)
%
%PURPOSE
%
%To compute principal coordinates (linear Metric Multi-Dimensional Scaling)
%
%INPUTS
%
% D 	(matrix) NxN matrix of dissimilarities
%
%OUTPUTS
%
% X   (matrix) NxN matrix of coordinates
% s   (vector) Nx1 vector of singular values corresponding to the rows of matrix X
% U,V (matrix) U,V from MATLAB function SVD 

%COPYRIGHT NOTICE
%This function is a part of Icasso software library
%Copyright (C) 2003-2005 Johan Himberg
%
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% ver 1.2 100105 johan

n=size(D,1);

% Squared dissimilarities
D2=D.^2;

% Centering matrix
Z=eye(n)-ones(n)*(1./n);

% Double centered inner product 
B=-0.5*Z*D2*Z;

% SVD
[U,S,V]=svd(B);

% Singular values
s=diag(S);

% Coordinates
X=U*diag(sqrt(s));

