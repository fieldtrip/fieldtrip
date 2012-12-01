function i=centrotype(S)
%function i=centrotype(S)
%
%PURPOSE
%
%To compute the centrotype of objects whose relations are
%defined by similarity matrix S. Centrotype means the object that
%is most similar to all the others, i.e., it should be close
%to the centroid of the data set.
%
%INPUT
% 
% S (matrix) NxN similarity matrix (between objects 1,2,...,N)
%
%OUTPUT
%
% i (scalar) object i is the centrotype of the objects
%
%SEE ALSO
% icassoIdx2Centrotype

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

% ver 1.2 johan 100105

[tmp,i]=max(sum(S));
