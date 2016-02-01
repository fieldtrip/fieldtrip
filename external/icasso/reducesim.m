function S=reducesim(S,minLimit,partition,scatter,scatterLimit)
%function S=reducesim(S,minLimit)
% or
%function S=reducesim(S,minLimit,partition,scatter,scatterLimit)
% 
%PURPOSE
%
%To reduce the number of lines that are drawn when the similarity
%matrix is visualized. That is, to set similarities below certain
%threshold to zero, and optionally, also within-cluster
%similarities above certain threshold to zero.
%
%INPUT
%
% Two input arguments:
%
% S        (matrix) NxN similarity matrix
% minLimit (scalar) all values below this threshold in S are set to
%            zero
%
% Five input arguments:
%
% S         (matrix) NxN similarity matrix
% minLimit  (scalar) all values below this threshold in S are set to
%             zero
% partition (vector) Nx1 partition vector into K clusters (see
%            explanation for 'partition vector' in function hcluster)
% scatter   (vector) Kx1 vector that contains within-scatter for
%             each cluster i=1,2,...,K implied by vector partition 
% scatterLimit (scalar) threshold for within-cluster scatter.
%
%
%OUTPUT
%
% S (matrix) NxN similarity matrix with some entries cut to zero.
%
%SEE ALSO
% clusterstat
% icassoShow

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

if nargin==2|nargin==5,
  ;
else
  error('You must specify 2 or 5 input arguments!');
end

%%% Maximum number of lines
MAX_NLINES=5000;

%%% Set within-cluster similarities above certain threshold to zero
if nargin==5,
  Ncluster=max(partition);
  for i=1:Ncluster,
    index=partition==i;
    if scatter(i)>scatterLimit;
      S(index,index)=0;
    end
  end
end

%%% We assume symmetric matrix and don't care about
%self-similarities upper diagonal is enough

S=tril(S); S(eye(size(S))==1)=0;

%% Number of lines to draw
Nlines=sum(sum(S>minLimit));

%%% If too many lines, try to set better minLimit
if Nlines>MAX_NLINES,
  warning('Creates overwhelming number of lines');
  warning('Tries to change the limit...');
  minLimit=reduce(S,MAX_NLINES);
  warning(['New limit =' num2str(minLimit)]);
end

%%% Set values below minLimit to zero
S(S<minLimit)=0;

function climit=reduce(c,N)

c(c==0)=[];
c=c(:);
if N==0,
  warning('N=0 not allowed; setting N=1');
end

if N>length(c),
  N=length(c);
end

climit=sort(-c); climit=-climit;
climit=climit(N);





