function out=icassoGet(sR,field,index)
%function out=icassoGet(sR,field,[index])
%
%PURPOSE 
%
%Auxiliary function for obtaining various information from the
%Icasso result data structure. Using icassoGet, you can return
%information related to original data, e.g.: data matrix,
%(de)whitening matrix, total number of estimates, number of
%estimates on each round. You can also return specified estimated
%independent components (sources), and rows of demixing matrices
%from. (However, it is easier to use function icassoResult to return
%the final estimation results from the complete Icasso procedure.) 
%
%EXAMPLES OF BASIC USAGE
%
%   M=icassoGet(sR,'m');        
% 
%returns total number of estimates.
%
%   nic=icassoGet(sR,'numOfIC')
%
%returns number of estimated IC components on each round in an Nx1
%vector (The number may differ in deflatory estimation mode due to
%convergence problems.)
%
%   r=icassoGet(sR,'round',[25 17 126]); 
% 
%workspace variable r contains now a 3x2 matrix where the first 
%column shows the number of estimation cycle. The second column 
%is the order of appearance of the estimate within the cycle: 
%e.g. 2  5 : estimate 25 was 5th estimate in cycle 2
%     1 17 :          17     1st                   1
%     7  6 :         126     6th                   7
%
%   W=icassoGet(sR,'demixingmatrix',[25 17 126]); 
%
%returns the demixing matrix rows for estimates 25, 17, and 126.
% 
%   s=icassoGet(sR,'source'); 
%
%returns _all_ M estimated independent components (sources), i.e.,
%estimates from all resampling cycles.
%
%INPUT
%
%[An argument in brackets is optional. If it isn't  given or it's
% an empty matrix/string, the function will use a default value.] 
%
% sR    (struct) Icasso result data structure (from icassoEst).
% field (string) (case insensitive) determines what is returned: 
%    'data'           the original data centered 
%    'N'              the number of estimation cycles
%    'M'              the total number of estimates available;
%                     equal to sum(icassoGet(sR,'numofic'));
%    'numofic'        number of ICs extracted on each cycle (Nx1
%                     vector) 
%    'odim'           original data dimension
%    'rdim'           reduced data dimension (after PCA)
%    'round'          output is an Mx2 matrix that identifies from
%                     which estimation  round each estimate
%                     originates:  
%                     out(i,1) is the estimation round for 
%                     estimate i, 
%                     out(i,2) is the ordinal number of the
%                     estimate within the round.
%    'source'         estimated source signals 
%    'dewhitemat'     dewhitening matrix for the original data
%    'whitemat'       dewhitening matrix for the original data
%    'demixingmatrix' demixing matrix rows ('W' works as well)
%
% [index] (vector of integers) specifies which estimates to return 
%          Applies only for 'demixingmatrix', and 'source'
%          Default: leaving this argument out corresponds to
%          selecting all estimates, i.e., giving [1:M] as index
%          vector. 
%
%OUTPUT
%
% out   (type and meaning vary according to input argument 'field')
%
%SEE ALSO
% icassoResult


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

% ver 1.2 100105

M=size(sR.index,1);

if nargin<3,
   index=1:M;
end

switch lower(field)
 case 'm'
   out=size(sR.index,1);
 case 'n'
  out=length(sR.A);
 case 'odim'
  out=size(sR.whiteningMatrix,2);
 case 'rdim'
  out=size(sR.whiteningMatrix,1);
 case 'numofic'
  for i=1:length(sR.A),
    out(i,1)=size(sR.A{i},2);
  end
 case {'w','demixingmatrix'}
  out=getDeMixing(sR,index);   
 case 'data'
  out=sR.signal;   
 case {'dewhitemat'}
  out=sR.dewhiteningMatrix;
 case {'whitemat'}
  out=sR.whiteningMatrix;
 case 'source'
  out=getSource(sR,index);
 case 'round'
  out=sR.index(index,:)
 otherwise  
   error('Unknown operation.');
end

%function A=getMixing(sR,index);
%
%function A=getMixing(sR,index);
%

% reindex
%index2=sR.index(index,:);
%N=size(index2,1); A=[];
%for i=1:N,
%   A(:,i)=sR.A{index2(i,1)}(:,index2(i,2));
%end
   
function W=getDeMixing(sR,index);
%
%function W=getDeMixing(sR,index);
%

% reindex
index2=sR.index(index,:);
N=size(index2,1); W=[];
for i=1:N,
   W(i,:)=sR.W{index2(i,1)}(index2(i,2),:);
end
   
function S=getSource(sR,index)
%function S=getSource(sR,index)
%

X=sR.signal;
W=getDeMixing(sR,index);
S=W*X;

%% Old stuff 
%function dWh=getDeWhitening(sR,index);
%
%function dWh=getDeWhitening(sR,index);
%
%
%index=sR.index(index,:);
%Nindex=size(index,1); W=[];
%for i=1:Nindex,
%   dWh(:,i)=sR.dewhiteningMatrix{1}(:,index(i,2));
%end


%function W=getWhitening(sR,index);
%
%function W=getWhitening(sR,index);
%
%
%index=sR.index(index,:);
%Nindex=size(index,1); W=[];
%for i=1:Nindex,
%   W(:,i)=sR.whiteningMatrix{index(i,1)}(:,index(i,2));
%end

