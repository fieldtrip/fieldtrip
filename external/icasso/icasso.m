function [Iq,A,W,S,sR]=icasso(X,M,varargin)
%function [Iq,A,W,S,sR]=icasso(X,M,['argid1',value1,'argid2',value2,...])
%
%PURPOSE
%
%Demonstrates the basic Icasso procedure: 
%1. Runs FastICA with given parameters M times on data X. in mode
%'both'. See function icassoEst.
%2. Clusters the estimates and computes other statistics: See
%function icassoExp.
%3. Returns (and visualizes) the best estimates. See functions
%icassoResult and icassoShow. 
%
%BASIC USAGE EXAMPLES
%
%Load MEG example dataset. 30 resampling cycles. Reduce
%data dimension to 20: use 'tanh' as FastICA
%nonlinearity. Visualize the results. 
%
%  load megdata; 
%  [iq, A, W, S, sR]=icasso(megdata,30,'g','tanh','lastEig',20);
%
%If you wish to skip the visualization and just return the output
%variables, use 
%
%  [iq, A, W, S, sR]=icasso(megdata,30,'g','tanh','lastEig',20,'vis','off');
%
%INPUT
%
% X    (dxN matrix) data where d=dimension, N=number of vectors 
% M    (scalar) number of randomizations (estimation cycles) 
%
%Optional input arguments are given as argument identifier - value
%pairs: 'identifier1', value1, 'identifier2', value2,... 
%(case insensitive)  
%
% FastICA parameters apply here (see function fastica) 
%
%Default: 
% 'approach', 'symm', 'g', 'pow3', 'maxNumIterations', 100, 'vis', 'basic'
%
%In addition to the FastICA paramters, there is an additional
%argument identifier 'vis': 
%
%'vis' (string) 'basic' (default) | 'off' 
% String 'basic': Visualizes the results. See details in help text
%of function icassoShow, and further in icassoGraph,
%icassoDendrogram, icassoRindex, icassoStability.  
% String 'off': Only returns Icasso results and the estimates. Does
%not visualize the estimate space at all and _does not_ compute the
%projection for visualization (icassoProjection is not run; this
%saves some time). 
%
%OUTPUT
%
% Iq    (vector) quality index of each estimate
% A     (matrix) estimated columns of the mixing matrix (A) = pinv(W) 
% W     (matrix) estimated rows of the demixing matrix (W) 
% S     (matrix) estimated independent components 
% sR    (struct) Icasso result data structure 
%
%Output variable sR, the Icasso result data structure, contains all
%relevant information from the resampling and clustering process
%for further analyses where you can utilize other functions of
%Icasso toolbox. Function 'icasso' readily returns estimates of A
%(mixing matrix), W (demixing matrix), and  S (independent
%components), and the stability index of these (Iq ;see function
%icassoStability).       
%
%(Rows of) W correspond to the centroids (actually centrotype; see help
%in function 'centrotype' for details) of the
%estimate-clusters. These estimates represents better the "true"
%estimate than an arbitrary estimate from a single run. However,
%the resulting W does not present a strictly orthogonal base in the
%whitened space. S is computed by W*X where X is the original data
%centered. A is computed as a pseudoinverse of W (A=pinv(W) using
%MATLAB notation). Iq is the heuristic quality (stability) of the
%estimates. Iq(2) corresponds to S(2,:) and W(2,:). Ideally, Iq
%should be 1 for each estimate. 

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

% ver 1.21 johan 040305

%%%%%% Initiate & check some varaibles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options=varargin;
fasticaoptions={};
L=[];
vis='basic';

j=1;
for i=1:2:length(varargin)
  switch varargin{i}
   case 'vis'
    vis=varargin{i+1};
    switch lower(vis),
     case {'off','basic'}
      ; % ok
     otherwise
      error('Option ''vis'' must be ''basic'' or ''off''.');
    end
   otherwise
    fasticaoptions(j:j+1)=varargin(i:i+1);
    j=j+2;
  end
end

%%% 1. Estimate ICA: icassoEst

if isempty(fasticaoptions)
  sR=icassoEst('both',X,M)
else
  sR=icassoEst('both',X,M,fasticaoptions{:});
end

%%% 2. Compute & visualize/return results

L=icassoGet(sR,'rdim');

switch vis
 case 'basic'
  sR=icassoExp(sR);
  [Iq,A,W,S]=icassoShow(sR,'L',L);
 case 'off'
  sR=icassoCluster(sR);
  [Iq,A,W,S]=icassoResult(sR,L);
 %case 'interactive' % not yet ready
 % sR=icassoExp(sR);
 % icassoViz(sR);
 otherwise
  error('Option ''vis'' must be ''basic'' or ''off''.');
end

