function ir=icassoRindex(sR,L)
%function ir=icassoRindex(sR,[L])
%
%PURPOSE
%
%To return and/or plot a relative clustering validity index.
%R-index is originally defined in Levine & Dormay (2001),
%"Resampling method for unsupervised estimation of cluster
%validity. Neural Computation 13(11)) and redescribed in
%publication Himberg et al. (2004), "Validating the independent
%components of neuroimaging time-series via clustering and
%visualization". NeuroImage, 22:3(1214-1222).    
%
%INPUT 
%
%[An argument in brackets is optional. If it isn't  given or it's
% an empty matrix/string, the function will use a default value.] 
%
% sR  (struct) Icasso result data struct; estimate-clusters have to
%      be readily computed (icassoCluster)
% [L] (integer) number of estimate-clusters up to which the R-index
%      is plotted, default: at most the total number of estimates, at
%      least 1.5 times the reduced data dimension
%
%OUTPUT
%
% ir (vector) the computed NaN means that the index has not (or
%      cannot) be computed for that number of clusters
%
%DETAILS
%
%Local minimums of the index suggest a "good" number of clusters;
%the index is not well-established and acts only as a guideline!  
%
%The plot shows other information in addition to R-index: original
%(and reduced) data dimension and the (maximum) number of
%independent components extracted in a single resampling round. 
%
%R-index is not defined for a cluster consisting of just one item;
%also the function that computes R-index (rindex) limit the
%number of clusters for which the R-index is computed in order to
%reduce the computational load. 
%
%SEE ALSO
% rindex
% icassoCluster

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

% ver 1.21 johan 070305

% Get number of original signals and after PCA
%
ODim=icassoGet(sR,'odim');
RDim=icassoGet(sR,'rdim');
nic=icassoGet(sR,'numOfIC');

if nargin<2|isempty(L)
  L=RDim;
end

N=max(L,min(icassoGet(sR,'M'),ceil(RDim*1.5)));

%% Get R index
R=sR.cluster.index.R; R=R(:);

%% Get the point up to which R is computed
n=max(find(isfinite(R)));

if isempty(n),
  warning('No values for the R-index');
  n=1;
end

%%% normalize max to 1
R=R(1:n)./max(R);

%%%% Plot R-index

cla reset; hold on;

if RDim<N,
  set(patch([RDim RDim N+1 N+1],...
	    [0 1.1 1.1 0],...
	    [.9 .9 .9]),...
      'edgecolor','none');
end

set(plot(R(1:n),'o-'),'linewidth',3,'color',[.7 0 0]);

xlabel('Number of clusters');
title(['R-index computed for 2-' num2str(n) ' clusters']);
axis([0.5 N+.5 0 1.05]);

ax=axis; 

set(gca,'yticklabel',[]); 
ylabel('(Absolute scale has no relevance)');

set(plot([L L],[0 1],'r-'),'linewidth',3);
set(text(L,1,[num2str(L) ' clusters selected']),...
    'color','r','vert','bot','fontsize',12);

%% Plot number of ICA components 
if max(nic)==min(nic),
  set(plot([min(nic) min(nic)],[0 1],'b-'),'linewidth',2);
  set(text(min(nic),0.05,{['Number of ICs = ' num2str(min(nic))]}),...
      'color','b','rot',90,'vert','bot');
else
  
  set(plot([min(nic) min(nic)],[0 1],'b-'),'linewidth',2);
  set(text(min(nic),0.05,['Min. number of ICs = ' num2str(min(nic))]),...
		    'color','b','rot',90,'vert','top');

  set(plot([max(nic) max(nic)],[0 1],'b-'),'linewidth',2);
  set(text(max(nic),0.05,['Max. number of ICs = ' num2str(max(nic))]),...
      'color','b','rot',90,'vert','top');
end

%%% Plot data dimension
if ODim==RDim,
  set(plot([ODim ODim],[0 1],'k-'),'linewidth',2);
  set(text(ODim,0.5,['Data dimension = ' num2str(ODim)]),'color','k','rot',90,'vert','bot');
else
  set(plot([ODim ODim],[0 1],'k-'),'linewidth',2);
  set(text(ODim,0.5,['Orig. data dimension = ' num2str(ODim)]),'color','k','rot',90, ...
		    'vert','bot');
  
  set(plot([RDim RDim],[0 1],'k-'),'linewidth',2);
  set(text(RDim,0.5,['Data dim. after PCA = ' num2str(RDim)]),'color','k','rot',90,'vert','bot');
end

x=get(gca,'xtick'); x(round(x)~=x)=[]; set(gca,'xtick',x);
set(gcf,'name','Icasso: Clustering Quality');

if nargout>0,
  ir=sR.cluster.index.R;
end




