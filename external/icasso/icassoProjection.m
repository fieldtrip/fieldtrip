function sR=icassoProjection(sR, method, varargin)
%function sR=icassoProjection(sR,[method],['identifier1',val1,'identifier2',val2,...]))
%
%PURPOSE
%
%To project points on plane so that Euclidean distances between the
%projected points correspond to the similarity matrix between IC
%estimates in the Icasso result structure. 
%
%EXAMPLES OF BASIC USAGE
%
%   sR=icassoProjection(sR); 
%
%makes a CCA projection using default parameters. This is equivalent
%of giving command:  
%
%   sR=icassoProjection(sR,'cca', 's2d','sim2dis2','epochs',70,'alpha',0.7); 
%
%In the following example 'sim2dis' (i.e., D=1-S) is used to make
%the similarity-to-dissimilarity transformation and a longer,
%sammon's projection is used and a longer training sequence is engaged:
%
%   sR=icassoProjection(sR,'sammon','s2d','sim2dis','epochs',200):
%
%INPUT
%
%[An argument in brackets is optional. If it isn't  given or it's
% an empty matrix/string, the function will use a default value.] 
%
% sR       (struct) Icasso result data structure
% [method] (string) 'cca' (defalut) | 'mmds' | 'sammon' 
%
%Optional input arguments are given as argument identifier - value
%pairs: 'identifier1', value1, 'identifier2', value2,... 
%(case insensitive)
%
%  'epochs'     (scalar) see details 
%  'alpha',     (scalar) see details
%  'radius'     (scalar) see details
%  's2d'        (string) see details 
%
%OUTPUT
%
% sR (struct) updated Icasso result data structure, 
%
%The function updates _only_ the following fields:
%sR.projection.method, sR.projection.parameters, and
%sR.projection.coordinates. See icassoStruct. 
%
%DETAILS
%
%The function transforms the similarities S in field
%sR.cluster.similarities into dissimilarities using function
%sim2dis2 as default. Note that this may be different from
%the transformation that was used for making clustering. 
%Input argument pair 's2d',<string> gives the name of function
%that is used to make the transformation from similarities
%S=sR.cluster.similarity. There are two ready made functions
%sim2dis.m and sim2dis2.m that can be used: 
%   function      makes transformation      
%  'sim2dis2'     D=sqrt(1-S)   (default)
%  'sim2dis'      D=1-S 
%or you can specify your own function.
%
%The function can use three methods to do the projection on D 
%1. Curvilinear Component Analysis (CCA) (preferred)
%2. Principal Coordinates (linear Metric Multi-Dimensional Scaling,MMDS) 
%3. Sammon's projection (Sammon) 
%CCA and Sammon require some parameters. The default values are set
%according to experience and can be altered if necessary. MMDS is
%automatically used as an initial projection for both Sammon and CCA.   
%
%'epochs'  (scalar) the number of epochs for CCA or Sammon.
%           Ignored for MMDS.
%           default(s): CCA: 75, Sammon: 100
%            
%'radius'  (scalar) CCA initial radius. Ignored for Sammon and MMDS
%            default(s) CCA: max(10,M/20) where M is the number of
%            estimates
%             
%'alpha'   (scalar) learning rate factor for CCA or Sammon.
%            ignored for MMDS
%            default(s): CCA: 0.7, Sammon: 0.7
%
%SEE ALSO
% mmds
% cca (in SOM Toolbox)
% sammon (in SOM Toolbox)

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

% ver 1.21 030305 johan

% Set default projection method
if nargin<2|isempty(method),
  method='cca';
end

% Check the method
method=lower(method);
switch method
 case {'sammon','cca','mmds'}
  ;
 otherwise
  error('Unknown projection.');
end

% We project onto plane
outputDimension=2;

% Set default parameters for proj, methods
switch method 
 case 'sammon'
  default={'alpha',0.7,'epochs',100,'s2d','sqrtsim2dis'};
 case 'cca'
  default={'alpha',0.7,'epochs',75,...
	   'radius',max(icassoGet(sR,'M')/20,10),'s2d','sqrtsim2dis'};
 case 'mmds'
  default={'s2d','sqrtsim2dis'};
end

%% Check optional arguments and add defaults
projectionparameters=processvarargin(varargin,default);
num_of_args=length(projectionparameters);

%% check arguments
for i=1:2:num_of_args;
  switch lower(projectionparameters{i})
   case 's2d'
    sim2dis=projectionparameters{i+1};
   case 'epochs'
    epochs=projectionparameters{i+1};
   case 'alpha'
    alpha=projectionparameters{i+1};
   case 'radius'
    CCAradius=projectionparameters{i+1};
   otherwise
    error(['Indentifier ' projectionparameters{i} ' not recognized.']);
  end
end

% Make similarity-to-dissimilarity transformation

D=feval(sim2dis,sR.cluster.similarity);

disp([char(13) 'Projection, using ' upper(method) char(13)]);

switch method 
 case 'mmds'
  P=mmds(D); 
  P=P(:,1:outputDimension);
 otherwise
  % Start from MMDS
  initialProjection=mmds(D); initialProjection=initialProjection(:,1:2);
  % 
  dummy=rand(size(D,1),outputDimension);
  % rand. init projection: set 
  % initialProjection=dummy;
  switch method
   case 'sammon' % Use SOM Toolbox Sammon 
    P=sammon(dummy,initialProjection,epochs,'steps',alpha,D);
   case 'cca'    % Use SOM Toolbox CCA 
    P=cca(dummy,initialProjection,epochs,D,alpha,CCAradius);
  end
end

sR.projection.method=method;
sR.projection.parameters=projectionparameters;
sR.projection.coordinates=P;


