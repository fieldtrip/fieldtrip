function [cluster, num] = findcluster(onoff, spatdimneighbstructmat, varargin)

% FINDCLUSTER returns all connected clusters in a 3 dimensional matrix
% with a connectivity of 6.
%
% Use as
%   [cluster, num] = findcluster(onoff, spatdimneighbstructmat, minnbchan)
% or ar
%   [cluster, num] = findcluster(onoff, spatdimneighbstructmat, spatdimneighbselmat, minnbchan)
% where 
%   onoff                   is a 3D boolean matrix with size N1xN2xN3
%   spatdimneighbstructmat  defines the neighbouring channels/combinations, see below
%   minnbchan               the minimum number of neighbouring channels/combinations 
%   spatdimneighbselmat     is a special neighbourhood matrix that is used for selecting
%                           channels/combinations on the basis of the minnbchan criterium
%
% The neighbourhood structure for the first dimension is specified using 
% spatdimneighbstructmat, which is a 2D (N1xN1) matrix. Each row and each column corresponds
% to a channel (combination) along the first dimension and along that row/column, elements
% with "1" define the neighbouring channel(s) (combinations). The first dimension of
% onoff should correspond to the channel(s) (combinations).
% The lower triangle of spatdimneighbstructmat, including the diagonal, is
% assumed to be zero. 
%
% See also BWSELECT, BWLABELN (image processing toolbox) 
% and SPM_CLUSTERS (spm2 toolbox).

% Copyright (C) 2004, Robert Oostenveld
%
% $Log: findcluster.m,v $
% Revision 1.4  2006/06/12 08:22:03  erimar
% Added changes to allow processing with cfg.minnbchan~=0.
%
% Revision 1.3  2006/03/21 15:02:42  roboos
% restructured the help, no functional change
%
% Revision 1.2  2005/05/17 17:50:49  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.1  2004/11/04 15:30:26  erimar
% Erimar's version of findcluster. Previously, findcluster was part of Roboos' misc.
%
% Revision 1.7 10/27/04  1:13 PM  erimar
% allow for the minnbchan option, allow for different matrices for
% selecting and clustering, some variable name changes
%
% Revision 1.6  2004/03/15 17:02:34  roberto
% fixed stupid bug at end of code due to renaming of variable
%
% Revision 1.5  2004/03/10 16:50:15  roberto
% added speed improvements of Eric (which accidentally were NOT in the previous revision)
% updated help, some cosmetical changes
%
% Revision 1.3  2004/02/25 10:34:10  roberto
% implemented smarter handling of real 3d data, improved help
%
% Revision 1.2  2004/02/23 15:41:03  roberto
% added global fb flag to determine amount of online feedback
% fixed bug in checking size of input matrices
%
% Revision 1.1  2004/02/17 14:30:37  roberto
% new efficient implementation to find clusters according to Eric Maris
%

spatdimlength = size(onoff, 1);
nfreq = size(onoff, 2);
ntime = size(onoff, 3);

if length(size(spatdimneighbstructmat))~=2 || ~all(size(spatdimneighbstructmat)==spatdimlength)
  error('invalid dimension of spatdimneighbstructmat');
end

minnbchan=0;
if length(varargin)==1
    minnbchan=varargin{1};
end;
if length(varargin)==2
    spatdimneighbselmat=varargin{1};
    minnbchan=varargin{2};
end;

if minnbchan>0
	% For every (time,frequency)-element, it is calculated how many significant
	% neighbours this channel has. If a significant channel has less than minnbchan
	% significant neighbours, then this channel is removed from onoff.
    
    if length(varargin)==1
        selectmat = single(spatdimneighbstructmat | spatdimneighbstructmat');
    end;
    if length(varargin)==2
        selectmat = single(spatdimneighbselmat | spatdimneighbselmat');
    end;
    nremoved=1;
	while nremoved>0
		nsigneighb=reshape(selectmat*reshape(single(onoff),[spatdimlength (nfreq*ntime)]),[spatdimlength nfreq ntime]);
		remove=(onoff.*nsigneighb)<minnbchan;
		nremoved=length(find(remove.*onoff));
		onoff(remove)=0;
	end;
end;

% for each channel (combination), find the connected time-frequency clusters
labelmat = zeros(size(onoff));
total = 0;
for spatdimlev=1:spatdimlength
  [labelmat(spatdimlev, :, :), num] = bwlabeln(reshape(onoff(spatdimlev, :, :), nfreq, ntime), 4);
  labelmat(spatdimlev, :, :) = labelmat(spatdimlev, :, :) + (labelmat(spatdimlev, :, :)~=0)*total;
  total = total + num;
end

% combine the time and frequency dimension for simplicity
labelmat = reshape(labelmat, spatdimlength, nfreq*ntime);

% combine clusters that are connected in neighbouring channel(s)
% (combinations).
replaceby=1:total;
for spatdimlev=1:spatdimlength
  neighbours=find(spatdimneighbstructmat(spatdimlev,:));
  for nbindx=neighbours
    indx = find((labelmat(spatdimlev,:)~=0) & (labelmat(nbindx,:)~=0));
    for i=1:length(indx)
      a = labelmat(spatdimlev, indx(i));
      b = labelmat(nbindx, indx(i));
      if replaceby(a)==replaceby(b)
        % do nothing
        continue;
      elseif replaceby(a)<replaceby(b)
        % replace all entries with content replaceby(b) by replaceby(a).
        replaceby(find(replaceby==replaceby(b))) = replaceby(a); 
      elseif replaceby(b)<replaceby(a)
        % replace all entries with content replaceby(a) by replaceby(b).
        replaceby(find(replaceby==replaceby(a))) = replaceby(b); 
      end
    end
  end
end

% renumber all the clusters
num = 0;
cluster = zeros(size(labelmat));
for uniquelabel=unique(replaceby(:))'
  num = num+1;
  cluster(ismember(labelmat(:),find(replaceby==uniquelabel))) = num;
end

% reshape the output to the original format of the data
cluster = reshape(cluster, spatdimlength, nfreq, ntime);

