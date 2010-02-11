function [centers,sigmas] = subclustering(X,radii,xBounds,options)
%SUBCLUST Locates data cluster centers using subtractive clustering.
%
%   SUBCLUST operates by finding the optimal data point to define a cluster
%   center, based on the density of surrounding data points. All data points
%   within the distance RADII of this point are then removed, in order to
%   determine the next data cluster and its center. This process is repeated
%   until all of the data is within the distance RADII of a cluster center.
%
%   [C] = SUBCLUST(X,RADII) clusters the data points in the M-by-N matrix X,
%   where M is the number of data points and N is the number of coordinates
%   (data dimensions) that specify each point. RADII has a value between 0 and 1
%   and specifies the size of the cluster in each of the data dimensions,
%   assuming the data fall within a unit hyperbox (range [0 1]). Specifying a
%   smaller cluster radius will usually yield more (smaller) clusters in the
%   data. When RADII is a scalar it is applied to all data dimensions, when it
%   is a vector, it has one entry for each data dimension. The cluster centers
%   are returned as rows of the matrix C. C has size J-by-N, if J clusters are
%   required to cluster the data.
%
%   [C] = SUBCLUST(...,XBOUNDS) also specifies a matrix XBOUNDS, of size
%   2-by-N, used to normalize the data X into a unit hyperbox (range [0 1]).
%   Each column of XBOUNDS provides the minimum and maximum values for the
%   corresponding data dimension. If XBOUNDS is an empty matrix or not provided,
%   the minimum and maximum data values found in X, are used as defaults.
%
%   [C] = SUBCLUST(...,OPTIONS) specifies a vector for changing the default
%   algorithm parameters:
%      OPTIONS(1):  The squash factor, is used to multiply the RADII values to
%                   determine the neighborhood of a cluster center within which
%                   the existence of other cluster centers are discouraged.
%      OPTIONS(2):  The accept ratio, sets the potential, as a fraction of the
%                   potential of the first cluster center, above which
%                   another
%                   data point will be accepted as a cluster center.
%      OPTIONS(3):  The reject ratio sets the potential, as a fraction of the
%                   potential of the first cluster center, below which a data
%                   point will be rejected as a cluster center.
%      OPTIONS(4):  Displays progress information unless it is set to zero.
%   The default values for the OPTIONS vector are [1.25 0.5 0.15 0].
%
%   Examples
%       X1 = 10*rand(50,1);
%       X2 =  5*rand(50,1);
%       X3 = 20*rand(50,1)-10;
%       X = [X1 X2 X3];
%      [C] = subclust(X,0.5) specifies a range of influence of 0.5 for all data
%      dimensions.
%
%      [C] = subclust(X,[0.5 0.25 0.3],[],[2.0 0.8 0.7 0]) specifies a range of
%      influence of 0.5, 0.25, and 0.3 for the first, second and third data
%      dimensions. The scaling factors for mapping the data into a unit hyperbox
%      will be obtained from the minimum and maximum data values. The squash
%      factor is set to 2.0, indicating that we want to only find clusters that
%      are far from each other, the accept ratio is set to 0.8, indicating that
%      we will only accept data points that have very strong potential of being
%      cluster centers, the reject ratio is set to 0.7, indicating that we want
%      to reject all data points without a strong potential.
%
%   See also GENFIS2.

%   Steve Chiu, 1-25-95
%   Roger Jang, 2-26-96, replace some for-loops with vectorized codes
%   N. Hickey, 04-16-01
%   Copyright 1994-2002 The MathWorks, Inc. 
%   $Revision: 1.19 $  $Date: 2002/04/14 22:21:08 $

%	Reference 
%	S. Chiu, "Fuzzy Model Identification Based on Cluster Estimation," J. of
%	Intelligent & Fuzzy Systems, Vol. 2, No. 3, 1994.


[numPoints,numParams] = size(X);

if nargin < 4
	options = [1.25 0.5 0.15 0];
end

if nargin < 3
    xBounds = [];
end

% if only one value is given as the range of influence, then apply that
% value to all data dimensions 
if length(radii) == 1 & numParams ~= 1
    radii = radii * ones(1,numParams);
end

sqshFactor = options(1);
acceptRatio = options(2);
rejectRatio = options(3);
verbose = options(4);

% distance multipliers for accumulating and squashing cluster potentials
accumMultp = 1.0 ./ radii;
sqshMultp = 1.0 ./ (sqshFactor * radii);

if verbose
    disp('Normalizing data...');
end

if length(xBounds) == 0
    % No data scaling range values are specified, use the actual minimum and maximum values
    minX = min(X);
    maxX = max(X);
    % If the actual min and max values have a range of zero, calculate a small range
    % relative to the data, to allow a sigma value to be calculated.
    index = find(maxX == minX);
    minX(index) = minX(index) - 0.0001*(1 + abs(minX(index)));
    maxX(index) = maxX(index) + 0.0001*(1 + abs(maxX(index)));
else
    % Use the user supplied data range values in xBounds
    minX = xBounds(1,:);
    maxX = xBounds(2,:);
    % Verify correct dimensions and values for xBounds were supplied
    if length(minX) ~=  size(X,2)
        error('xBounds contains the wrong dimensions for input data X');
    elseif any(maxX == minX)
        error('xBounds has a data range of zero');
    end
end

% Normalize the data into a unit hyperbox using the verified minX and maxX
for id = 1:numParams
    X(:,id) = (X(:,id) - minX(id)) / (maxX(id) - minX(id));
end
X = min(max(X,0),1);

if verbose
    disp('Computing potential for each data point...');
end

% potVals = the potential of each data point to be a cluster center
potVals = zeros(1,numPoints);

% compute the initial potentials for each data point
%for i=1:numPoints
%	thePoint = X(i,:);
%	% potVals(i) is the potential of the i'th data point
%	potVals(i) = potVals(i) + 1.0;  % add 1.0 for being close to yourself
%
%	for j=(i+1):numPoints
%		nextPoint = X(j,:);
%		dx = (thePoint - nextPoint) .* accumMultp;
%		dxSq = dx * dx';
%		mu = exp(-4.0 * dxSq);
%
%		potVals(i) = potVals(i) + mu;
%		potVals(j) = potVals(j) + mu;
%	end	% endfor j=(i+1):numdata
%
%end	% end for i=1:numdata

% compute the initial potentials for each data point
% Vectorized version by Roger Jang, 2-26-96
new_accumMultp = accumMultp(ones(1,numPoints),:);
UseWaitBar = (numPoints > 500);
if UseWaitBar
	h = waitbar(0,'Please wait...');
end

for i=1:numPoints,
	thePoint = X(i,:);
	thePoint = thePoint(ones(1,numPoints),:);
	dx = (thePoint - X) .* new_accumMultp;
	if numParams == 1,
		potVals(i) = sum(exp(-4*dx.^2));
	else
		potVals(i) = sum(exp(-4*sum(dx.^2,2)));
	end
	if UseWaitBar & mod(i,100)==0
		waitbar(i/numPoints);
	end 
end

if UseWaitBar
	close(h)    
end

% Find the data point with highest potential value.  refPotVal is the
% highest potential value, used as a reference for accepting/rejecting
% other data points as cluster centers.
[refPotVal,maxPotIndex] = max(potVals);

% Start iteratively finding cluster centers and subtracting potential
% from neighboring data points.  maxPotVal is the current highest
% potential value and maxPotIndex is the associated data point's index.
maxPotVal = refPotVal;

% centers = the cluster centers that has been found
centers = [];
numClusters = 0;
findMore = 1;

while findMore & maxPotVal
	findMore = 0;
	maxPoint = X(maxPotIndex,:);
	maxPotRatio = maxPotVal/refPotVal;

	if maxPotRatio > acceptRatio
		% the new peak value is significant, accept it
		findMore = 1;
	elseif maxPotRatio > rejectRatio
		% accept this data point only if it achieves a good balance between having
		% a reasonable potential and being far from all existing cluster centers
		minDistSq = -1;

		for i=1:numClusters
			dx = (maxPoint - centers(i,:)) .* accumMultp;
			dxSq = dx * dx';

			if minDistSq < 0 | dxSq < minDistSq
				minDistSq = dxSq;
			end
		end	% end for i=1:numClusters

		minDist = sqrt(minDistSq);
		if (maxPotRatio + minDist) >= 1
			findMore = 1;	% tentatively accept this data point as a cluster center
		else
			findMore = 2;	% remove this point from further consideration, and continue
		end
	end	% end if maxPotRatio > acceptRatio

	if findMore == 1
		% add the data point to the list of cluster centers
		centers = [centers ; maxPoint];
		numClusters = numClusters + 1;

		if verbose
			msg = sprintf('Found cluster %g, potential = %g',numClusters,maxPotRatio);
			disp(msg);
		end

		% subtract potential from data points near the new cluster center
%		for i=1:numPoints
%			nextPoint = X(i,:);
%			potVal = potVals(i);
%			dx = (maxPoint - nextPoint) .* sqshMultp;
%   			dxSq = dx * dx';
%
%			potVal = potVal - (maxPotVal * exp(-4.0 * dxSq));
%			if potVal < 0
%				potVal = 0;
%			end
%
%			potVals(i) = potVal;
%		end % end for i=1:numdata

		% subtract potential from data points near the new cluster center
		% Vectorized version by Roger Jang, 2-26-96
		new_sqshMultp = sqshMultp(ones(1,numPoints),:);
		tmp = maxPoint(ones(1,numPoints), :);
		dx = (tmp - X) .* new_sqshMultp;
		if numParams == 1,
			deduct = maxPotVal*exp(-4*dx.^2);
		else
			deduct = maxPotVal*exp(-4*sum(dx.^2,2));
		end
		potVals = potVals - deduct';
		potVals(potVals<0) = 0; 

		% find the data point with the highest remaining potential
		[maxPotVal,maxPotIndex] = max(potVals);

	elseif findMore == 2
		potVals(maxPotIndex) = 0;
		[maxPotVal,maxPotIndex] = max(potVals);
	end % end if findMore == 1
end % end while findMore & maxPotVal

% Scale the cluster centers from the normalized values back to values in
% the original range
for i=1:numParams
	centers(:,i) = (centers(:,i) * (maxX(i) - minX(i))) + minX(i);
end

% Compute the sigma values for the clusters
sigmas = (radii .* (maxX - minX)) / sqrt(8.0);

