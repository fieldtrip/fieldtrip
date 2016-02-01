function dist = dist_euclidean(metric, x)
% DIST_EUCLIDEAN   is a mex file that returns the euclidean distance matrix
% for metric_euclidean.  If we end up here the mex file has not been
% compiled. In this case the function returns NaN and the distance is
% evaluated in metric_euclidean.

% Copyright (c) 2008 Jarno Vanhatalo     

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

dist = NaN;