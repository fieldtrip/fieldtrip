function list = peerschedule(list, memreq, timreq, cpureq)

% PEERSCHEDULE sorts the list of avaialble peers according to a number of heuristic 
% reules that try to optimize the use of the available resources.
%
% Use as
%   list = peerschedule(list, memreq, timreq, cpureq)

% Copyright (C) 2011, Robert Oostenveld
%
% $Id$

% the first penalty measure is based on the excess memory
memavail = [list.memavail];
penalty = memavail - memreq;

% increase the penalty for slaves with less than 1GB
penalty = penalty + (memavail<1e9)*1e9;

% select the slave peer that has the best match with the job requirements
[penalty, indx] = sort(penalty);

% sort the list according to the penalty
list = list(indx);
