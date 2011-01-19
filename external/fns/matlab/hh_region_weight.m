function G = hh_region_weight(data,compress,node_sizes,dipoles,weight)
% hh_region_weight - This function finds the weight matrix for the given
% brain region
%
% Usage: G = hh_region_weight(data,compress,node_sizes,region)
% where 
%   data        is the reciprocity data matrix
%   compress    is the compress vector
%   node_sizes  is the model grid sizes
%   gridlocs    is a Nx3 matrix which keeps the grid locations of nodes
%   belong to the region of interrest
%
% $Author: Hung Dang$
% $Id: hh_region_weight.m$
% 
% Revision 1.0 Thu Aug  5 14:01:36 MDT 2010, hungptit
% Created to test the fMRI weighted beamformer
%
% Revision 1.1 Thu Aug 19 02:07:43 MDT 2010, hungptit
% Dipoles is a Nx3 integer matrix which stores the grid locations
% of dipoles.
%
% $$$ Update parameters
K = size(dipoles,1); % The number of nodes belong to the region of interrest
N = size(data,1);   % This is the number of sensors

% $$$ Preallocate memory for G
G = zeros(N,3*K);    

% $$$ Obtain the lead field matrix G for the region of interrest
if isempty(weight)
    for k = 1:K
        L = hh_leadfield(data,compress,node_sizes,...
                         dipoles(k,1), ...
                         dipoles(k,2),...
                         dipoles(k,3));
        if ~isempty(L)
            G(:,(3*k - 2):(3*k)) = L;
        end
    end
else
    for k = 1:K
        L = hh_leadfield(data,compress,node_sizes,...
                         dipoles(k,1), ...
                         dipoles(k,2),...
                         dipoles(k,3));
        if ~isempty(L)
            G(:,(3*k - 2):(3*k)) = L * weight(k);
        end
    end
end