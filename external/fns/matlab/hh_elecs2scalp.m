function electrodes = hh_elecs2scalp(X,scalp,P,Q)
% hh_elecs2scalp - This function fit the electrode to the scalp
% given the locations of electrode, scalp data structure and the
% matrices P and Q. To find out more information about P and Q
% please consult function fiducial_coregister, which is distributed
% in the Bioelectromagnetism package.
% 
% $Author: Hung Dang$
% $Id: hh_elecs2scalp.m$
% $Log$
% Revision 1.1 Thu Aug 26 17:08:37 MDT 2010, hungptit
% First update and checked.

% Obtain the transformation matrix
T = fiducial_coregister(P,Q);
P2Q = [P,ones(size(P,1),1)] * T;
gridlocs = [X,ones(size(X,1),1)] * T;

% Fit electrode to the scalp
Y = gridlocs([1,4:end],1:3);
Y(:,1) = Y(:,1) * scalp.voxel_sizes(1);
Y(:,2) = Y(:,2) * scalp.voxel_sizes(2);
Y(:,3) = Y(:,3) * scalp.voxel_sizes(3);
N = size(Y,1);
electrodes.locations = zeros(N,3);
electrodes.gridlocs = int32(zeros(N,3));
electrodes.node_sizes = int32(scalp.node_sizes);
electrodes.voxel_sizes = scalp.voxel_sizes;
electrodes.status = int32(ones(N,1));
for n = 1:N
    dx = scalp.locations(:,1) - Y(n,1);
    dy = scalp.locations(:,2) - Y(n,2);
    dz = scalp.locations(:,3) - Y(n,3);
    dis = dx.^2 + dy.^2 + dz.^2;
    pos = find(dis == min(dis));
    electrodes.gridlocs(n,:) = scalp.gridlocs(pos,:);
    electrodes.locations(n,:) = double(scalp.gridlocs(pos,:)) .* ...
        electrodes.voxel_sizes';
    electrodes.status(n) = int32(1);
    electrodes.values(n) = 0.0;
end
