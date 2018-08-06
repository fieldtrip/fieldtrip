function  [fc,FC,total_flo] = flow_coef_bd(CIJ)
%FLOW_COEF_BD        Node-wise flow coefficients
%
%   [hc,HC,total_flo] = flow_coef_bd(CIJ)
%
%   Computes the flow coefficient for each node and averaged over the
%   network, as described in Honey et al. (2007) PNAS. The flow coefficient
%   is similar to betweenness centrality, but works on a local
%   neighborhood. It is mathematically related to the clustering
%   coefficient  (cc) at each node as, fc+cc <= 1.
%
%   input:      CIJ,	connection/adjacency matrix (binary, directed)
%   output:     fc,     flow coefficient for each node
%               FC,     average flow coefficient over the network
%        total_flo,     number of paths that "flow" across the central node
%
%   Reference:  Honey et al. (2007) Proc Natl Acad Sci U S A
%
%   Olaf Sporns, Indiana University, 2007/2010/2012

N = size(CIJ,1);

% initialize ...
fc        = zeros(1,N);
total_flo = fc;
max_flo   = fc;

% loop over nodes
for v=1:N
    % find neighbors - note: treats incoming and outgoing connections as equal
    [nb] = find(CIJ(v,:) + CIJ(:,v)');
    fc(v) = 0;
    if (~isempty(nb))
        CIJflo = -CIJ(nb,nb);
        for i=1:length(nb)
            for j=1:length(nb)
                if((CIJ(nb(i),v))==1)&&(CIJ(v,nb(j))==1);
                    CIJflo(i,j) = CIJflo(i,j) + 1;
                end;
            end;
        end;
        total_flo(v) = sum(sum(double(CIJflo==1).*~eye(length(nb))));
        max_flo(v) = length(nb)^2-length(nb);
        fc(v) = total_flo(v)/max_flo(v);
    end;
end;

% handle nodes that are NaNs
fc(isnan(fc)) = 0;

FC = mean(fc);
