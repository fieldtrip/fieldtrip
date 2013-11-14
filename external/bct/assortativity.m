function   r = assortativity(CIJ,flag)
%ASSORTATIVITY      Assortativity coefficient
%
%   r = assortativity(CIJ,flag);
%
%   The assortativity coefficient is a correlation coefficient between the 
%   degrees of all nodes on two opposite ends of a link. A positive 
%   assortativity coefficient indicates that nodes tend to link to other 
%   nodes with the same or similar degree.
%
%   Inputs:     CIJ,        binary directed/undirected connection matrix
%               flag,       1 = directed graph; 0 = non-directed graph
%
%   Outputs:    r,          assortativity
%
%   Notes: The function accepts weighted networks, but all connection
%   weights are ignored. The main diagonal should be empty.
%
%   Reference:  Newman (2002) Phys Rev Lett 89:208701.
%
%
%   Olaf Sporns, Indiana University, 2007/2008
%   Vassilis Tsiaras, University of Crete, 2009

if (flag==0)
    [deg] = degrees_und(CIJ);
    [i,j] = find(triu(CIJ,1)>0);
    K = length(i);
    for k=1:K
        degi(k) = deg(i(k));
        degj(k) = deg(j(k));
    end;
end;
if (flag==1)
    [id,od,deg] = degrees_dir(CIJ);
    [i,j] = find(CIJ>0);
    K = length(i);
    for k=1:K
        degi(k) = deg(i(k));
        degj(k) = deg(j(k));
    end;
end;

% compute assortativity
r = (sum(degi.*degj)/K - (sum(0.5*(degi+degj))/K)^2)/(sum(0.5*(degi.^2+degj.^2))/K - (sum(0.5*(degi+degj))/K)^2);

