function   r = assortativity_bin(CIJ,flag)
%ASSORTATIVITY      Assortativity coefficient
%
%   r = assortativity(CIJ,flag);
%
%   The assortativity coefficient is a correlation coefficient between the
%   degrees of all nodes on two opposite ends of a link. A positive
%   assortativity coefficient indicates that nodes tend to link to other
%   nodes with the same or similar degree.
%
%   Inputs:     CIJ,    binary directed/undirected connection matrix
%               flag,   0, undirected graph: degree/degree correlation
%                       1, directed graph: out-degree/in-degree correlation
%                       2, directed graph: in-degree/out-degree correlation
%                       3, directed graph: out-degree/out-degree correlation
%                       4, directed graph: in-degree/in-degree correlation
%
%   Outputs:    r,      assortativity coefficient
%
%   Notes: The function accepts weighted networks, but all connection
%   weights are ignored. The main diagonal should be empty. For flag 1
%   the function computes the directed assortativity described in Rubinov
%   and Sporns (2010) NeuroImage.
%
%   Reference:  Newman (2002) Phys Rev Lett 89:208701
%               Foster et al. (2010) PNAS 107:10815–10820
%
%   Olaf Sporns, Indiana University, 2007/2008
%   Vassilis Tsiaras, University of Crete, 2009
%   Murray Shanahan, Imperial College London, 2012
%   Mika Rubinov, University of Cambridge, 2012

if (flag==0)                        % undirected version
    deg = degrees_und(CIJ);
    [i,j] = find(triu(CIJ,1)>0);
    K = length(i);
    degi = deg(i);
    degj = deg(j);

else                                % directed versions
    [id,od] = degrees_dir(CIJ);
    [i,j] = find(CIJ>0);
    K = length(i);

    switch flag
        case 1
            degi = od(i);
            degj = id(j);
        case 2
            degi = id(i);
            degj = od(j);
        case 3
            degi = od(i);
            degj = od(j);
        case 4
            degi = id(i);
            degj = id(j);
    end
end

% compute assortativity
r = ( sum(degi.*degj)/K - (sum(0.5*(degi+degj))/K)^2 ) / ...
    ( sum(0.5*(degi.^2+degj.^2))/K - (sum(0.5*(degi+degj))/K)^2 );

