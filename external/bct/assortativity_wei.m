function   r = assortativity_wei(CIJ,flag)
%ASSORTATIVITY      Assortativity coefficient
%
%   r = assortativity_wei(CIJ,flag);
%
%   The assortativity coefficient is a correlation coefficient between the
%   strengths (weighted degrees) of all nodes on two opposite ends of a link.
%   A positive assortativity coefficient indicates that nodes tend to link to
%   other nodes with the same or similar strength.
%
%   Inputs:     CIJ,    weighted directed/undirected connection matrix
%               flag,   0, undirected graph: strength/strength correlation
%                       1, directed graph: out-strength/in-strength correlation
%                       2, directed graph: in-strength/out-strength correlation
%                       3, directed graph: out-strength/out-strength correlation
%                       4, directed graph: in-strength/in-strength correlation
%
%   Outputs:    r,      assortativity coefficient
%
%   Notes: The main diagonal should be empty. For flag 1 the function computes 
%   the directed assortativity described in Rubinov and Sporns (2010) NeuroImage.
%
%   Reference:  Newman (2002) Phys Rev Lett 89:208701
%               Foster et al. (2010) PNAS 107:10815-10820
%
%   Olaf Sporns, Indiana University, 2007/2008
%   Vassilis Tsiaras, University of Crete, 2009
%   Murray Shanahan, Imperial College London, 2012
%   Mika Rubinov, University of Cambridge, 2012

if (flag==0)                        % undirected version
    str = strengths_und(CIJ);
    [i,j] = find(triu(CIJ,1)>0);
    K = length(i);
    stri = str(i);
    strj = str(j);

else                                % directed versions
    [is,os] = strengths_dir(CIJ);
    [i,j] = find(CIJ>0);
    K = length(i);

    switch flag
        case 1
            stri = os(i);
            strj = is(j);
        case 2
            stri = is(i);
            strj = os(j);
        case 3
            stri = os(i);
            strj = os(j);
        case 4
            stri = is(i);
            strj = is(j);
    end
end

% compute assortativity
r = ( sum(stri.*strj)/K - (sum(0.5*(stri+strj))/K)^2 ) / ...
    ( sum(0.5*(stri.^2+strj.^2))/K - (sum(0.5*(stri+strj))/K)^2 );

