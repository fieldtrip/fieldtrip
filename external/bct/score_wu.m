function [CIJscore,sn] = score_wu(CIJ,s)
%SCORE_WU       S-score
%
%   [CIJscore,sn] = score_wu(CIJ,s);
%
%   The s-core is the largest subnetwork comprising nodes of strength at
%   least s. This function computes the s-core for a given weighted
%   undirected connection matrix. Computation is analogous to the more
%   widely used k-core, but is based on node strengths instead of node
%   degrees. 
%
%   input:          CIJ,	connection/adjacency matrix (weighted, undirected)
%                     s,    level of s-core. Note: s can take on any fractional value
%
%   output:    CIJscore,    connection matrix of the s-core.  This matrix 
%                           contains only nodes with a strength of at least s.
%                    sn,    size of s-score
%
%   Olaf Sporns, Indiana University, 2007/2008/2010/2012

while 1

    % get strengths of matrix
    [str] = strengths_und(CIJ);

    % find nodes with strength <s
    ff = find((str<s)&(str>0));
    
    % if none found -> stop
    if (isempty(ff)) break; end;            %#ok<SEPEX>

    % peel found nodes
    CIJ(ff,:) = 0;
    CIJ(:,ff) = 0;

end;

CIJscore = CIJ;
sn = sum(str>0);

