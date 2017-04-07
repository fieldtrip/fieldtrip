function [EC,ec,degij] = edge_nei_overlap_bu(CIJ)
%EDGE_NEI_OVERLAP_BU        overlap amongst neighbors of two adjacent nodes
%
%   [EC,ec,degij] = edge_nei_bu(CIJ);
%
%   This function determines the neighbors of two nodes that are linked by 
%   an edge, and then computes their overlap.  Connection matrix must be
%   binary and directed.  Entries of 'EC' that are 'inf' indicate that no
%   edge is present.  Entries of 'EC' that are 0 denote "local bridges", i.e.
%   edges that link completely non-overlapping neighborhoods.  Low values
%   of EC indicate edges that are "weak ties".
%
%   If CIJ is weighted, the weights are ignored.
%
%   Inputs:     CIJ,    undirected (binary/weighted) connection matrix
%  
%   Outputs:    EC,     edge neighborhood overlap matrix
%               ec,     edge neighborhood overlap per edge, in vector format
%               degij,  degrees of node pairs connected by each edge
%
%   Reference: Easley and Kleinberg (2010) Networks, Crowds, and Markets. 
%              Cambridge University Press, Chapter 3.
%
%   Olaf Sporns, Indiana University, 2012

[ik,jk,ck] = find(CIJ);
lel = length(ck);
N = size(CIJ,1);

[deg] = degrees_und(CIJ);

ec = zeros(1,lel);
degij = zeros(2,lel);
for e=1:lel
    neiik = setdiff(union(find(CIJ(ik(e),:)),find(CIJ(:,ik(e))')),[ik(e) jk(e)]);
    neijk = setdiff(union(find(CIJ(jk(e),:)),find(CIJ(:,jk(e))')),[ik(e) jk(e)]);
    ec(e) = length(intersect(neiik,neijk))/length(union(neiik,neijk));
    degij(:,e) = [deg(ik(e)) deg(jk(e))];
end;

ff = find(CIJ);
EC = 1./zeros(N);
EC(ff) = ec;                        %#ok<FNDSB>


