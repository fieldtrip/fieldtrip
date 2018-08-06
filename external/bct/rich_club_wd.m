function   [Rw] = rich_club_wd(CIJ,varargin)
%RICH_CLUB_WD	Rich club coefficients curve (weighted directed graph)
%
%   Rw = rich_club_wd(CIJ,varargin)
%
%   The weighted rich club coefficient, Rw, at level k is the fraction of
%   edge weights that connect nodes of degree k or higher out of the
%   maximum edge weights that such nodes might share.
%
%   Inputs:
%       CIJ:        weighted directed connection matrix
%
%       k-level:    (optional) max level of RC(k).
%                   (by default k-level quals the maximal degree of CIJ)
%                
%   Output:
%       Rw:         rich-club curve
%
%
%   References:     
%       T Opsahl et al. Phys Rev Lett, 2008, 101(16)
%       M van den Heuvel, O Sporns, J Neurosci 2011 31(44)
%
%   Martijn van den Heuvel, University Medical Center Utrecht, 2011

%   Modification History:
%   2011: Original
%   2015: Expanded documentation (Mika Rubinov)


NofNodes = size(CIJ,2); %#ok<NASGU>        %number of nodes
NodeDegree = sum((CIJ~=0))+sum((CIJ'~=0)); %define degree of each node (indegree + outdegree)

%define to which level rc should be computed
if size(varargin,2)==1
    klevel = varargin{1};
elseif isempty(varargin)
   klevel = max(NodeDegree);   
else
    error('number of inputs incorrect. Should be [CIJ], or [CIJ, klevel]')
end


%wrank contains the ranked weights of the network, with strongest connections on top

wrank = sort(CIJ(:), 'descend');
    
%loop over all possible k-levels 
for kk = 1:klevel

    SmallNodes=find(NodeDegree<kk);

    if isempty(SmallNodes);
        Rw(kk)=NaN;         %#ok<*AGROW>
        continue
    end
    
    %remove small nodes with NodeDegree<kk
    CutoutCIJ=CIJ;
    CutoutCIJ(SmallNodes,:)=[];
    CutoutCIJ(:,SmallNodes)=[];

    %total weight of connections in subset E>r
    Wr = sum(CutoutCIJ(:));

    %total number of connections in subset E>r
    Er = length(find(CutoutCIJ~=0));

    %E>r number of connections with max weight in network
    wrank_r = wrank(1:1:Er);

    %weighted rich-club coefficient
    Rw(kk)=Wr / sum(wrank_r);

end
 
    
    
    
    
        
