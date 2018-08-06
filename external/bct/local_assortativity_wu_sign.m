function [loc_assort_pos,loc_assort_neg] = local_assortativity_wu_sign(W)
%LOCAL_ASSORTATIVITY_WU_SIGN     Local Assortativity
%
%   [loc_assort_pos,loc_assort_neg] = local_assortativity_wu_sign(W);
%
%   Local Assortativity measures the extent to which nodes are connected to
%   nodes of similar strength (vs. higher or lower strength). Adapted from
%   Thedchanamoorthy et al. (2014)'s formula to allow weighted/signed 
%   networks (node degree replaced with node strength). Note, output values 
%   sum to total assortativity. 
%
%   Inputs:     W,        undirected connection matrix with positive and
%                         negative weights
%
%   Output:     loc_assort_pos, local assortativity from positive weights
%
%               loc_assort_neg, local assortativity from negative weights
%
%   Reference: Thedchanamoorthy G, Piraveenan M, Kasthuriratna D, 
%              Senanayake U. Proc Comp Sci (2014) 29:2449-2461.
%
%
%   Jeff Spielberg, Boston University

%   Modification History:
%   May 2015: Original

W(1:(size(W,1)+1):end) = 0;
r_pos = assortativity_wei(W.*(W>0),0);
r_neg = assortativity_wei(-W.*(W<0),0);
[str_pos,str_neg] = strengths_und_sign(W);
loc_assort_pos = nan(size(W,1),1);
loc_assort_neg = nan(size(W,1),1);

for curr_node = 1:size(W,1)
    [~,j_pos] = find(W(curr_node,:)>0);
    loc_assort_pos(curr_node,1) = sum(abs(str_pos(j_pos)-str_pos(curr_node)))/str_pos(curr_node);
    
    [~,j_neg] = find(W(curr_node,:)<0);
    loc_assort_neg(curr_node,1) = sum(abs(str_neg(j_neg)-str_neg(curr_node)))/str_neg(curr_node);
end

loc_assort_pos = ((r_pos+1)/size(W,1))-(loc_assort_pos/sum(loc_assort_pos));
loc_assort_neg = ((r_neg+1)/size(W,1))-(loc_assort_neg/sum(loc_assort_neg));
