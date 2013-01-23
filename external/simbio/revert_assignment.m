function [node_ele_assignment] = revert_assignment(elements,sparse_flag)

% REVERT_ASSIGNMENT
%
% $Id$
    
N_nodes = max(max(elements));
N_elements = max(size(elements));
[sort_elements_vec, sort_ind] = sort(elements(:));
[dummy, m, dummy] = unique(sort_elements_vec,'last');
occurence = diff([0;m]);

node_ele_assignment = zeros(max(occurence),N_nodes);
sort_ind = mod(sort_ind,N_elements);
sort_ind(sort_ind == 0) = N_elements;
flip = diff([sort_elements_vec;sort_elements_vec(end)]);
col = 1;
row = 1;

for i=1:length(flip)
    node_ele_assignment(row,col) = sort_ind(i);
    row = row+1;
    if(flip(i))
        col = col+1;
        row = 1;
    end
end 


s1 = whos('node_ele_assignment');
sp = sparse(node_ele_assignment);
s2 = whos('sp');
if(s2.bytes/s1.bytes < sparse_flag)
    node_ele_assignment = sp;
end

end
