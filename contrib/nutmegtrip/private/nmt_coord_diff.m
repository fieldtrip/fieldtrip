function C = nmt_coord_diff(A,B)
% C = NMT_COORD_DIFF(A,B)
%
% subtract single 1 x 3 coordinate B from list of coordinates A (N x 3)
if size(B,1)>1
    error('size of second input incorrect')
    return
end

C = [A(:,1)-B(1) A(:,2)-B(2) A(:,3)-B(3)];
