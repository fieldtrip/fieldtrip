function [d] = inv2x2(x)

%computes inverse of matrix x, using explicit analytic definition if
%size(x,1) < 4, otherwise use matlab inv-function

siz = size(x);
if all(siz(1:2)==2),
  adjx  = [x(2,2,:,:) -x(1,2,:,:); -x(2,1,:,:) x(1,1,:,:)];
  denom = det2x2(x);
  d     = adjx./denom([1 1],[1 1],:,:);
elseif all(siz(1:2)==3),
  adjx = [ det2x2(x([2 3],[2 3],:,:)) -det2x2(x([1 3],[2 3],:,:))  det2x2(x([1 2],[2 3],:,:)); ...
          -det2x2(x([2 3],[1 3],:,:))  det2x2(x([1 3],[1 3],:,:)) -det2x2(x([1 2],[1 3],:,:)); ...
	   det2x2(x([2 3],[1 2],:,:)) -det2x2(x([1 3],[1 2],:,:))  det2x2(x([1 2],[1 2],:,:))];
  denom = det2x2(x);
  d     = adjx./denom([1 1 1],[1 1 1],:,:);
elseif numel(siz)==2,
  d = inv(x);
else
  error('cannot compute slicewise inverse');
  %write for loop for the higher dimensions, using normal inv
end
