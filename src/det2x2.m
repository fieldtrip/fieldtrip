function [d] = det2x2(x)

%computes determinant of matrix x, using explicit analytic definition if
%size(x,1) < 4, otherwise use matlab det-function

siz = size(x);
if all(siz(1:2)==2),
  d = x(1,1,:,:).*x(2,2,:,:) - x(1,2,:,:).*x(2,1,:,:);
elseif all(siz(1:2)==3),
  d = x(1,1,:,:).*x(2,2,:,:).*x(3,3,:,:) - ...
      x(1,1,:,:).*x(2,3,:,:).*x(3,2,:,:) - ...
      x(1,2,:,:).*x(2,1,:,:).*x(3,3,:,:) + ...
      x(1,2,:,:).*x(2,3,:,:).*x(3,1,:,:) + ...
      x(1,3,:,:).*x(2,1,:,:).*x(3,2,:,:) - ...
      x(1,3,:,:).*x(2,2,:,:).*x(3,1,:,:);
elseif numel(siz)==2,
  d = det(x);
else
  %error   
  %write for loop
  %for
  %end
end
