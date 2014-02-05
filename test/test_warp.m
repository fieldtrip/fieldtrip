function test_warp

% MEM 1500mb
% WALLTIME 00:10:00

% function to test the warp functionality and its robustness

% create a set of points
[x,y,z] = ndgrid(1:7,1:10,1:13);
input   = [x(:) y(:) z(:)];
clear x y z;

T = [eye(3) [-1.5 0 2.5]'; 0 0 0 1];
R = [1 0 0 0;0 cos(0.2) -sin(0.2) 0; 0 sin(0.2) cos(0.2) 0; 0 0 0 1];
S = eye(4).*diag([1.1;0.9;1;1]);

target = ft_warp_apply(T*R*S, input);
target = target + 0.1.*randn(size(input,1),3);

global fb;
fb = 1;

try
  [result1, M1] = ft_warp_optim(input, target, 'rigidbody');
catch
  error('rigidbody coregistration failed');
end
try
  [result2, M2] = ft_warp_optim(input, target, 'globalrescale');
catch
  error('rigidbody coregistration failed');
end
try
  [result3, M3] = ft_warp_optim(input, target, 'traditional');
catch
  error('rigidbody coregistration failed');
end
try
  [result4, M4] = ft_warp_optim(input, target, 'nonlin1');
catch
  error('rigidbody coregistration failed');
end
try
  [result5, M5] = ft_warp_optim(input, target, 'nonlin2');
catch
  error('rigidbody coregistration failed');
end
try
  [result6, M6] = ft_warp_optim(input, target, 'nonlin3');
catch
  error('rigidbody coregistration failed');
end
try
  [result7, M7] = ft_warp_optim(input, target, 'nonlin4');
catch
  error('rigidbody coregistration failed');
end
try
  [result8, M8] = ft_warp_optim(input, target, 'nonlin5');
catch
  error('rigidbody coregistration failed');
end
