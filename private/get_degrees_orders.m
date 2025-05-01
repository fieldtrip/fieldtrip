function [degrees, orders] = get_degrees_orders(order)

% this helper function assumes the order of the basis vectors, as computed
% in the MNE-Python implementation is the same as the order of the basis
% vectors as computed in the SPM implementation

% def _get_degrees_orders(order):
%     """Get the set of degrees used in our basis functions."""

degrees = zeros(1,0);
orders  = zeros(1,0);
for k = 1:order
  orders  = cat(2, orders, -k:k);
  degrees = cat(2, degrees, ones(1,2.*k+1).*k); 
end
