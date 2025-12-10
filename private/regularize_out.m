function [out_remove] = regularize_out(int_order, ext_order, S, ismag, extended_remove)

if nargin<4 || isempty(ismag)
  ismag = true(size(S,1),1);
end
if nargin<5 || isempty(extended_remove)
  extended_remove = [];
end

n_in  = (int_order+2).*int_order;
remove_homog = ext_order>0 && ~any(ismag);
if remove_homog
  out_remove = [n_in + (1:3) extended_remove];
else
  out_remove = extended_remove;
end
