function res = isdigit(s)
% res = isdigit(s);
% res(i) = 1 is s(1) is a digit between 0 and 9 (inclusive)
%
id = find(s >= '0' & s <= '9');
res = zeros(size(s));
res(id) = 1;

return;
