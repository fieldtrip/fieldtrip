function res = afni_isdigit(s)
% res = afni_isdigit(s);
% res(i) = 1 is s(1) is a digit between 0 and 9 (inclusive)
%
digit_msk=s >= '0' & s <= '9';
res = zeros(size(s));
res(digit_msk) = 1;

return;
