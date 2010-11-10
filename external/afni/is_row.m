function ans = is_row (v);
% ans = is_row (v);
%returns 1 if v is a row vector
%returns 0 if v is a column vector
%returns -1 if v is a matrix
%
% Ziad Saad July 20 97

if (size (v,1) ~= 1 & size (v,2) ~= 1),
	ans = -1;
end;

if (size (v,1) == 1),
	ans = 1;
	elseif (size (v,2) == 1),
		ans = 0;
end;

return;
