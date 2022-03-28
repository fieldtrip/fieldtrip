function [ant] = iseven (x);
% [ans] = ISEVEN (x)
% returns 1 if x is even and 0 if x is odd
%
%			Ziad Saad May 8 95.


	d = x ./ 2;
	if fix(d) == d,
	  ant = 1;
 	else
		ant = 0;
	end;


