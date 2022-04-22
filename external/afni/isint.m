function [res, int_part, dec_part] = isint (v)
%
%   [res, int_part, dec_part] = isint (v)
%
%Purpose:
%   obvious
%
%
%Input Parameters:
%
%
%
%Output Parameters:
%   res : 0 it's not an integer
%       : 1
%
%
%
%Key Terms:
%
%More Info :
%
%   [res, int_part, dec_part] = isint (4.07)
%   [res, int_part, dec_part] = isint (4.00)
%   [res, int_part, dec_part] = isint (-4.07)
%   [res, int_part, dec_part] = isint (-4.00)
%
%     Author : Ziad Saad
%     Date : Wed Aug 11 16:42:55 CDT 1999


%Define the function name for easy referencing
FuncName = 'isint';

%Debug Flag
DBG = 1;

%initailize return variables
fv = fix(v);
df = v - fv;
if (~df),
	res = 1;
	dec_part = 0;
	int_part = v;
else
	res = 0;
	dec_part = df.*sign(v);
	int_part = fv.*sign(v);
end





err = 0;
return;

