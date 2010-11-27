function S = pad_strn (strg,padchar,n,loc)
%	S = PAD_STRN (strg,padchar,n,loc)
%This function pad the string 'strg' with 'padchar' characters
%  such that 'S' is 'n' characters long.
% if loc = 1 the 'padchar' are added to the left of 'strg'
% if loc = -1 the 'padchar' are added to the right of 'strg'


if (~ischar(strg) | ~ischar(padchar)),
	err = ErrEval(mfilename, 'Err_First two arguments must be characters');
	S = '';
	return;
end

S = strg;
while (length(S) < n)

	if (loc == 1),
		S = sprintf ('%s%s',padchar,S);
	elseif (loc == -1),
		S = sprintf ('%s%s',S,padchar);
	end

end;

return;
