function  s = strtrim_improve(s0)
%             
% strtrim handles only a limited number of whitespace types.
% But does NOT handle 0 for example. This function turns all
% whitespaces including 0 to plain old spaces which strtrim
% does handle.
%
%
s = s0;
s(s<33) = ' ';     % Replace all non-printable chracters with spaces
    
% Feed the digestable new string to strtrim.
s = strtrim(s);


