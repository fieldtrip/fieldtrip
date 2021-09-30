function rgb = skin_dark

% This returns a predefined color as [red green blue] values
%   red               = [255   0   0]/255;
%   green             = [  0 192   0]/255;
%   blue              = [  0   0 255]/255;
%   magenta           = [255 255   0]/255;
%   cyan              = [  0 255 255]/255;
%   yellow            = [255 255   0]/255;
%   white             = [255 255 255]/255;
%   black             = [  0   0   0]/255;
%
%   skull             = [140  85  85]/255
%   cortex            = [255 213 119]/255;
%   cortex_light      = [199 194 169]/255;
%   cortex_dark       = [100  97  85]/255;
%   skin              = [249 223 192]/255;
%   skin_light        = [249 223 192]/255;
%   skin_medium_light = [225 194 158]/255;
%   skin_medium       = [188 142 106]/255;
%   skin_medium_dark  = [155 102	65]/255;
%   skin_dark         = [ 91  71  61]/255;
%
% The different skin-based colors follow the Fitzpatrick scale with type I and II
% combined, and return RGB values that approximate those used by Apple in the emoji
% skin tones. See also https://emojipedia.org/emoji-modifier-sequence/
%
% If no specific skin tone is specified, this function returns a light skin color.
% This corresponds with that of one of the developers who approximated his own skin
% color more than 15 years ago upon the first implementation of this function.

% returns a predefined color as [red green blue] values
rgb = [ 91	 71	 61]/255;