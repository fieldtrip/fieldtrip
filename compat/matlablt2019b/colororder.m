function rgb = colororder

% This is a very limited overloaded version of the MathWorks function COLORORDER 
% which has been introduced in MATLAB 2019b. This version does not offer full 
% backwards compatibility, it is only capable of returning the default
% colororder, and not capable of taking any input arguments (as the later
% versions of colororder > 2023b seem to take)

% Copyright (C) 2025 Jan-Mathijs Schoffelen

figure('visible','off');
rgb = get(gca,'colororder');

