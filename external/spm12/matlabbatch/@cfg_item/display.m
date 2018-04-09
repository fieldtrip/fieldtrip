function display(item)

% function display(item)
% Display a configuration object
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: display.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

disp(' ');
disp([inputname(1),' = ']);
disp(' ');
disp(item);
disp(' ');
