function openmeeg_license

% OPENMEEG_LICENSE prints the license only once upon the first call to this
% function. If the user does a "clear all", the license will again be shown.
% This function should be included in every openmeeg function to ensure
% that the license is displayed at least once.

% Copyright (C) 2009, Robert Oostenveld
%
% $Log: openmeeg_license.m,v $
% Revision 1.3  2009/08/01 11:51:22  alegra
% adding license and bibtex citation
%
% Revision 1.2  2009/03/30 15:16:39  roboos
% show the content of the COPYING file
%
% Revision 1.1  2009/03/27 14:45:24  roboos
% added first skeleton and license function
%

persistent status

if isempty(status)
  % this is for the first time
  status = false;
end

if ~status
  % show the license
  [p, f] = fileparts(mfilename);
  clc
  fprintf('==============================================================================\n');
  type(fullfile(p, 'COPYING'));
  fprintf('==============================================================================\n');
  % remember that it has been shown
  status = true;
  pause(3)
end

