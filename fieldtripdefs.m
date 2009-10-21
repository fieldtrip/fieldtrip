function fieldtripdefs

% FIELDTRIPDEFS is called at the begin of all FieldTrip functions and
% contains some defaults and path settings
%
% Note that this should be a function and not a script, otherwise the
% hastoolbox function appears not be found in fieldtrip/private.

% $Log: fieldtripdefs.m,v $
% Revision 1.14  2009/04/02 19:46:24  roboos
% added plotting module
%
% Revision 1.13  2009/03/26 09:28:17  roboos
% removed fixpath from the matlab path
%
% Revision 1.12  2008/11/12 11:41:54  sashae
% changed trackconfig options
%
% Revision 1.11  2008/11/11 13:13:44  roboos
% added/improved size checking
%
% Revision 1.10  2008/11/10 21:17:28  roboos
% added public subdirectory
%
% Revision 1.9  2008/10/08 10:21:57  roboos
% added realtime module for addpath
%
% Revision 1.8  2008/10/02 14:04:17  roboos
% added global ft_default, see also checkconfig
%
% Revision 1.7  2008/10/01 10:18:30  roboos
% changed from script into function to fix private path
%
% Revision 1.6  2008/10/01 08:22:40  roboos
% added specest module
%
% Revision 1.5  2008/09/23 07:54:09  roboos
% added template directory, moved all layouts and cortical meshes there
%
% Revision 1.4  2008/09/23 07:40:10  roboos
% changed order of the path additioons, fixpath should be the first (i.e. lowest on the path)
%
% Revision 1.3  2008/09/23 07:31:13  roboos
% Moved all statfuns and trialfuns to their own directories, where they will be easier to find for the end-user. Also updated fieldtripdefs accordingly.
%
% Revision 1.2  2008/09/22 20:58:57  roboos
% added modules: fileio, forwinv, preproc
%
% Revision 1.1  2008/09/22 20:02:50  roboos
% initial version, contains fixpath
%


% set the global defaults, the checkconfig function will copy these into the local configurations
global ft_default
if ~isfield(ft_default, 'trackconfig'), ft_default.trackconfig = 'off';   end % cleanup, report, off
if ~isfield(ft_default, 'checkconfig'), ft_default.checkconfig = 'loose'; end % pedantic, loose, silent
if ~isfield(ft_default, 'checksize'),   ft_default.checksize   = 1e5;     end % number in bytes, can be inf

try
  % this contains general usefull functions, that do not have to be private
  hastoolbox('public', 1, 1);
end

try
  % this contains layouts and cortical meshes
  hastoolbox('template', 1, 1);
end

try
  % this is used in statistics
  hastoolbox('statfun', 1, 1);
end

try
  % this is used in definetrial
  hastoolbox('trialfun', 1, 1);
end

try
  % numerous functions depend on this module
  hastoolbox('fileio', 1, 1);
end

try
  % numerous functions depend on this module
  hastoolbox('forwinv', 1, 1);
end

try
  % numerous functions depend on this module
  hastoolbox('preproc', 1, 1);
end

try
  % numerous functions depend on this module
  hastoolbox('plotting', 1, 1);
end

try
  % only very few functions depent on this module
  % it is not yet included in the FTP release version
  hastoolbox('specest', 1, 1);
end

try
  % this contains some examples for realtime processing
  hastoolbox('realtime', 1, 1);
end

