function ft_defaults

% FT_DEFAULTS (ending with "s") sets some general settings in the global variable
% ft_default (without the "s") and takes care of the required path settings. You can
% call this function in your startup.m script. This function is also called at the
% begin of all FieldTrip functions.
%
% The global configuration defaults are stored in the global "ft_default" structure.
% The ft_checkconfig function that is called by many FieldTrip functions will merge
% these global configuration defaults with the cfg ctructure that you pass to
% the FieldTrip function that you are calling.
%
% The global options and their default values are
%   ft_default.checkconfig       = string, can be 'pedantic', 'loose', 'silent' (default = 'loose')
%   ft_default.checkpath         = string, can be 'pedantic', 'once', 'no' (default = 'pedantic')
%   ft_default.checksize         = number in bytes, can be inf (default = 1e5)
%   ft_default.showcallinfo      = string, can be 'yes' or 'no' (default = 'yes')
%   ft_default.trackcallinfo     = string, can be 'yes' or 'no' (default = 'yes')
%   ft_default.trackconfig       = string, can be 'cleanup', 'report', 'off' (default = 'off')
%   ft_default.trackusage        = false, or string with salt for one-way encryption of identifying information (by default this is enabled and an automatic salt is created)
%   ft_default.trackdatainfo     = string, can be 'yes' or 'no' (default = 'no')
%   ft_default.outputfilepresent = string, can be 'keep', 'overwrite', 'error' (default = 'overwrite')
%   ft_default.debug             = string, can be 'display', 'displayonerror', 'displayonsuccess', 'save', 'saveonerror', saveonsuccess' or 'no' (default = 'no')
%   ft_default.toolbox.signal    = string, can be 'compat' or 'matlab' (default = 'compat')
%   ft_default.toolbox.stats     = string, can be 'compat' or 'matlab' (default = 'compat')
%   ft_default.toolbox.images    = string, can be 'compat' or 'matlab' (default = 'compat')
%   ft_default.reproducescript   = string, directory to which the script and intermediate data are written (default = [])
%
% If you want to overrule these default settings, you can add something like this in your startup.m script
%   ft_defaults
%   global ft_default
%   ft_default.option1 = value1
%   ft_default.option2 = value2
%
% The toolbox option for signal, stats and images allows you to specify whether you
% want to use a compatible drop-in to be used for these MathWorks toolboxes, or the
% original version from MathWorks.  The default is 'compat', which has the advantage
% that you do not need a license for these toolboxes.
%
% See also FT_HASTOOLBOX, FT_CHECKCONFIG, FT_TRACKUSAGE

% undocumented options
%   ft_default.siunits        = 'yes' or 'no'

% Copyright (C) 2009-2022, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

global ft_default

persistent initialized
persistent checkpath

if isempty(initialized)
  initialized = false;
end

if isempty(checkpath)
  checkpath = false;
end

% ft_warning is located in fieldtrip/utilities, which may not be on the path yet
if ~exist('ft_warning', 'file')
  ft_warning = @warning;
end

% locate the file with the persistent FieldTrip preferences
fieldtripprefs = fullfile(prefdir, 'fieldtripprefs.mat');
if exist(fieldtripprefs, 'file')
  prefs       = load(fieldtripprefs); % the file contains multiple fields
  ft_default  = mergeconfig(ft_default, prefs);
end

% Set the defaults in a global variable, ft_checkconfig will copy these over into the local configuration.
% NOTE ft_getopt might not be available on the path at this moment and can therefore not yet be used.
% NOTE all options here should be explicitly listed as allowed in ft_checkconfig

if ~isfield(ft_default, 'trackconfig'),       ft_default.trackconfig    = 'off';      end % cleanup, report, off
if ~isfield(ft_default, 'checkconfig'),       ft_default.checkconfig    = 'loose';    end % pedantic, loose, silent
if ~isfield(ft_default, 'checkpath'),         ft_default.checkpath      = 'pedantic'; end % pedantic, once, no
if ~isfield(ft_default, 'checksize'),         ft_default.checksize      = 1e5;        end % number in bytes, can be inf
if ~isfield(ft_default, 'showcallinfo'),      ft_default.showcallinfo   = 'yes';      end % yes or no, this is used in ft_pre/postamble_provenance
if ~isfield(ft_default, 'debug'),             ft_default.debug          = 'no';       end % no, save, saveonerror, display, displayonerror, this is used in ft_pre/postamble_debug
if ~isfield(ft_default, 'outputfilepresent'), ft_default.outputfilepresent = 'overwrite'; end % can be keep, overwrite, error

% These options allow to disable parts of the provenance
if ~isfield(ft_default, 'trackcallinfo'),  ft_default.trackcallinfo  = 'yes';    end % yes or no
if ~isfield(ft_default, 'trackdatainfo'),  ft_default.trackdatainfo  = 'no';     end % yes or no
if ~isfield(ft_default, 'tracktimeinfo'),  ft_default.tracktimeinfo  = 'yes';    end % yes or no
if ~isfield(ft_default, 'trackmeminfo')
  if ispc()
    % don't track memory usage info under Windows; this does not work (yet)
    ft_default.trackmeminfo = 'no';
  else
    ft_default.trackmeminfo = 'yes';
  end
end

% These options allow to prefer the MATLAB toolbox implementations ('matlab') over the drop-in replacements ('compat').
if ~isfield(ft_default, 'toolbox'), ft_default.toolbox  = []; end
if ~isfield(ft_default.toolbox, 'images'), ft_default.toolbox.images  = 'compat'; end % matlab or compat
if ~isfield(ft_default.toolbox, 'stats') , ft_default.toolbox.stats   = 'compat'; end % matlab or compat
if ~isfield(ft_default.toolbox, 'signal'), ft_default.toolbox.signal  = 'compat'; end % matlab or compat

% Some people mess up their path settings and then have stuff on the path that should not be there.
% The following will issue a warning
switch ft_default.checkpath
  case 'pedantic'
    % check every time
    checkIncorrectPath();
  case 'once'
    % check only once
    if ~checkpath
      checkIncorrectPath();
      checkpath = true;
    end
  case 'no'
    % do not check
end % case

% Check whether this ft_defaults function has already been executed. Note that we
% should not use ft_default itself directly, because the user might have set stuff
% in that struct already prior to ft_defaults being called for the first time.
if initialized && exist('ft_hastoolbox', 'file')
  return;
end

if isfield(ft_default, 'toolbox') && isfield(ft_default.toolbox, 'cleanup')
  prevcleanup = ft_default.toolbox.cleanup;
else
  prevcleanup = {};
end

% show a welcome message
fprintf([
  '-------------------------------------------------------------------------------------------' ...
  '\nFieldTrip is developed by members and collaborators of the Donders Institute for Brain,' ...
  '\nCognition and Behaviour at Radboud University, Nijmegen, the Netherlands.' ...
  '\n'...
  '\n                          --------------------------'...
  '\n                        /                            \\'...
  '\n                     ------------------------------------'...
  '\n                    /                                    \\'...
  '\n          -------------------------------------------------'...
  '\n         /                            /\\/\\/\\/\\/\\ '...
  '\n         ---------------------------------------------------'...
  '\n                  |        F  i  e  l  d  T  r  i  p       |'...
  '\n                  ------------------------------------------'...
  '\n                   \\                                      /'...
  '\n                     ------------------------------------'...
  '\n                          \\            /'...
  '\n                            ----------'...
  '\n' ...
  '\nPlease cite the FieldTrip reference paper when you have used FieldTrip in your study.' ...
  '\nRobert Oostenveld, Pascal Fries, Eric Maris, and Jan-Mathijs Schoffelen. FieldTrip: Open' ...
  '\nSource Software for Advanced Analysis of MEG, EEG, and Invasive Electrophysiological Data.' ...
  '\nComputational Intelligence and Neuroscience, vol. 2011, Article ID 156869, 9 pages, 2011.' ...
  '\ndoi:10.1155/2011/156869.' ...
  '\n-------------------------------------------------------------------------------------------\n'
  ]);

% Ensure that the path containing ft_defaults is on the path.
% This allows people to do "cd path_to_fieldtrip; ft_defaults"
ftPath = fileparts(mfilename('fullpath')); % get the full path to this function, strip away 'ft_defaults'
ftPath = strrep(ftPath, '\', '\\');
if isempty(regexp(path, [ftPath pathsep '|' ftPath '$'], 'once'))
  ft_warning('FieldTrip is not yet on your MATLAB path, adding %s', strrep(ftPath, '\\', '\'));
  addpath(ftPath);
end

if ~isdeployed
  
  if isempty(which('ft_hastoolbox')) || isempty(which('ft_warning'))
    % the fieldtrip/utilities directory contains the ft_hastoolbox and ft_warning
    % functions, which are required for the remainder of this script
    addpath(fullfile(fileparts(which('ft_defaults')), 'utilities'));
  end
  
  % Some people mess up their path settings and then have different versions of certain toolboxes on the path.
  % The following will issue a warning
  checkMultipleToolbox('FieldTrip',           'ft_defaults.m');
  checkMultipleToolbox('spm',                 'spm.m');
  checkMultipleToolbox('mne',                 'fiff_copy_tree.m');
  checkMultipleToolbox('eeglab',              'eeglab2fieldtrip.m');
  checkMultipleToolbox('dipoli',              'write_tri.m');
  checkMultipleToolbox('eeprobe',             'read_eep_avr.mexa64');
  checkMultipleToolbox('yokogawa',            'GetMeg160ChannelInfoM.p');
  checkMultipleToolbox('simbio',              'sb_compile_vista.m');
  checkMultipleToolbox('fns',                 'fns_region_read.m');
  checkMultipleToolbox('bemcp',               'bem_Cii_cst.mexa64');
  checkMultipleToolbox('bci2000',             'load_bcidat.m');
  checkMultipleToolbox('openmeeg',            'openmeeg_helper.m');
  checkMultipleToolbox('freesurfer',          'vox2ras_ksolve.m');
  checkMultipleToolbox('fastica',             'fastica.m');
  checkMultipleToolbox('besa',                'readBESAmul.m');
  checkMultipleToolbox('neuroshare',          'ns_GetAnalogData.m');
  checkMultipleToolbox('ctf',                 'setCTFDataBalance.m');
  checkMultipleToolbox('afni',                'WriteBrikHEAD.m');
  checkMultipleToolbox('gifti',               '@gifti/display.m');
  checkMultipleToolbox('sqdproject',          'sqdread.m');
  checkMultipleToolbox('xml4mat',             'xml2mat.m');
  checkMultipleToolbox('cca',                 'ccabss.m');
  checkMultipleToolbox('bsmart',              'armorf.m');
  checkMultipleToolbox('iso2mesh',            'iso2meshver.m');
  checkMultipleToolbox('bct',                 'degrees_und.m');
  checkMultipleToolbox('yokogawa_meg_reader', 'getYkgwHdrEvent.p');
  checkMultipleToolbox('biosig',              'sopen.m');
  checkMultipleToolbox('icasso',              'icassoEst.m');
  
  try
    % external/signal contains alternative implementations of some signal processing functions
    if ~ft_platform_supports('signal') || ~strcmp(ft_default.toolbox.signal, 'matlab') || ~ft_hastoolbox('signal')
      addpath(fullfile(fileparts(which('ft_defaults')), 'external', 'signal'));
    end
  end
  
  try
    % external/stats contains alternative implementations of some statistics functions
    if ~ft_platform_supports('stats') || ~strcmp(ft_default.toolbox.stats, 'matlab') || ~ft_hastoolbox('stats')
      addpath(fullfile(fileparts(which('ft_defaults')), 'external', 'stats'));
    end
  end
  
  try
    % external/images contains alternative implementations of some image processing functions
    if ~ft_platform_supports('images') || ~strcmp(ft_default.toolbox.images, 'matlab') || ~ft_hastoolbox('images')
      addpath(fullfile(fileparts(which('ft_defaults')), 'external', 'images'));
    end
  end
  
  % this directory contains various functions that were obtained from elsewere, e.g. MATLAB file exchange
  addtopath = {'fileexchange'};

  % these directories deal with compatibility with older MATLAB versions and with Octave
  if ft_platform_supports('matlabversion', -inf, '2008a'), addtopath{end+1} = 'compat/matlablt2008b'; end
  if ft_platform_supports('matlabversion', -inf, '2008b'), addtopath{end+1} = 'compat/matlablt2009a'; end
  if ft_platform_supports('matlabversion', -inf, '2009a'), addtopath{end+1} = 'compat/matlablt2009b'; end
  if ft_platform_supports('matlabversion', -inf, '2009b'), addtopath{end+1} = 'compat/matlablt2010a'; end
  if ft_platform_supports('matlabversion', -inf, '2010a'), addtopath{end+1} = 'compat/matlablt2010b'; end
  if ft_platform_supports('matlabversion', -inf, '2010b'), addtopath{end+1} = 'compat/matlablt2011a'; end
  if ft_platform_supports('matlabversion', -inf, '2011a'), addtopath{end+1} = 'compat/matlablt2011b'; end
  if ft_platform_supports('matlabversion', -inf, '2011b'), addtopath{end+1} = 'compat/matlablt2012a'; end
  if ft_platform_supports('matlabversion', -inf, '2012a'), addtopath{end+1} = 'compat/matlablt2012b'; end
  if ft_platform_supports('matlabversion', -inf, '2012b'), addtopath{end+1} = 'compat/matlablt2013a'; end
  if ft_platform_supports('matlabversion', -inf, '2013a'), addtopath{end+1} = 'compat/matlablt2013b'; end
  if ft_platform_supports('matlabversion', -inf, '2013b'), addtopath{end+1} = 'compat/matlablt2014a'; end
  if ft_platform_supports('matlabversion', -inf, '2014a'), addtopath{end+1} = 'compat/matlablt2014b'; end
  if ft_platform_supports('matlabversion', -inf, '2014d'), addtopath{end+1} = 'compat/matlablt2015a'; end
  if ft_platform_supports('matlabversion', -inf, '2015a'), addtopath{end+1} = 'compat/matlablt2015b'; end
  if ft_platform_supports('matlabversion', -inf, '2015b'), addtopath{end+1} = 'compat/matlablt2016a'; end
  if ft_platform_supports('matlabversion', -inf, '2016a'), addtopath{end+1} = 'compat/matlablt2016b'; end
  if ft_platform_supports('matlabversion', -inf, '2016b'), addtopath{end+1} = 'compat/matlablt2017a'; end
  if ft_platform_supports('matlabversion', -inf, '2017a'), addtopath{end+1} = 'compat/matlablt2017b'; end
  if ft_platform_supports('matlabversion', -inf, '2017b'), addtopath{end+1} = 'compat/matlablt2018a'; end
  if ft_platform_supports('matlabversion', -inf, '2018a'), addtopath{end+1} = 'compat/matlablt2018b'; end
  if ft_platform_supports('matlabversion', -inf, '2018b'), addtopath{end+1} = 'compat/matlablt2019a'; end
  if ft_platform_supports('matlabversion', -inf, '2019a'), addtopath{end+1} = 'compat/matlablt2019b'; end
  if ft_platform_supports('matlabversion', -inf, '2019b'), addtopath{end+1} = 'compat/matlablt2020a'; end
  if ft_platform_supports('matlabversion', -inf, '2020a'), addtopath{end+1} = 'compat/matlablt2020b'; end
  if ft_platform_supports('matlabversion', -inf, '2020b'), addtopath{end+1} = 'compat/matlablt2021a'; end
  if ft_platform_supports('matlabversion', -inf, '2021a'), addtopath{end+1} = 'compat/matlablt2021b'; end
  if ft_platform_supports('matlabversion', -inf, '2021b'), addtopath{end+1} = 'compat/matlablt2022a'; end
  if ft_platform_supports('matlabversion', -inf, '2022a'), addtopath{end+1} = 'compat/matlablt2022b'; end
  if ft_platform_supports('matlabversion', -inf, '2022b'), addtopath{end+1} = 'compat/matlablt2023a'; end
  if ft_platform_supports('matlabversion', -inf, '2023a'), addtopath{end+1} = 'compat/matlablt2023b'; end
  if ft_platform_supports('matlabversion', -inf, '2023b'), addtopath{end+1} = 'compat/matlablt2024a'; end
  if ft_platform_supports('matlabversion', -inf, '2024a'), addtopath{end+1} = 'compat/matlablt2024b'; end
  if ft_platform_supports('matlabversion', -inf, '2024b'), addtopath{end+1} = 'compat/matlablt2025a'; end
  % this deals with compatibility with all OCTAVE versions
  if ft_platform_supports('octaveversion', -inf, +inf),    addtopath{end+1} = 'compat/octave'; end

  try
    ft_hastoolbox(addtopath, 3, 1); % not required
  end
  
  addtopath        = {'template/layout' 'template/anatomy' 'template/headmodel' 'template/electrode' 'template/neighbours' 'template/sourcemodel'}; % these contains template layouts, neighbour structures, MRIs and cortical meshes
  addtopath{end+1} = 'statfun';        % this is used in ft_statistics
  addtopath{end+1} = 'trialfun';       % this is used in ft_definetrial
  addtopath{end+1} = 'fileio';         % this contains the low-level reading functions
  addtopath{end+1} = 'preproc';        % this is for filtering etc. on time-series data
  addtopath{end+1} = 'forward';        % this contains forward models for the EEG and MEG volume conductor
  addtopath{end+1} = 'inverse';        % this contains inverse source estimation methods
  addtopath{end+1} = 'plotting';       % this contains intermediate-level plotting functions, e.g. multiplots and 3-d objects
  addtopath{end+1} = 'specest';        % this contains intermediate-level functions for spectral analysis
  addtopath{end+1} = 'connectivity';   % this contains the functions to compute connectivity metrics
  addtopath{end+1} = 'test';           % this contains test scripts
  addtopath{end+1} = 'contrib/spike';  % this contains the functions for spike and spike-field analysis
  addtopath{end+1} = 'contrib/misc';   % this contains user contributed functions
  
  try 
    ft_hastoolbox(addtopath, 1, 1); % required
  end
   
  % try to add a few more, but it is not a show stopper if these fail
  try
    % this contains specific code and examples for realtime processing
    ft_hastoolbox({'realtime/example' 'realtime/online_mri' 'realtime/online_meg' 'realtime/online_eeg'}, 3, 1); % not required
  end
  
end % if not deployed

% the toolboxes added by this function should not be removed by FT_POSTAMBLE_HASTOOLBOX
ft_default.toolbox.cleanup = prevcleanup;

% track the usage of this function, this only happens once at startup
ft_trackusage('startup');

% remember that the function has executed in a persistent variable
initialized = true;

end % function ft_default


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkMultipleToolbox(toolbox, keyfile)

persistent warned
if isempty(warned)
  warned = false;
end

if ~ft_platform_supports('which-all')
  return;
end

list = which(keyfile, '-all');
if length(list)>1
  ft_warning('Multiple versions of %s on your path will confuse FieldTrip', toolbox);
  if ~warned % only throw the following warnings once
    warned = true;
    for i=1:length(list)
      ft_warning('one version of %s is found here: %s', toolbox, list{i});
    end
  end
  ft_warning('You probably used addpath(genpath(''path_to_fieldtrip'')), this can lead to unexpected behavior. See http://www.fieldtriptoolbox.org/faq/should_i_add_fieldtrip_with_all_subdirectories_to_my_matlab_path');
end
end % function checkMultipleToolbox

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkIncorrectPath
p = fileparts(mfilename('fullpath'));
incorrect = fullfile(p, 'compat', 'incorrect');
if contains(path, incorrect)
  ft_warning('Your path is set up incorrectly. You probably used addpath(genpath(''path_to_fieldtrip'')), this can lead to unexpected behavior. See http://www.fieldtriptoolbox.org/faq/should_i_add_fieldtrip_with_all_subdirectories_to_my_matlab_path');
end
end % function checkIncorrectPath
