function ft_defaults

% FT_DEFAULTS contains some general settings and takes care of the required path
% settings. This function is called at the begin of all FieldTrip functions.
%
% The configuration defaults are stored in the global "ft_default" structure.
% The ft_checkconfig function that is called by many FieldTrip functions will
% merge this global ft_default structure with the cfg ctructure that you pass to
% the FieldTrip function that you are calling.
%
% The global options and their default values are
%   ft_default.trackconfig    string, can be cleanup, report, off (default = 'off')
%   ft_default.checkconfig    string, can be pedantic, loose, silent (default = 'loose')
%   ft_default.checksize      number in bytes, can be inf (default = 1e5)
%   ft_default.showcallinfo   string, can be yes or no (default = 'yes')
%
% See also FT_HASTOOLBOX, FT_CHECKCONFIG

% Note that this should be a function and not a script, otherwise the
% ft_hastoolbox function appears not be found in fieldtrip/private.

% Copyright (C) 2009-2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

% set the defaults in a global variable, ft_checkconfig will copy these over into the local configuration
% note that ft_getopt might not be available on the path at this moment
global ft_default
if ~isfield(ft_default, 'trackconfig'),    ft_default.trackconfig    = 'off';    end % cleanup, report, off
if ~isfield(ft_default, 'checkconfig'),    ft_default.checkconfig    = 'loose';  end % pedantic, loose, silent
if ~isfield(ft_default, 'checksize'),      ft_default.checksize      = 1e5;      end % number in bytes, can be inf
if ~isfield(ft_default, 'showcallinfo'),   ft_default.showcallinfo   = 'yes';    end % yes or no, this is used in ft_pre/postamble_provenance
if ~isfield(ft_default, 'debug'),          ft_default.debug          = 'no';     end % yes or no, this is used in ft_pre/postamble_debug

% these options allow to disable parts of the provenance
if ~isfield(ft_default, 'trackcallinfo'),  ft_default.trackcallinfo  = 'yes';    end % yes or no
if ~isfield(ft_default, 'trackdatainfo'),  ft_default.trackdatainfo  = 'no';     end % yes or no, this is still under development
if ~isfield(ft_default, 'trackparaminfo'), ft_default.trackparaminfo = 'no';     end % yes or no, this is still under development

% Ensure that the path containing ft_defaults is on the path.
% This allows people to do "cd path_to_fieldtrip; ft_defaults"
ftPath = fileparts(mfilename('fullpath')); % get path, strip away 'ft_defaults'
ftPath = strrep(ftPath, '\', '\\');
if isempty(regexp(path, [ftPath pathsep '|' ftPath '$'], 'once'))
  warning('FieldTrip is not yet on your MATLAB path, adding %s', ftPath);
  addpath(ftPath);
end

% Some people mess up their path settings and then have
% different versions of certain toolboxes on the path.
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

if ~isdeployed
  
  if isempty(which('ft_hastoolbox'))
    % the fieldtrip/utilities directory contains the ft_hastoolbox function
    % which is required for the remainder of this script
    addpath(fullfile(fileparts(which('ft_defaults')), 'utilities'));
  end
  
  try
    % this directory contains the backward compatibility wrappers for the ft_xxx function name change
    ft_hastoolbox('compat', 3, 1); % not required
  end
  
  try
    % this directory contains the backward compatibility wrappers for the fieldtrip/utilities functions
    ft_hastoolbox('utilities/compat', 3, 1);
  end
  
  try
    % these contains template layouts, neighbour structures, MRIs and cortical meshes
    ft_hastoolbox('template/layout', 1, 1);
    ft_hastoolbox('template/anatomy', 1, 1);
    ft_hastoolbox('template/headmodel', 1, 1);
    ft_hastoolbox('template/electrode', 1, 1);
    ft_hastoolbox('template/neighbours', 1, 1);
    ft_hastoolbox('template/sourcemodel', 1, 1);
  end
  
  try
    % this is used in statistics
    ft_hastoolbox('statfun', 1, 1);
  end
  
  try
    % this is used in definetrial
    ft_hastoolbox('trialfun', 1, 1);
  end
  
  try
    % this contains the low-level reading functions
    ft_hastoolbox('fileio', 1, 1);
    ft_hastoolbox('fileio/compat', 3, 1); % not required
  end
  
  try
    % this is for filtering time-series data
    ft_hastoolbox('preproc', 1, 1);
    ft_hastoolbox('preproc/compat', 3, 1); % not required
  end
  
  try
    % this contains forward models for the EEG and MEG volume conduction problem
    ft_hastoolbox('forward', 1, 1);
    ft_hastoolbox('forward/compat', 3, 1); % not required
  end
  
  try
    % numerous functions depend on this module
    ft_hastoolbox('inverse', 1, 1);
  end
  
  try
    % this contains intermediate-level plotting functions, e.g. multiplots and 3-d objects
    ft_hastoolbox('plotting', 1, 1);
    ft_hastoolbox('plotting/compat', 1, 1);
  end
  
  try
    % this contains the functions to compute connecitivy metrics
    ft_hastoolbox('connectivity', 1,1);
  end
  
  try
    % this contains the functions for spike and spike-field analysis
    ft_hastoolbox('spike', 1,1);
  end
  
  try
    % this contains specific code and examples for realtime processing
    ft_hastoolbox('realtime', 3, 1);                    % not required
    ft_hastoolbox('realtime/acquisition/matlab', 3, 1); % not required
  end
  
  try
    % this contains intermediate-level functions for spectral analysis
    ft_hastoolbox('specest', 1, 1);
  end
  
end

end % function ft_default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkMultipleToolbox(toolbox, keyfile)
list = which(keyfile, '-all');
if length(list)>1
  [ws, warned] = warning_once(sprintf('Multiple versions of %s on your path will confuse FieldTrip', toolbox));
  if warned % only throw the warning once
    for i=1:length(list)
      warning('one version of %s is found here: %s', toolbox, list{i});
    end
  end
  warning_once('You probably used addpath(genpath(''path_to_fieldtrip'')), this can lead to unexpected behaviour. See http://fieldtrip.fcdonders.nl/faq/should_i_add_fieldtrip_with_all_subdirectories_to_my_matlab_path');
end
end % function checkMultipleToolbox
