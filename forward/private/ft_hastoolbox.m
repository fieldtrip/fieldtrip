function [status] = ft_hastoolbox(toolbox, autoadd, silent)

% FT_HASTOOLBOX tests whether an external toolbox is installed. Optionally
% it will try to determine the path to the toolbox and install it
% automatically.
%
% Use as
%   [status] = ft_hastoolbox(toolbox, autoadd, silent)
%
% autoadd = 0 means that it will not be added
% autoadd = 1 means that give an error if it cannot be added
% autoadd = 2 means that give a warning if it cannot be added
% autoadd = 3 means that it remains silent if it cannot be added
%
% silent = 0 means that it will give some feedback about adding the toolbox
% silent = 1 means that it will not give feedback

% Copyright (C) 2005-2013, Robert Oostenveld
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

% this function is called many times in FieldTrip and associated toolboxes
% use efficient handling if the same toolbox has been investigated before
% persistent previous previouspath
%
% if ~isequal(previouspath, path)
%   previous = [];
% end
%
% if isempty(previous)
%   previous = struct;
% elseif isfield(previous, fixname(toolbox))
%   status = previous.(fixname(toolbox));
%   return
% end

if isdeployed
  % it is not possible to check the presence of functions or change the path in a compiled application
  status = 1;
  return
end

% this points the user to the website where he/she can download the toolbox
url = {
  'AFNI'       'see http://afni.nimh.nih.gov'
  'DSS'        'see http://www.cis.hut.fi/projects/dss'
  'EEGLAB'     'see http://www.sccn.ucsd.edu/eeglab'
  'NWAY'       'see http://www.models.kvl.dk/source/nwaytoolbox'
  'SPM99'      'see http://www.fil.ion.ucl.ac.uk/spm'
  'SPM2'       'see http://www.fil.ion.ucl.ac.uk/spm'
  'SPM5'       'see http://www.fil.ion.ucl.ac.uk/spm'
  'SPM8'       'see http://www.fil.ion.ucl.ac.uk/spm'
  'SPM12'      'see http://www.fil.ion.ucl.ac.uk/spm'
  'MEG-PD'     'see http://www.kolumbus.fi/kuutela/programs/meg-pd'
  'MEG-CALC'   'this is a commercial toolbox from Neuromag, see http://www.neuromag.com'
  'BIOSIG'     'see http://biosig.sourceforge.net'
  'EEG'        'see http://eeg.sourceforge.net'
  'EEGSF'      'see http://eeg.sourceforge.net'  % alternative name
  'MRI'        'see http://eeg.sourceforge.net'  % alternative name
  'NEUROSHARE' 'see http://www.neuroshare.org'
  'BESA'         'see http://www.besa.de/downloads/matlab/ and get the "BESA MATLAB Readers"'
  'MATLAB2BESA'  'see http://www.besa.de/downloads/matlab/ and get the "MATLAB to BESA Export functions"'
  'EEPROBE'    'see http://www.ant-neuro.com, or contact Maarten van der Velde'
  'YOKOGAWA'   'this is deprecated, please use YOKOGAWA_MEG_READER instead'
  'YOKOGAWA_MEG_READER' 'see http://www.yokogawa.com/me/me-login-en.htm'
  'BEOWULF'    'see http://robertoostenveld.nl, or contact Robert Oostenveld'
  'MENTAT'     'see http://robertoostenveld.nl, or contact Robert Oostenveld'
  'SON2'       'see http://www.kcl.ac.uk/depsta/biomedical/cfnr/lidierth.html, or contact Malcolm Lidierth'
  '4D-VERSION' 'contact Christian Wienbruch'
  'COMM'       'see http://www.mathworks.com/products/communications'
  'SIGNAL'     'see http://www.mathworks.com/products/signal'
  'OPTIM'      'see http://www.mathworks.com/products/optim'
  'IMAGE'      'see http://www.mathworks.com/products/image'  % Mathworks refers to this as IMAGES
  'SPLINES'    'see http://www.mathworks.com/products/splines'
  'DISTCOMP'   'see http://www.mathworks.nl/products/parallel-computing/'
  'COMPILER'   'see http://www.mathworks.com/products/compiler'
  'FASTICA'    'see http://www.cis.hut.fi/projects/ica/fastica'
  'BRAINSTORM' 'see http://neuroimage.ucs.edu/brainstorm'
  'FILEIO'     'see http://www.fieldtriptoolbox.org'
  'PREPROC'    'see http://www.fieldtriptoolbox.org'
  'FORWARD'    'see http://www.fieldtriptoolbox.org'
  'INVERSE'    'see http://www.fieldtriptoolbox.org'
  'SPECEST'    'see http://www.fieldtriptoolbox.org'
  'REALTIME'   'see http://www.fieldtriptoolbox.org'
  'PLOTTING'   'see http://www.fieldtriptoolbox.org'
  'SPIKE'      'see http://www.fieldtriptoolbox.org'
  'CONNECTIVITY' 'see http://www.fieldtriptoolbox.org'
  'PEER'          'see http://www.fieldtriptoolbox.org'
  'PLOTTING'      'see http://www.fieldtriptoolbox.org'
  'DENOISE'       'see http://lumiere.ens.fr/Audition/adc/meg, or contact Alain de Cheveigne'
  'BCI2000'       'see http://bci2000.org'
  'NLXNETCOM'     'see http://www.neuralynx.com'

  'DIPOLI'        'see ftp://ftp.fcdonders.nl/pub/fieldtrip/external'
  'MNE'           'see http://www.nmr.mgh.harvard.edu/martinos/userInfo/data/sofMNE.php'
  'TCP_UDP_IP'    'see http://www.mathworks.com/matlabcentral/fileexchange/345, or contact Peter Rydesaeter'
  'BEMCP'         'contact Christophe Phillips'
  'OPENMEEG'      'see http://gforge.inria.fr/projects/openmeeg and http://gforge.inria.fr/frs/?group_id=435'
  'PRTOOLS'       'see http://www.prtools.org'
  'ITAB'          'contact Stefania Della Penna'
  'BSMART'        'see http://www.brain-smart.org'
  'PEER'          'see http://www.fieldtriptoolbox.org/development/peer'
  'FREESURFER'    'see http://surfer.nmr.mgh.harvard.edu/fswiki'
  'SIMBIO'        'see https://www.mrt.uni-jena.de/simbio/index.php/Main_Page'
  'VGRID'         'see http://www.rheinahrcampus.de/~medsim/vgrid/manual.html'
  'FNS'           'see http://hhvn.nmsu.edu/wiki/index.php/FNS'
  'GIFTI'         'see http://www.artefact.tk/software/matlab/gifti'
  'XML4MAT'       'see http://www.mathworks.com/matlabcentral/fileexchange/6268-xml4mat-v2-0'
  'SQDPROJECT'    'see http://www.isr.umd.edu/Labs/CSSL/simonlab'
  'BCT'           'see http://www.brain-connectivity-toolbox.net/'
  'CCA'           'see http://www.imt.liu.se/~magnus/cca or contact Magnus Borga'
  'EGI_MFF'       'see http://www.egi.com/ or contact either Phan Luu or Colin Davey at EGI'
  'TOOLBOX_GRAPH' 'see http://www.mathworks.com/matlabcentral/fileexchange/5355-toolbox-graph or contact Gabriel Peyre'
  'NETCDF'        'see http://www.mathworks.com/matlabcentral/fileexchange/15177'
  'MYSQL'         'see http://www.mathworks.com/matlabcentral/fileexchange/8663-mysql-database-connector'
  'ISO2MESH'      'see http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Home or contact Qianqian Fang'
  'DATAHASH'      'see http://www.mathworks.com/matlabcentral/fileexchange/31272'
  'IBTB'          'see http://www.ibtb.org'
  'ICASSO'        'see http://www.cis.hut.fi/projects/ica/icasso'
  'XUNIT'         'see http://www.mathworks.com/matlabcentral/fileexchange/22846-matlab-xunit-test-framework'
  'PLEXON'        'available from http://www.plexon.com/assets/downloads/sdk/ReadingPLXandDDTfilesinMatlab-mexw.zip'
  'MISC'          'various functions that were downloaded from http://www.mathworks.com/matlabcentral/fileexchange and elsewhere'
  '35625-INFORMATION-THEORY-TOOLBOX'      'see http://www.mathworks.com/matlabcentral/fileexchange/35625-information-theory-toolbox'
  '29046-MUTUAL-INFORMATION'              'see http://www.mathworks.com/matlabcentral/fileexchange/35625-information-theory-toolbox'
  '14888-MUTUAL-INFORMATION-COMPUTATION'  'see http://www.mathworks.com/matlabcentral/fileexchange/14888-mutual-information-computation'
  'PLOT2SVG'      'see http://www.mathworks.com/matlabcentral/fileexchange/7401-scalable-vector-graphics-svg-export-of-figures'
  'BRAINSUITE'    'see http://brainsuite.bmap.ucla.edu/processing/additional-tools/'
  'BRAINVISA'     'see http://brainvisa.info'
  'FILEEXCHANGE'  'see http://www.mathworks.com/matlabcentral/fileexchange/'
  'NEURALYNX_V6'  'see http://neuralynx.com/research_software/file_converters_and_utilities/ and take the version from Neuralynx (windows only)'
  'NEURALYNX_V3'  'see http://neuralynx.com/research_software/file_converters_and_utilities/ and take the version from Ueli Rutishauser'
  'NPMK'          'see https://github.com/BlackrockMicrosystems/NPMK'
  'VIDEOMEG'      'see https://github.com/andreyzhd/VideoMEG'
  'WAVEFRONT'     'see http://mathworks.com/matlabcentral/fileexchange/27982-wavefront-obj-toolbox'
  'NEURONE'       'see http://www.megaemg.com/support/unrestricted-downloads'
  };

if nargin<2
  % default is not to add the path automatically
  autoadd = 0;
end

if nargin<3
  % default is not to be silent
  silent = 0;
end

% determine whether the toolbox is installed
toolbox = upper(toolbox);

% In case SPM8 or higher not available, allow to use fallback toolbox
fallback_toolbox='';
switch toolbox
  case 'AFNI'
    dependency={'BrikLoad', 'BrikInfo'};
  case 'DSS'
    dependency={'denss', 'dss_create_state'};
  case 'EEGLAB'
    dependency = 'runica';
  case 'NWAY'
    dependency = 'parafac';
  case 'SPM'
    dependency = 'spm'; % any version of SPM is fine
  case 'SPM99'
    dependency = {'spm', get_spm_version()==99};
  case 'SPM2'
    dependency = {'spm', get_spm_version()==2};
  case 'SPM5'
    dependency = {'spm', get_spm_version()==5};
  case 'SPM8'
    dependency = {'spm', get_spm_version()==8};
  case 'SPM8UP' % version 8 or later, but not SPM 9X
    dependency = {'spm', get_spm_version()>=8, get_spm_version()<95};

    %This is to avoid crashes when trying to add SPM to the path
    fallback_toolbox = 'SPM8';
  case 'SPM12'
    dependency = {'spm', get_spm_version()==12};
  case 'MEG-PD'
    dependency = {'rawdata', 'channames'};
  case 'MEG-CALC'
    dependency = {'megmodel', 'megfield', 'megtrans'};
  case 'BIOSIG'
    dependency = {'sopen', 'sread'};
  case 'EEG'
    dependency = {'ctf_read_res4', 'ctf_read_meg4'};
  case 'EEGSF'  % alternative name
    dependency = {'ctf_read_res4', 'ctf_read_meg4'};
  case 'MRI'    % other functions in the mri section
    dependency = {'avw_hdr_read', 'avw_img_read'};
  case 'NEUROSHARE'
    dependency = {'ns_OpenFile', 'ns_SetLibrary', ...
                            'ns_GetAnalogData'};
  case 'ARTINIS'
    dependency = {'read_artinis_oxy3'};
  case 'BESA'
    dependency = {'readBESAavr', 'readBESAelp', 'readBESAswf'};
  case 'MATLAB2BESA'
    dependency = {'besa_save2Avr', 'besa_save2Elp', 'besa_save2Swf'};
  case 'EEPROBE'
    dependency = {'read_eep_avr', 'read_eep_cnt'};
  case 'YOKOGAWA'
    dependency = @()hasyokogawa('16bitBeta6');
  case 'YOKOGAWA12BITBETA3'
    dependency = @()hasyokogawa('12bitBeta3');
  case 'YOKOGAWA16BITBETA3'
    dependency = @()hasyokogawa('16bitBeta3');
  case 'YOKOGAWA16BITBETA6'
    dependency = @()hasyokogawa('16bitBeta6');
  case 'YOKOGAWA_MEG_READER'
    dependency = @()hasyokogawa('1.4');
  case 'BEOWULF'
    dependency = {'evalwulf', 'evalwulf', 'evalwulf'};
  case 'MENTAT'
    dependency = {'pcompile', 'pfor', 'peval'};
  case 'SON2'
    dependency = {'SONFileHeader', 'SONChanList', 'SONGetChannel'};
  case '4D-VERSION'
    dependency  = {'read4d', 'read4dhdr'};
  case {'STATS', 'STATISTICS'}
    dependency = has_license('statistics_toolbox');         % check the availability of a toolbox license
  case {'OPTIM', 'OPTIMIZATION'}
    dependency = has_license('optimization_toolbox');       % check the availability of a toolbox license
  case {'SPLINES', 'CURVE_FITTING'}
    dependency = has_license('curve_fitting_toolbox');      % check the availability of a toolbox license
  case 'COMM'
    dependency = {has_license('communication_toolbox'), 'de2bi'}; % also check the availability of a toolbox license
  case 'SIGNAL'
    dependency = {has_license('signal_toolbox'), 'window'}; % also check the availability of a toolbox license
  case 'IMAGE'
    dependency = has_license('image_toolbox');              % check the availability of a toolbox license
  case {'DCT', 'DISTCOMP'}
    dependency = has_license('distrib_computing_toolbox');  % check the availability of a toolbox license
  case 'COMPILER'
    dependency = has_license('compiler');                   % check the availability of a toolbox license
  case 'FASTICA'
    dependency = 'fpica';
  case 'BRAINSTORM'
    dependency = 'bem_xfer';
  case 'DENOISE'
    dependency = {'tsr', 'sns'};
  case 'CTF'
    dependency = {'getCTFBalanceCoefs', 'getCTFdata'};
  case 'BCI2000'
    dependency  = {'load_bcidat'};
  case 'NLXNETCOM'
    dependency = {'MatlabNetComClient', 'NlxConnectToServer', ...
                    'NlxGetNewCSCData'};
  case 'DIPOLI'
    dependency = {'dipoli.maci', 'file'};
  case 'MNE'
    dependency = {'fiff_read_meas_info', 'fiff_setup_read_raw'};
  case 'TCP_UDP_IP'
    dependency = {'pnet', 'pnet_getvar', 'pnet_putvar'};
  case 'BEMCP'
    dependency = {'bem_Cij_cog', 'bem_Cij_lin', 'bem_Cij_cst'};
  case 'OPENMEEG'
    dependency = {'om_save_tri'};
  case 'PRTOOLS'
    dependency = {'prversion', 'dataset', 'svc'};
  case 'ITAB'
    dependency = {'lcReadHeader', 'lcReadData'};
  case 'BSMART'
    dependency = 'bsmart';
  case 'FREESURFER'
    dependency = {'MRIread', 'vox2ras_0to1'};
  case 'FNS'
    dependency = 'elecsfwd';
  case 'SIMBIO'
    dependency = {'calc_stiff_matrix_val', 'sb_transfer'};
  case 'VGRID'
    dependency = 'vgrid';
  case 'GIFTI'
    dependency = 'gifti';
  case 'XML4MAT'
    dependency = {'xml2struct', 'xml2whos'};
  case 'SQDPROJECT'
    dependency = {'sqdread', 'sqdwrite'};
  case 'BCT'
    dependency = {'macaque71.mat', 'motif4funct_wei'};
  case 'CCA'
    dependency = {'ccabss'};
  case 'EGI_MFF'
    dependency = {'mff_getObject', 'mff_getSummaryInfo'};
  case 'TOOLBOX_GRAPH'
    dependency = 'toolbox_graph';
  case 'NETCDF'
    dependency = {'netcdf'};
  case 'MYSQL'
    % not sure if 'which' would work fine here, so use 'exist'
    dependency = has_mex('mysql'); % this only consists of a single mex file
  case 'ISO2MESH'
    dependency = {'vol2surf', 'qmeshcut'};
  case 'QSUB'
    dependency = {'qsubfeval', 'qsubcellfun'};
  case 'ENGINE'
    dependency = {'enginefeval', 'enginecellfun'};
  case 'DATAHASH'
    dependency = {'DataHash'};
  case 'IBTB'
    dependency = {'make_ibtb','binr'};
  case 'ICASSO'
    dependency = {'icassoEst'};
  case 'XUNIT'
    dependency = {'initTestSuite', 'runtests'};
  case 'PLEXON'
    dependency = {'plx_adchan_gains', 'mexPlex'};
  case '35625-INFORMATION-THEORY-TOOLBOX'
    dependency = {'conditionalEntropy', 'entropy', 'jointEntropy',...
                    'mutualInformation' 'nmi' 'nvi' 'relativeEntropy'};
  case '29046-MUTUAL-INFORMATION'
    dependency = {'MI', 'license.txt'};
  case '14888-MUTUAL-INFORMATION-COMPUTATION'
    dependency = {'condentropy', 'demo_mi', 'estcondentropy.cpp',...
                    'estjointentropy.cpp', 'estpa.cpp', ...
                    'findjointstateab.cpp', 'makeosmex.m',...
                    'mutualinfo.m', 'condmutualinfo.m',...
                    'entropy.m', 'estentropy.cpp',...
                    'estmutualinfo.cpp', 'estpab.cpp',...
                    'jointentropy.m' 'mergemultivariables.m' };
  case 'PLOT2SVG'
    dependency = {'plot2svg.m', 'simulink2svg.m'};
  case 'BRAINSUITE'
    dependency = {'readdfs.m', 'writedfc.m'};
  case 'BRAINVISA'
    dependency = {'loadmesh.m', 'plotmesh.m', 'savemesh.m'};
  case 'NEURALYNX_V6'
    dependency = has_mex('Nlx2MatCSC');
  case 'NEURALYNX_V3'
    dependency = has_mex('Nlx2MatCSC_v3');
  case 'NPMK'
    dependency = {'OpenNSx' 'OpenNEV'};
  case 'VIDEOMEG'
    dependency = {'comp_tstamps' 'load_audio0123', 'load_video123'};
  case 'WAVEFRONT'
    dependency = {'write_wobj' 'read_wobj'};
  case 'NEURONE'
    dependency = {'readneurone' 'readneuronedata' 'readneuroneevents'};

    % the following are fieldtrip modules/toolboxes
  case 'FILEIO'
    dependency = {'ft_read_header', 'ft_read_data', ...
                    'ft_read_event', 'ft_read_sens'};
  case 'FORWARD'
    dependency = {'ft_compute_leadfield', 'ft_prepare_vol_sens'};
  case 'PLOTTING'
    dependency = {'ft_plot_topo', 'ft_plot_mesh', 'ft_plot_matrix'};
  case 'PEER'
    dependency = {'peerslave', 'peermaster'};
  case 'CONNECTIVITY'
    dependency = {'ft_connectivity_corr', 'ft_connectivity_granger'};
  case 'SPIKE'
    dependency = {'ft_spiketriggeredaverage', 'ft_spiketriggeredspectrum'};
  case 'FILEEXCHANGE'
    dependency = is_subdir_in_fieldtrip_path('/external/fileexchange');
  case {'INVERSE', 'REALTIME', 'SPECEST', 'PREPROC', ...
          'COMPAT', 'STATFUN', 'TRIALFUN', 'UTILITIES/COMPAT', ...
          'FILEIO/COMPAT', 'PREPROC/COMPAT', 'FORWARD/COMPAT', ...
          'PLOTTING/COMPAT', 'TEMPLATE/LAYOUT', 'TEMPLATE/ANATOMY' ,...
          'TEMPLATE/HEADMODEL', 'TEMPLATE/ELECTRODE', ...
          'TEMPLATE/NEIGHBOURS', 'TEMPLATE/SOURCEMODEL'}
    dependency = is_subdir_in_fieldtrip_path(toolbox);
  otherwise
    if ~silent, warning('cannot determine whether the %s toolbox is present', toolbox); end
    dependency = false;
end

status = is_present(dependency);
if ~status && ~isempty(fallback_toolbox)
  % in case of SPM8UP
  toolbox = fallback_toolbox;
end

% try to determine the path of the requested toolbox
if autoadd>0 && ~status

  % for core fieldtrip modules
  prefix = fileparts(which('ft_defaults'));
  if ~status
    status = myaddpath(fullfile(prefix, lower(toolbox)), silent);
  end

  % for external fieldtrip modules
  prefix = fullfile(fileparts(which('ft_defaults')), 'external');
  if ~status
    status = myaddpath(fullfile(prefix, lower(toolbox)), silent);
    licensefile = [lower(toolbox) '_license'];
    if status && exist(licensefile, 'file')
      % this will execute openmeeg_license and mne_license
      % which display the license on screen for three seconds
      feval(licensefile);
    end
  end

  % for contributed fieldtrip extensions
  prefix = fullfile(fileparts(which('ft_defaults')), 'contrib');
  if ~status
    status = myaddpath(fullfile(prefix, lower(toolbox)), silent);
    licensefile = [lower(toolbox) '_license'];
    if status && exist(licensefile, 'file')
      % this will execute openmeeg_license and mne_license
      % which display the license on screen for three seconds
      feval(licensefile);
    end
  end

  % for linux computers in the Donders Centre for Cognitive Neuroimaging
  prefix = '/home/common/matlab';
  if ~status && isdir(prefix)
    status = myaddpath(fullfile(prefix, lower(toolbox)), silent);
  end

  % for windows computers in the Donders Centre for Cognitive Neuroimaging
  prefix = 'h:\common\matlab';
  if ~status && isdir(prefix)
    status = myaddpath(fullfile(prefix, lower(toolbox)), silent);
  end

  % use the MATLAB subdirectory in your homedirectory, this works on linux and mac
  prefix = fullfile(getenv('HOME'), 'matlab');
  if ~status && isdir(prefix)
    status = myaddpath(fullfile(prefix, lower(toolbox)), silent);
  end

  if ~status
    % the toolbox is not on the path and cannot be added
    sel = find(strcmp(url(:,1), toolbox));
    if ~isempty(sel)
      msg = sprintf('the %s toolbox is not installed, %s', toolbox, url{sel, 2});
    else
      msg = sprintf('the %s toolbox is not installed', toolbox);
    end
    if autoadd==1
      error(msg);
    elseif autoadd==2
      ft_warning(msg);
    else
      % fail silently
    end
  end
end

% this function is called many times in FieldTrip and associated toolboxes
% use efficient handling if the same toolbox has been investigated before
if status
  previous.(fixname(toolbox)) = status;
end

% remember the previous path, allows us to determine on the next call
% whether the path has been modified outise of this function
previouspath = path;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = myaddpath(toolbox, silent)
if isdeployed
  warning('cannot change path settings for %s in a compiled application', toolbox);
  status = 1;
elseif exist(toolbox, 'dir')
  if ~silent,
    ws = warning('backtrace', 'off');
    warning('adding %s toolbox to your MATLAB path', toolbox);
    warning(ws); % return to the previous warning level
  end
  addpath(toolbox);
  status = 1;
elseif (~isempty(regexp(toolbox, 'spm5$', 'once')) || ~isempty(regexp(toolbox, 'spm8$', 'once')) || ~isempty(regexp(toolbox, 'spm12$', 'once'))) && exist([toolbox 'b'], 'dir')
  status = myaddpath([toolbox 'b'], silent);
else
  status = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function path = unixpath(path)
%path(path=='\') = '/'; % replace backward slashes with forward slashes
path = strrep(path,'\','/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = hasfunction(funname, toolbox)
try
  % call the function without any input arguments, which probably is inapropriate
  feval(funname);
  % it might be that the function without any input already works fine
  status = true;
catch
  % either the function returned an error, or the function is not available
  % availability is influenced by the function being present and by having a
  % license for the function, i.e. in a concurrent licensing setting it might
  % be that all toolbox licenses are in use
  m = lasterror;
  if strcmp(m.identifier, 'MATLAB:license:checkouterror')
    if nargin>1
      warning('the %s toolbox is available, but you don''t have a license for it', toolbox);
    else
      warning('the function ''%s'' is available, but you don''t have a license for it', funname);
    end
    status = false;
  elseif strcmp(m.identifier, 'MATLAB:UndefinedFunction')
    status = false;
  else
    % the function seems to be available and it gave an unknown error,
    % which is to be expected with inappropriate input arguments
    status = true;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = is_subdir_in_fieldtrip_path(toolbox_name)
  fttrunkpath = unixpath(fileparts(which('ft_defaults')));
  fttoolboxpath = fullfile(fttrunkpath, lower(toolbox_name));

  needle=[pathsep fttoolboxpath pathsep];
  haystack = [pathsep path() pathsep];

  status = ~isempty(findstr(needle, haystack));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = has_mex(name)
  full_name=[name '.' mexext];
  status = (exist(full_name, 'file')==3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = get_spm_version()
  if ~is_present('spm')
    v=NaN;
    return
  end

  version_str = spm('ver');
  token = regexp(version_str,'(\d*)','tokens');
  v = str2num([token{:}{:}]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = has_license(toolbox_name)
  status = license('checkout', toolbox_name)==1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = is_present(dependency)
  if iscell(dependency)
    % use recursion
    status = all(cellfun(@is_present,dependency));
  elseif islogical(dependency)
    % boolean
    status = all(dependency);
  elseif ischar(dependency)
    % name of a function
    status = is_function_present_in_search_path(dependency);
  elseif isa(dependency, 'function_handle')
    status = dependency();
  else
    assert(false,'this should not happen');
  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = is_function_present_in_search_path(function_name)
  w = which(function_name);

  % must be in path and not a variable
  status = ~isempty(w) && ~isequal(w, 'variable');
