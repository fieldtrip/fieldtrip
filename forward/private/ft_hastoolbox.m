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

% Copyright (C) 2005-2010, Robert Oostenveld
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

% this function is called many times in FieldTrip and associated toolboxes
% use efficient handling if the same toolbox has been investigated before
persistent previous previouspath

if ~isequal(previouspath, matlabpath)
  previous = [];
end

if isempty(previous)
  previous = struct;
elseif isfield(previous, fixname(toolbox))
  status = previous.(fixname(toolbox));
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
  'MEG-PD'     'see http://www.kolumbus.fi/kuutela/programs/meg-pd'
  'MEG-CALC'   'this is a commercial toolbox from Neuromag, see http://www.neuromag.com'
  'BIOSIG'     'see http://biosig.sourceforge.net'
  'EEG'        'see http://eeg.sourceforge.net'
  'EEGSF'      'see http://eeg.sourceforge.net'  % alternative name
  'MRI'        'see http://eeg.sourceforge.net'  % alternative name
  'NEUROSHARE' 'see http://www.neuroshare.org'
  'BESA'       'see http://www.megis.de, or contact Karsten Hoechstetter'
  'EEPROBE'    'see http://www.ant-neuro.com, or contact Maarten van der Velde'
  'YOKOGAWA'   'see http://www.yokogawa.co.jp, or contact Nobuhiko Takahashi'
  'BEOWULF'    'see http://oostenveld.net, or contact Robert Oostenveld'
  'MENTAT'     'see http://oostenveld.net, or contact Robert Oostenveld'
  'SON2'       'see http://www.kcl.ac.uk/depsta/biomedical/cfnr/lidierth.html, or contact Malcolm Lidierth'
  '4D-VERSION' 'contact Christian Wienbruch'
  'SIGNAL'     'see http://www.mathworks.com/products/signal'
  'OPTIM'      'see http://www.mathworks.com/products/optim'
  'IMAGE'      'see http://www.mathworks.com/products/image'
  'SPLINES'    'see http://www.mathworks.com/products/splines'
  'FASTICA'    'see http://www.cis.hut.fi/projects/ica/fastica'
  'BRAINSTORM' 'see http://neuroimage.ucs.edu/brainstorm'
  'FILEIO'     'see http://www.ru.nl/neuroimaging/fieldtrip'
  'FORWINV'    'see http://www.ru.nl/neuroimaging/fieldtrip'
  'PLOTTING'   'see http://www.ru.nl/neuroimaging/fieldtrip'
  'DENOISE'    'see http://lumiere.ens.fr/Audition/adc/meg, or contact Alain de Cheveigne'
  'BCI2000'    'see http://bci2000.org'
  'NLXNETCOM'  'see http://www.neuralynx.com'
  'DIPOLI'     'see ftp://ftp.fcdonders.nl/pub/fieldtrip/external'
  'MNE'        'see http://www.nmr.mgh.harvard.edu/martinos/userInfo/data/sofMNE.php'
  'TCP_UDP_IP' 'see http://www.mathworks.com/matlabcentral/fileexchange/345, or contact Peter Rydes?ter'
  'BEMCP'      'contact Christophe Phillips'
  'OPENMEEG'   'see http://gforge.inria.fr/projects/openmeeg and http://gforge.inria.fr/frs/?group_id=435'
  'PRTOOLS'    'see http://www.prtools.org'
  'ITAB'       'contact Stefania Della Penna'
  'BSMART'     'see http://www.brain-smart.org'
  'PEER'       'see http://fieldtrip.fcdonders.nl/development/peer'
  'FREESURFER' 'see http://surfer.nmr.mgh.harvard.edu/fswiki'
  'SIMBIO'     'see https://www.mrt.uni-jena.de/simbio/index.php/Main_Page'
  'FNS'        'see http://hhvn.nmsu.edu/wiki/index.php/FNS'
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
switch toolbox
  case 'AFNI'
    status = (exist('BrikLoad') && exist('BrikInfo'));
  case 'DSS'
    status = exist('denss', 'file') && exist('dss_create_state', 'file');
  case 'EEGLAB'
    status = exist('runica', 'file');
  case 'NWAY'
    status = exist('parafac', 'file');
  case 'SPM'
    status = exist('spm.m'); % any version of SPM is fine
  case 'SPM99'
    status = exist('spm.m') && strcmp(spm('ver'),'SPM99');
  case 'SPM2'
    status = exist('spm.m') && strcmp(spm('ver'),'SPM2');
  case 'SPM5'
    status = exist('spm.m') && strcmp(spm('ver'),'SPM5');
  case 'SPM8'
    status = exist('spm.m') && strncmp(spm('ver'),'SPM8', 4);
  case 'MEG-PD'
    status = (exist('rawdata') && exist('channames'));
  case 'MEG-CALC'
    status = (exist('megmodel') && exist('megfield') && exist('megtrans'));
  case 'BIOSIG'
    status = (exist('sopen') && exist('sread'));
  case 'EEG'
    status = (exist('ctf_read_res4') && exist('ctf_read_meg4'));
  case 'EEGSF'  % alternative name
    status = (exist('ctf_read_res4') && exist('ctf_read_meg4'));
  case 'MRI'    % other functions in the mri section
    status = (exist('avw_hdr_read') && exist('avw_img_read'));
  case 'NEUROSHARE'
    status  = (exist('ns_OpenFile') && exist('ns_SetLibrary') && exist('ns_GetAnalogData'));
  case 'BESA'
    status = (exist('readBESAtfc') && exist('readBESAswf'));
  case 'EEPROBE'
    status  = (exist('read_eep_avr') && exist('read_eep_cnt'));
  case 'YOKOGAWA'
      status = (exist('hasyokogawa') && strcmp(hasyokogawa, '16bitBeta6'));
  case 'YOKOGAWA16bitBeta3'
    status = (exist('hasyokogawa') && strcmp(hasyokogawa, '16bitBeta3'));
  case 'YOKOGAWA16bitBeta6'
    status = (exist('hasyokogawa') && strcmp(hasyokogawa, '16bitBeta6'));
  case 'BEOWULF'
    status = (exist('evalwulf') && exist('evalwulf') && exist('evalwulf'));
  case 'MENTAT'
    status  = (exist('pcompile') && exist('pfor') && exist('peval'));
  case 'SON2'
    status  = (exist('SONFileHeader') && exist('SONChanList') && exist('SONGetChannel'));
  case '4D-VERSION'
    status  = (exist('read4d') && exist('read4dhdr'));
  case {'STATS', 'STATISTICS'}
    status = license('checkout', 'statistics_toolbox');         % also check the availability of a toolbox license
  case {'OPTIM', 'OPTIMIZATION'}
    status = license('checkout', 'optimization_toolbox');       % also check the availability of a toolbox license
  case {'SPLINES', 'CURVE_FITTING'}
    status = license('checkout', 'curve_fitting_toolbox');      % also check the availability of a toolbox license
  case 'SIGNAL'
    status = license('checkout', 'signal_toolbox');             % also check the availability of a toolbox license
  case 'IMAGE'
    status = license('checkout', 'image_toolbox');              % also check the availability of a toolbox license
  case 'DCT'
    status = license('checkout', 'distrib_computing_toolbox');  % also check the availability of a toolbox license
  case 'FASTICA'
    status  = exist('fastica', 'file');
  case 'BRAINSTORM'
    status  = exist('bem_xfer');
  case 'FILEIO'
    status  = (exist('ft_read_header') && exist('ft_read_data') && exist('ft_read_event') && exist('ft_read_sens'));
  case 'FORMWARD'
    status  = (exist('ft_compute_leadfield') && exist('ft_prepare_vol_sens'));
  case 'DENOISE'
    status  = (exist('tsr') && exist('sns'));
  case 'CTF'
    status  = (exist('getCTFBalanceCoefs') && exist('getCTFdata'));
  case 'BCI2000'
    status  = exist('load_bcidat');
  case 'NLXNETCOM'
    status  = (exist('MatlabNetComClient') && exist('NlxConnectToServer') && exist('NlxGetNewCSCData'));
  case 'DIPOLI'
    status  = exist('dipoli.m', 'file');
  case 'MNE'
    status  = (exist('fiff_read_meas_info', 'file') && exist('fiff_setup_read_raw', 'file'));
  case 'TCP_UDP_IP'
    status  = (exist('pnet', 'file') && exist('pnet_getvar', 'file') && exist('pnet_putvar', 'file'));
  case 'BEMCP'
    status  = (exist('bem_Cij_cog', 'file') && exist('bem_Cij_lin', 'file') && exist('bem_Cij_cst', 'file'));
  case 'OPENMEEG'
    status = exist('openmeeg.m', 'file');
  case 'PLOTTING'
    status  = (exist('ft_plot_topo', 'file') && exist('ft_plot_mesh', 'file') && exist('ft_plot_matrix', 'file'));
  case 'PRTOOLS'
    status  = (exist('prversion', 'file') && exist('dataset', 'file') && exist('svc', 'file'));
  case 'ITAB'
    status  = (exist('lcReadHeader', 'file') && exist('lcReadData', 'file'));
  case 'BSMART' 
    status  = exist('bsmart'); 
  case 'PEER' 
    status  = exist('peerslave', 'file') && exist('peermaster', 'file');
  case 'CONNECTIVITY'
    status  = exist('ft_connectivity_corr', 'file') && exist('ft_connectivity_granger', 'file');
  case 'FREESURFER'
    status  = exist('MRIread', 'file') && exist('vox2ras_0to1', 'file');
  case 'FNS'
    status  = exist('elecsfwd', 'file') && exist('img_get_gray', 'file');
  case 'SIMBIO'
    status  = exist('ipm_linux_opt_Venant', 'file');
  otherwise
    if ~silent, warning('cannot determine whether the %s toolbox is present', toolbox); end
    status = 0;
end

% it should be a boolean value
status = (status~=0);

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

  % for linux computers in the F.C. Donders Centre
  prefix = '/home/common/matlab';
  if ~status && (strcmp(computer, 'GLNX86') || strcmp(computer, 'GLNXA64'))
    status = myaddpath(fullfile(prefix, lower(toolbox)), silent);
  end

  % for windows computers in the F.C. Donders Centre
  prefix = 'h:\common\matlab';
  if ~status && (strcmp(computer, 'PCWIN') || strcmp(computer, 'PCWIN64'))
    status = myaddpath(fullfile(prefix, lower(toolbox)), silent);
  end

  % use the matlab subdirectory in your homedirectory, this works on unix and mac
  prefix = [getenv('HOME') '/matlab'];
  if ~status
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
      warning(msg);
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
previouspath = matlabpath;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = myaddpath(toolbox, silent)
if exist(toolbox, 'dir')
  if ~silent, warning('adding %s toolbox to your Matlab path', toolbox); end
  addpath(toolbox);
  status = 1;
else
  status = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = fixname(toolbox)
out = lower(toolbox);
out(out=='-') = '_'; % fix dashes
out(out==' ') = '_'; % fix spaces
out(out=='/') = '_'; % fix forward slashes
out(out=='\') = '_'; % fix backward slashes

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

