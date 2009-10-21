function [status] = hastoolbox(toolbox, autoadd, silent)

% HASTOOLBOX tests whether an external toolbox is installed. Optionally
% it will try to determine the path to the toolbox and install it
% automatically.
%
% Use as
%   [status] = hastoolbox(toolbox, autoadd, silent)

% Copyright (C) 2005-2008, Robert Oostenveld
%
% $Log: hastoolbox.m,v $
% Revision 1.37  2009/10/13 10:11:06  roboos
% added lc-libs
%
% Revision 1.36  2009/09/08 14:34:01  roboos
% also detect 64 bit windows version (thanks to arno)
%
% Revision 1.35  2009/04/21 09:54:15  roboos
% added prtools
%
% Revision 1.34  2009/04/02 19:47:33  roboos
% added plotting module
%
% Revision 1.33  2009/03/30 15:06:14  roboos
% added the patch from Alexandre to support openmeeg
%
% Revision 1.32  2009/03/12 10:40:21  roboos
% added the splines toolbox (mainly for testing) and changed the warning message related to the license
%
% Revision 1.31  2009/03/12 10:33:35  roboos
% not only check that a function is available, also check whether a license for that function is available
%
% Revision 1.30  2009/03/11 21:26:27  roboos
% detect spm8b just as spm8
%
% Revision 1.29  2009/03/11 10:35:19  roboos
% spm detection was confused with function and directory, explicitely check for "spm.m" which is the function
%
% Revision 1.28  2009/03/11 08:49:04  roboos
% improved the detection of the various spm versions
%
% Revision 1.27  2009/02/11 11:03:08  roboos
% changed naming of the functions of Chris in accordance with SPM8
%
% Revision 1.26  2009/02/02 12:57:21  roboos
% added bemcp, image, tcp_udp_ip
%
% Revision 1.25  2009/01/19 15:02:21  roboos
% added mne for fiff access
%
% Revision 1.24  2009/01/08 17:00:02  roboos
% improved caching in case the toolbox is not present
%
% Revision 1.23  2008/12/24 09:10:46  roboos
% added dipoli
%
% Revision 1.22  2008/10/29 15:45:12  roboos
% fix dashes and spaces in directory names for caching
%
% Revision 1.21  2008/10/20 21:50:56  roboos
% added NlxNetCom
%
% Revision 1.20  2008/10/20 16:31:15  roboos
% fixed problem in case with dash "-" in the directory
%
% Revision 1.19  2008/09/29 09:00:19  roboos
% implemented smart handling of previously seen toolboxes using a persistent variable
% this should speed up fieldtrip and fileio (e.g. read_data checks the presence of ctf for every trial)
%
% Revision 1.18  2008/09/24 15:43:00  roboos
% added read_data and read_sens for fileio, should solve problem for MEEGfileio in spm5
%
% Revision 1.17  2008/09/22 19:42:09  roboos
% added option for silent processing
%
% Revision 1.16  2008/08/11 16:11:19  roboos
% also automatically add to path for fieldtrip code and external modules
%
% Revision 1.15  2008/06/20 07:25:56  roboos
% added check for presence of BCI2000 load_bcidat mex file
%
% Revision 1.14  2008/05/15 10:52:29  roboos
% added ctf
%
% Revision 1.13  2008/03/17 08:29:40  roboos
% changed some contact addresses
%
% Revision 1.12  2008/03/14 10:20:29  roboos
% added denoise
%
% Revision 1.11  2008/03/05 10:59:14  roboos
% added fileio and forwinv
%
% Revision 1.10  2007/05/06 09:10:07  roboos
% added spm5
%
% Revision 1.9  2007/02/26 13:41:07  roboos
% made small change to fastica detection (suggested by Sameer)
%
% Revision 1.8  2007/02/13 17:22:27  roboos
% added MRI from eeg.sf.net
%
% Revision 1.7  2007/02/13 14:01:26  roboos
% added brainstorm
%
% Revision 1.6  2007/02/12 19:43:23  roboos
% added fastica, optim
%
% Revision 1.5  2007/01/17 17:05:34  roboos
% added matlab signal processing toolbox
%
% Revision 1.4  2007/01/04 12:25:19  roboos
% added SON2
%
% Revision 1.3  2007/01/03 17:01:15  roboos
% added 4d-version toolbox
%
% Revision 1.2  2006/06/07 10:48:02  roboos
% changed the "see xxx" string
%
% Revision 1.1  2006/06/07 09:28:41  roboos
% renamed fieldtrip/private/checktoolbox into misc/hastoolbox
%
% Revision 1.8  2006/06/06 14:18:22  roboos
% added neuroshare, eeprobe, yokogawa
%
% Revision 1.7  2006/05/31 08:56:24  roboos
% implemented detection of toolbox in users ~/matlab/toolboxname
%
% Revision 1.6  2006/05/23 10:20:46  roboos
% added beowulf and mentat toolboxes
%
% Revision 1.5  2006/04/26 11:37:22  roboos
% added besa toolbox
%
% Revision 1.4  2006/02/07 20:01:39  roboos
% aded biosig and meg-pd (neuromag)
%
% Revision 1.3  2006/01/17 14:05:54  roboos
% added GLNA64 for mentat000
%
% Revision 1.2  2006/01/06 11:39:23  roboos
% added copyrigth and cvs logging, changed some comments
%

% this function is called many times in FieldTrip and associated toolboxes
% use efficient handling if the same toolbox has been investigated before
persistent previous
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
  'OPENMEEG'   'see http://gforge.inria.fr/projects/openmeeg'
  'PRTOOLS'    'see http://www.prtools.org'
  'LC-LIBS'    'contact Stefania Della Penna'
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
    status = exist('dss', 'file') && exist('dss_create_state', 'file');
  case 'EEGLAB'
    status = exist('runica', 'file');
  case 'NWAY'
    status = exist('parafac', 'file');
  case 'SPM99'
    status = exist('spm.m') && strcmp(spm('ver'),'SPM99');
  case 'SPM2'
    status = exist('spm.m') && strcmp(spm('ver'),'SPM2');
  case 'SPM5'
    status = exist('spm.m') && strcmp(spm('ver'),'SPM5');
  case 'SPM8'
    status = exist('spm.m') && strncmp(spm('ver'),'SPM8', 3);
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
    status  = (exist('GetMeg160ChannelInfoM') && exist('GetMeg160ContinuousRawDataM'));
  case 'BEOWULF'
    status = (exist('evalwulf') && exist('evalwulf') && exist('evalwulf'));
  case 'MENTAT'
    status  = (exist('pcompile') && exist('pfor') && exist('peval'));
  case 'SON2'
    status  = (exist('SONFileHeader') && exist('SONChanList') && exist('SONGetChannel'));
  case '4D-VERSION'
    status  = (exist('read4d') && exist('read4dhdr'));
  case 'SIGNAL'
    status = hasfunction('medfilt1', toolbox); % also check the availability of a toolbox license
  case 'OPTIM'
    status  = hasfunction('fmincon', toolbox) && hasfunction('fminunc', toolbox); % also check the availability of a toolbox license
  case 'SPLINES'
    status  = hasfunction('bspline', toolbox) && hasfunction('csape', toolbox); % also check the availability of a toolbox license
  case 'IMAGE'
    status = hasfunction('bwlabeln', toolbox); % also check the availability of a toolbox license
  case 'FASTICA'
    status  = exist('fastica', 'file');
  case 'BRAINSTORM'
    status  = exist('bem_xfer');
  case 'FILEIO'
    status  = (exist('read_header') && exist('read_data') && exist('read_event') && exist('read_sens'));
  case 'FORWINV'
    status  = (exist('compute_leadfield') && exist('prepare_vol_sens'));
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
    status  = (exist('plot_topo', 'file') && exist('plot_mesh', 'file') && exist('plot_matrix', 'file'));
  case 'PRTOOLS'
    status  = (exist('prversion', 'file') && exist('dataset', 'file') && exist('svc', 'file'));
  case 'LC-LIBS'
    status  = (exist('lcReadHeader', 'file') && exist('lcReadData', 'file'));
  otherwise
    if ~silent, warning(sprintf('cannot determine whether the %s toolbox is present', toolbox)); end
    status = 0;
end

% it should be a boolean value
status = (status~=0);

% try to determine the path of the requested toolbox
if autoadd && ~status

  % for core fieldtrip modules
  prefix = fileparts(which('preprocessing'));
  if ~status
    status = myaddpath(fullfile(prefix, lower(toolbox)), silent);
  end

  % for external fieldtrip modules
  prefix = fullfile(fileparts(which('preprocessing')), 'external');
  if ~status
    status = myaddpath(fullfile(prefix, lower(toolbox)), silent);
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
    error(msg);
  end
end

% this function is called many times in FieldTrip and associated toolboxes
% use efficient handling if the same toolbox has been investigated before
if status
  previous.(fixname(toolbox)) = status;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = myaddpath(toolbox, silent)
if exist(toolbox, 'dir')
  if ~silent, warning(sprintf('adding %s toolbox to your Matlab path', toolbox)); end
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

