function [status] = hastoolbox(toolbox, autoadd, silent)

% HASTOOLBOX tests whether an external toolbox is installed. Optionally
% it will try to determine the path to the toolbox and install it
% automatically.
%
% Use as
%   [status] = hastoolbox(toolbox, autoadd, silent)

% Copyright (C) 2005-2008, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

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
  'OPENMEEG'   'see http://gforge.inria.fr/projects/openmeeg'
  'PRTOOLS'    'see http://www.prtools.org'
  'LC-LIBS'    'contact Stefania Della Penna'
  'BSMART'     'see http://www.brain-smart.org'
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
  case 'BSMART' 
    status  = exist('bsmart'); 
  otherwise
    if ~silent, warning(sprintf('cannot determine whether the %s toolbox is present', toolbox)); end
    status = 0;
end

% it should be a boolean value
status = (status~=0);

% try to determine the path of the requested toolbox
if autoadd && ~status

  % for core fieldtrip modules
  prefix = fileparts(which('fieldtripdefs'));
  if ~status
    status = myaddpath(fullfile(prefix, lower(toolbox)), silent);
  end

  % for external fieldtrip modules
  prefix = fullfile(fileparts(which('fieldtripdefs')), 'external');
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

% remember the previous path, allows us to determine on the next call
% whether the path has been modified outise of this function
previouspath = matlabpath;

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

