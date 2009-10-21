function [type] = filetype(filename, desired, varargin)

% FILETYPE determines the filetype of many EEG/MEG/MRI data files by
% looking at the name, extension and optionally (part of) its contents.
% It tries to determine the global type of file (which usually
% corresponds to the manufacturer, the recording system or to the
% software used to create the file) and the particular subtype (e.g.
% continuous, average).
%
% Use as
%   type = filetype(filename)
%   type = filetype(dirname)
%
% This gives you a descriptive string with the data type, and can be
% used in a switch-statement. The descriptive string that is returned
% usually is something like 'XXX_YYY'/ where XXX refers to the
% manufacturer and YYY to the type of the data.
%
% Alternatively, use as
%   flag = filetype(filename, type)
%   flag = filetype(dirname, type)
% This gives you a boolean flag (0 or 1) indicating whether the file
% is of the desired type, and can be used to check whether the
% user-supplied file is what your subsequent code expects.
%
% Alternatively, use as
%   flag = filetype(dirlist, type)
% where the dirlist contains a list of files contained within one
% directory. This gives you a boolean vector indicating for each file
% whether it is of the desired type.
%
% Most filetypes of the following manufacturers and/or software programs are recognized
%  - VSMMedtech/CTF
%  - Elektra/Neuromag
%  - Yokogawa
%  - 4D/BTi
%  - EDF
%  - Neuroscan
%  - Analyse
%  - EEProbe
%  - BrainVision
%  - BESA
%  - Curry
%  - ASA
%  - LORETA
%  - Analyze/SPM
%  - MINC
%  - AFNI
%  - Neuralynx
%  - Plexon
%
% See also READ_XXX_YYY where XXX=manufacturer and YYY=subtype

% Copyright (C) 2003-2007 Robert Oostenveld
%
% $Log: filetype.m,v $
% Revision 1.101  2009/10/16 07:31:18  roboos
% renamed chieti into itab for consistency with other formats
%
% Revision 1.100  2009/10/13 10:12:51  roboos
% added support for chieti_raw
%
% Revision 1.99  2009/10/07 09:48:00  jansch
% added some lines to be able to deal with memory mapped files
%
% Revision 1.98  2009/07/07 10:38:38  roboos
% added header check for dicom
%
% Revision 1.97  2009/05/06 15:42:06  roboos
% also remember/check previous directory and don't do filetype caching in case type=unknown
%
% Revision 1.96  2009/04/01 06:53:23  roboos
% implemented caching for the type detection in case the same input is given multiple times
% this uses a persistent variable
% renamed the internal variable "ftype" into "type"
% removed the second output argument with the details, because sofar they have never been used in fieldtrip and therefore don't seem needed
%
% Revision 1.95  2009/02/12 11:47:23  vlalit
% Added support for neuro prax (eldith) EEG format based on functions from the manufacturer
%  used with permission from the company's representative Mr. Klaus Schellhorn.
%
% Revision 1.94  2009/02/09 14:21:00  roboos
% added inport of micromed_trc data
%
% Revision 1.93  2009/01/07 10:37:55  roboos
% changed description for neuralynx_sdma, should not have any functional consequences
%
% Revision 1.92  2008/12/24 10:33:48  roboos
% fixed header for mbfys_ama
%
% Revision 1.91  2008/12/19 14:39:25  marvger
% added support for udp, tcp and fifo
%
% Revision 1.90  2008/12/18 11:24:30  vlalit
% Fixed detection of spike6 matlab file
%
% Revision 1.89  2008/12/09 16:50:41  roboos
% disabled detection of ced_spike6mat, because it contains an error. The error should be fixed by Vladimir.
%
% Revision 1.88  2008/11/28 10:24:24  roboos
% added ctf mri version 4, thanks to Ivar
%
% Revision 1.87  2008/11/12 16:50:26  roboos
% improved detection for BCI2000
%
% Revision 1.86  2008/11/02 10:38:44  roboos
% added another check for ctf_ds (to deal with '.')
%
% Revision 1.85  2008/10/28 16:08:14  roboos
% additional check on yokogawa raw (first 4 bytes should be zero)
%
% Revision 1.84  2008/07/24 08:43:40  roboos
% added nimh_cortex
%
% Revision 1.83  2008/07/09 12:16:38  roboos
% added mclust_t
%
% Revision 1.82  2008/07/01 13:35:25  roboos
% moved long sequence for mat files onto single line (2x)
%
% Revision 1.81  2008/06/18 08:24:33  roboos
% added support for BCI2000
%
% Revision 1.80  2008/05/27 16:12:26  vlalit
% Changed type name to ced_spike6mat
%
% Revision 1.79  2008/05/27 11:58:20  vlalit
% Added support of Matlab files exported from Spike 6
%
% Revision 1.78  2008/05/21 13:31:34  vlalit
% Added check for existence of the file before the other spmeeg checks.
%
% Revision 1.77  2008/05/09 17:04:41  vlalit
% Added even more stringent criteria for recognition of SPM EEG mat files
%
% Revision 1.76  2008/05/06 14:36:23  roboos
% Added readers for SPM5 and SPM8 EEG formats (fix cvs conflict between jansch and vlalit)
%
% Revision 1.75  2008/05/06 12:00:38  jansch
% changed strcmp into strfind for 4d-data
%
% Revision 1.74  2008/05/06 09:24:36  jansch
% removed check for . in filename of 4D-files. necessary to read in data that
% are not DC-recorded.
%
% Revision 1.72  2008/04/21 11:50:51  roboos
% added support for egi_sbin, thanks to Joseph Dien
%
% Revision 1.71  2008/04/18 14:08:03  roboos
% added eeglab_set
%
% Revision 1.70  2008/04/11 12:09:18  roboos
% added polhemus_fil
%
% Revision 1.69  2008/03/04 11:14:26  roboos
% added neuralynx_nst
%
% Revision 1.68  2008/02/19 10:07:34  roboos
% replaced tcpsocket with buffer
%
% Revision 1.67  2008/01/31 20:10:22  roboos
% It is necessary to identify when it is necessary to execute the new
% read_4d_hdr.m function.  To achieve this an extra IF statement was
% added from line 312.  This statement checks for a passed filename
% which has no extension (not the case for .xyz or .m4d) and that the
% second to fourth characters in the filename are ,rf.  If both of
% these conditions are met then the filetype and manufacturer are set
% to 4d. [thanks to Gavin]
%
% Revision 1.66  2008/01/25 07:41:57  roboos
% moved some code, no functional change
%
% Revision 1.65  2008/01/15 08:12:30  roboos
% more strict test for ctf_ds
%
% Revision 1.64  2007/12/19 07:17:26  roboos
% removed ftc, added txt
%
% Revision 1.63  2007/12/17 16:14:44  roboos
% added neuralynx_bin
%
% Revision 1.62  2007/12/17 08:25:04  roboos
% added nexstim_nxe
%
% Revision 1.61  2007/12/12 16:48:20  roboos
% deal with case of extension for nev/Nev
% fixed bug: correctly locate the events.nev file in the dataset directory
%
% Revision 1.60  2007/12/12 14:39:26  roboos
% added riff_wave
%
% Revision 1.59  2007/12/06 17:05:54  roboos
% added fcdc_ftc, only on file extension
%
% Revision 1.58  2007/10/25 12:47:31  roboos
% added *.nrd as neuralynx_dma, not sure yet whether this is correct
%
% Revision 1.57  2007/10/16 12:35:08  roboos
% implemented fcdc_global
%
% Revision 1.56  2007/10/15 16:02:46  roboos
% added rfb://<password>@<host>:<port>
%
% Revision 1.55  2007/10/04 11:55:40  roboos
% added the various spass file formats
%
% Revision 1.54  2007/09/13 09:46:51  roboos
% adedd ctf_wts and svl
%
% Revision 1.53  2007/08/06 09:07:33  roboos
% added 4d_hs
%
% Revision 1.52  2007/07/27 12:17:47  roboos
% added ctf_shm
%
% Revision 1.51  2007/07/04 13:20:51  roboos
% added support for egi_egis/egia, thanks to Joseph Dien
%
% Revision 1.50  2007/06/06 07:10:16  roboos
% changed detection of BCI streams into using URI scheme
% changed deterction of BESA source waveform files
%
% Revision 1.49  2007/06/04 18:24:49  roboos
% added nifti
%
% Revision 1.48  2007/05/31 09:13:28  roboos
% added tcpsocket ans serial port
%
% Revision 1.47  2007/05/15 14:59:54  roboos
% try to separate besa *.dat from brainvision *.dat
%
% Revision 1.46  2007/04/25 15:27:28  roboos
% rmoved besa_gen, instead added besa_sb (simple binary), which comes with the *.gen or *.generic ascii header file
%
% Revision 1.45  2007/03/21 17:21:30  roboos
% remove . and .. from the file listing in case of a directory as input
% removed neuralynx_nte, the correct extension is *.nts
% added header check to neuralynx_nts
% implemented subfunction most, c.f. any
% swiched from using any(...) to most(...) for determining content of dataset directory
% implemented plexon_ds for directory with nex files in it
% made some additional small changes
%
% Revision 1.44  2007/03/19 16:52:37  roboos
% added neuralynx_nte
%
% Revision 1.43  2007/01/09 09:29:25  roboos
% small change
%
% Revision 1.42  2007/01/04 08:12:00  roboos
% fixed bug for besa_avr, renamed an incorrect tag into plexon_plx
%

% these are for remembering the type on subsequent calls with the same input arguments
persistent previous_argin previous_argout previous_pwd

if nargin<2
  % ensure that all input arguments are defined
  desired = [];
end

current_argin = {filename, desired, varargin{:}};
current_pwd   = pwd;
if isequal(current_argin, previous_argin) && isequal(current_pwd, previous_pwd)
  % don't do the detection again, but return the previous value from cache
  type = previous_argout{1};
  return
end

if strcmp(class(filename), 'memmapfile'),
  filename = filename.Filename;
end

% % get the optional arguments
% checkheader = keyval('checkheader', varargin); if isempty(checkheader), checkheader=1; end
%
% if ~checkheader
%   % assume that the header is always ok, e.g when the file does not yet exist
%   % this replaces the normal function with a function that always returns true
%   filetype_check_header = @filetype_true;
% end

if iscell(filename)
  % perform the test for each filename, return a boolean vector
  type = false(size(filename));
  for i=1:length(filename)
    if strcmp(filename{i}(end), '.')
      % do not recurse into this directory or the parent directory
      continue
    else
      type(i) = filetype(filename{i}, desired);
    end
  end
  return
end

% start with unknown values
type         = 'unknown';
manufacturer = 'unknown';
content      = 'unknown';

[p, f, x] = fileparts(filename);

if isdir(filename)
  % the directory listing is needed below
  ls = dir(filename);
  % remove the parent directory and the directory itself from the list
  ls = ls(~strcmp({ls.name}, '.'));
  ls = ls(~strcmp({ls.name}, '..'));
  for i=1:length(ls)
    % make sure that the directory listing includes the complete path
    ls(i).name = fullfile(filename, ls(i).name);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start determining the filetype
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these are some streams for asynchronous BCI
if filetype_check_uri(filename, 'fifo')
  type        = 'fcdc_fifo';
  manufacturer = 'F.C. Donders Centre';
  content      = 'stream';
elseif filetype_check_uri(filename, 'buffer')
  type        = 'fcdc_buffer';
  manufacturer = 'F.C. Donders Centre';
  content      = 'stream';
elseif filetype_check_uri(filename, 'mysql')
  type        = 'fcdc_mysql';
  manufacturer = 'F.C. Donders Centre';
  content      = 'stream';
elseif filetype_check_uri(filename, 'tcp')
  type        = 'fcdc_tcp';
  manufacturer = 'F.C. Donders Centre';
  content      = 'stream';
elseif filetype_check_uri(filename, 'udp')
  type        = 'fcdc_udp';
  manufacturer = 'F.C. Donders Centre';
  content      = 'stream';
elseif filetype_check_uri(filename, 'rfb')
  type        = 'fcdc_rfb';
  manufacturer = 'F.C. Donders Centre';
  content      = 'stream';
elseif filetype_check_uri(filename, 'serial')
  type        = 'fcdc_serial';
  manufacturer = 'F.C. Donders Centre';
  content      = 'stream';
elseif filetype_check_uri(filename, 'global')
  type        = 'fcdc_global';
  manufacturer = 'F.C. Donders Centre';
  content      = 'global variable';
elseif filetype_check_uri(filename, 'shm')
  type        = 'ctf_shm';
  manufacturer = 'CTF';
  content      = 'real-time shared memory buffer';

  % known CTF file types
elseif isdir(filename) && filetype_check_extension(filename, '.ds') && exist(fullfile(filename, [f '.res4']))
  type = 'ctf_ds';
  manufacturer = 'CTF';
  content = 'MEG dataset';
elseif isdir(filename) && ~isempty(dir(fullfile(filename, '*.res4'))) && ~isempty(dir(fullfile(filename, '*.meg4')))
  type = 'ctf_ds';
  manufacturer = 'CTF';
  content = 'MEG dataset';
elseif filetype_check_extension(filename, '.res4') && (filetype_check_header(filename, 'MEG41RS') || filetype_check_header(filename, 'MEG42RS') || filetype_check_header(filename, 'MEG4RES'))
  type = 'ctf_res4';
  manufacturer = 'CTF';
  content = 'MEG/EEG header information';
elseif filetype_check_extension(filename, '.meg4') && filetype_check_header(filename, 'MEG41CP')
  type = 'ctf_meg4';
  manufacturer = 'CTF';
  content = 'MEG/EEG';
elseif filetype_check_extension(filename, '.mrk') && filetype_check_header(filename, 'PATH OF DATASET:')
  type = 'ctf_mrk';
  manufacturer = 'CTF';
  content = 'marker file';
elseif filetype_check_extension(filename, '.mri') && filetype_check_header(filename, 'CTF_MRI_FORMAT VER 2.2')
  type = 'ctf_mri';
  manufacturer = 'CTF';
  content = 'MRI';
elseif filetype_check_extension(filename, '.mri') && filetype_check_header(filename, 'CTF_MRI_FORMAT VER 4', 31)
  type = 'ctf_mri4';
  manufacturer = 'CTF';
  content = 'MRI';
elseif filetype_check_extension(filename, '.hdm')
  type = 'ctf_hdm';
  manufacturer = 'CTF';
  content = 'volume conduction model';
elseif filetype_check_extension(filename, '.hc')
  type = 'ctf_hc';
  manufacturer = 'CTF';
  content = 'headcoil locations';
elseif filetype_check_extension(filename, '.shape')
  type = 'ctf_shape';
  manufacturer = 'CTF';
  content = 'headshape points';
elseif filetype_check_extension(filename, '.shape_info')
  type = 'ctf_shapeinfo';
  manufacturer = 'CTF';
  content = 'headshape information';
elseif filetype_check_extension(filename, '.wts')
  type = 'ctf_wts';
  manufacturer = 'CTF';
  content = 'SAM coefficients, i.e. spatial filter weights';
elseif filetype_check_extension(filename, '.svl')
  type = 'ctf_svl';
  manufacturer = 'CTF';
  content = 'SAM (pseudo-)statistic volumes';

  % known Micromed file types
elseif filetype_check_extension(filename, '.trc') && filetype_check_header(filename, '* MICROMED')
  type = 'micromed_trc';
  manufacturer = 'Micromed';
  content = 'Electrophysiological data';

  % known Neuromag file types
elseif filetype_check_extension(filename, '.fif')
  type = 'neuromag_fif';
  manufacturer = 'Neuromag';
  content = 'MEG header and data';
elseif filetype_check_extension(filename, '.bdip')
  type = 'neuromag_bdip';
  manufacturer = 'Neuromag';
  content = 'dipole model';

  % known Yokogawa file types
elseif filetype_check_extension(filename, '.ave') || filetype_check_extension(filename, '.sqd')
  type = 'yokogawa_ave';
  manufacturer = 'Yokogawa';
  content = 'averaged MEG data';
elseif filetype_check_extension(filename, '.con')
  type = 'yokogawa_con';
  manufacturer = 'Yokogawa';
  content = 'continuous MEG data';
elseif filetype_check_extension(filename, '.raw') && filetype_check_header(filename, [char(0) char(0) char(0) char(0)])
  type = 'yokogawa_raw';
  manufacturer = 'Yokogawa';
  content = 'evoked/trialbased MEG data';
elseif filetype_check_extension(filename, '.mri') && filetype_check_header(filename, '####')  % FIXME, not correct
  type = 'yokogawa_mri';
  manufacturer = 'Yokogawa';
  content = 'anatomical MRI';

elseif filetype_check_extension(filename, '.pdf') && filetype_check_header(filename, 'E|lk') % I am not sure whether this header always applies
  type = '4d_pdf';
  manufacturer = '4D/BTI';
  content = 'raw MEG data (processed data file)';
elseif exist([filename '.m4d'], 'file') && exist([filename '.xyz'], 'file') % these two ascii header files accompany the raw data
  type = '4d_pdf';
  manufacturer = '4D/BTI';
  content = 'raw MEG data (processed data file)';
elseif filetype_check_extension(filename, '.m4d') && exist([filename(1:(end-3)) 'xyz'], 'file') % these come in pairs
  type = '4d_m4d';
  manufacturer = '4D/BTI';
  content = 'MEG header information';
elseif filetype_check_extension(filename, '.xyz') && exist([filename(1:(end-3)) 'm4d'], 'file') % these come in pairs
  type = '4d_xyz';
  manufacturer = '4D/BTI';
  content = 'MEG sensor positions';
elseif isequal(f, 'hs_file') % the filename is "hs_file"
  type = '4d_hs';
  manufacturer = '4D/BTI';
  content = 'head shape';
elseif length(filename)>=4 && ~isempty(strfind(filename,',rf'))
  type = '4d';
  manufacturer = '4D/BTi';
  content = '';

  % known EEProbe file types
elseif filetype_check_extension(filename, '.cnt') && filetype_check_header(filename, 'RIFF')
  type = 'eep_cnt';
  manufacturer = 'EEProbe';
  content = 'EEG';
elseif filetype_check_extension(filename, '.avr') && filetype_check_header(filename, char([38 0 16 0]))
  type = 'eep_avr';
  manufacturer = 'EEProbe';
  content = 'ERP';
elseif filetype_check_extension(filename, '.trg')
  type = 'eep_trg';
  manufacturer = 'EEProbe';
  content = 'trigger information';
elseif filetype_check_extension(filename, '.rej')
  type = 'eep_rej';
  manufacturer = 'EEProbe';
  content = 'rejection marks';

  % known ASA file types
elseif filetype_check_extension(filename, '.elc')
  type = 'asa_elc';
  manufacturer = 'ASA';
  content = 'electrode positions';
elseif filetype_check_extension(filename, '.vol')
  type = 'asa_vol';
  manufacturer = 'ASA';
  content = 'volume conduction model';
elseif filetype_check_extension(filename, '.bnd')
  type = 'asa_bnd';
  manufacturer = 'ASA';
  content = 'boundary element model details';
elseif filetype_check_extension(filename, '.msm')
  type = 'asa_msm';
  manufacturer = 'ASA';
  content = 'ERP';
elseif filetype_check_extension(filename, '.msr')
  type = 'asa_msr';
  manufacturer = 'ASA';
  content = 'ERP';
elseif filetype_check_extension(filename, '.dip')
  % FIXME, can also be CTF dipole file
  type = 'asa_dip';
  manufacturer = 'ASA';
elseif filetype_check_extension(filename, '.mri')
  % FIXME, can also be CTF mri file
  type = 'asa_mri';
  manufacturer = 'ASA';
  content = 'MRI image header';
elseif filetype_check_extension(filename, '.iso')
  type = 'asa_iso';
  manufacturer = 'ASA';
  content = 'MRI image data';

  % known BCI2000 file types
elseif filetype_check_extension(filename, '.dat') && (filetype_check_header(filename, 'BCI2000') || filetype_check_header(filename, 'HeaderLen='))
  type = 'bci2000_dat';
  manufacturer = 'BCI2000';
  content = 'continuous EEG';

  % known Neuroscan file types
elseif filetype_check_extension(filename, '.avg') && filetype_check_header(filename, 'Version 3.0')
  type = 'ns_avg';
  manufacturer = 'Neuroscan';
  content = 'averaged EEG';
elseif filetype_check_extension(filename, '.cnt') && filetype_check_header(filename, 'Version 3.0')
  type = 'ns_cnt';
  manufacturer = 'Neuroscan';
  content = 'continuous EEG';
elseif filetype_check_extension(filename, '.eeg') && filetype_check_header(filename, 'Version 3.0')
  type = 'ns_eeg';
  manufacturer = 'Neuroscan';
  content = 'epoched EEG';

elseif filetype_check_extension(filename, '.eeg') && filetype_check_header(filename, 'V3.0')
  type = 'neuroprax_eeg';
  manufacturer = 'eldith GmbH';
  content = 'continuous EEG';
elseif filetype_check_extension(filename, '.ee_')
  type = 'neuroprax_mrk';
  manufacturer = 'eldith GmbH';
  content = 'EEG markers';

  % known Analyze & SPM file types
elseif filetype_check_extension(filename, '.hdr')
  type = 'analyze_hdr';
  manufacturer = 'Mayo Analyze';
  content = 'PET/MRI image header';
elseif filetype_check_extension(filename, '.img')
  type = 'analyze_img';
  manufacturer = 'Mayo Analyze';
  content = 'PET/MRI image data';
elseif filetype_check_extension(filename, '.mnc')
  type = 'minc';
  content = 'MRI image data';
elseif filetype_check_extension(filename, '.nii')
  type = 'nifti';
  content = 'MRI image data';

  % known LORETA file types
elseif filetype_check_extension(filename, '.lorb')
  type = 'loreta_lorb';
  manufacturer = 'old LORETA';
  content = 'source reconstruction';
elseif filetype_check_extension(filename, '.slor')
  type = 'loreta_slor';
  manufacturer = 'sLORETA';
  content = 'source reconstruction';

  % known AFNI file types
elseif filetype_check_extension(filename, '.brik') || filetype_check_extension(filename, '.BRIK')
  type = 'afni_brik';
  content = 'MRI image data';
elseif filetype_check_extension(filename, '.head') || filetype_check_extension(filename, '.HEAD')
  type = 'afni_head';
  content = 'MRI header data';

  % known BrainVison file types
elseif filetype_check_extension(filename, '.vhdr')
  type = 'brainvision_vhdr';
  manufacturer = 'BrainProducts';
  content = 'EEG header';
elseif filetype_check_extension(filename, '.vmrk')
  type = 'brainvision_vmrk';
  manufacturer = 'BrainProducts';
  content = 'EEG markers';
elseif filetype_check_extension(filename, '.vabs')
  type = 'brainvision_vabs';
  manufacturer = 'BrainProducts';
  content = 'Brain Vison Analyzer macro';
elseif filetype_check_extension(filename, '.eeg')
  % FIXME, can also be Neuroscan epoched EEG data
  type = 'brainvision_eeg';
  manufacturer = 'BrainProducts';
  content = 'continuous EEG data';
elseif filetype_check_extension(filename, '.seg')
  type = 'brainvision_seg';
  manufacturer = 'BrainProducts';
  content = 'segmented EEG data';
elseif filetype_check_extension(filename, '.dat') && ~filetype_check_header(filename, 'HeaderLen=') && ~filetype_check_header(filename, 'BESA_SA_IMAGE') && ~(exist(fullfile(p, [f '.gen']), 'file') || exist(fullfile(p, [f '.generic']), 'file'))
  % WARNING this is a very general name, it could be exported BrainVision
  % data but also a BESA beamformer source reconstruction or BCI2000
  type = 'brainvision_dat';
  manufacturer = 'BrainProducts';
  content = 'exported EEG data';
elseif filetype_check_extension(filename, '.marker')
  type = 'brainvision_marker';
  manufacturer = 'BrainProducts';
  content = 'rejection markers';

  % known Polhemus file types
elseif filetype_check_extension(filename, '.pos')
  type = 'polhemus_pos';
  manufacturer = 'BrainProducts/CTF/Polhemus?'; % actually I don't know whose software it is
  content = 'electrode positions';

  % known Neuralynx file types
elseif filetype_check_extension(filename, '.nev') || filetype_check_extension(filename, '.Nev')
  type = 'neuralynx_nev';
  manufacturer = 'Neuralynx';
  content = 'event information';
elseif filetype_check_extension(filename, '.ncs') && filetype_check_header(filename, '####')
  type = 'neuralynx_ncs';
  manufacturer = 'Neuralynx';
  content = 'continuous single channel recordings';
elseif filetype_check_extension(filename, '.nse') && filetype_check_header(filename, '####')
  type = 'neuralynx_nse';
  manufacturer = 'Neuralynx';
  content = 'spike waveforms';
elseif filetype_check_extension(filename, '.nts')  && filetype_check_header(filename, '####')
  type = 'neuralynx_nts';
  manufacturer = 'Neuralynx';
  content = 'timestamps only';
elseif filetype_check_extension(filename, '.nvt')
  type = 'neuralynx_nvt';
  manufacturer = 'Neuralynx';
  content = 'video tracker';
elseif filetype_check_extension(filename, '.nst')
  type = 'neuralynx_nst';
  manufacturer = 'Neuralynx';
  content = 'continuous stereotrode recordings';
elseif filetype_check_extension(filename, '.ntt')
  type = 'neuralynx_ntt';
  manufacturer = 'Neuralynx';
  content = 'continuous tetrode recordings';
elseif strcmpi(f, 'logfile') && strcmpi(x, '.txt')  % case insensitive
  type = 'neuralynx_log';
  manufacturer = 'Neuralynx';
  content = 'log information in ASCII format';
elseif ~isempty(strfind(lower(f), 'dma')) && strcmpi(x, '.log')  % this is not a very strong detection
  type = 'neuralynx_dma';
  manufacturer = 'Neuralynx';
  content = 'raw aplifier data directly from DMA';
elseif filetype_check_extension(filename, '.nrd') % see also above, since Cheetah 5.x the file extension has changed
  type = 'neuralynx_dma';
  manufacturer = 'Neuralynx';
  content = 'raw aplifier data directly from DMA';
elseif isdir(filename) && (any(filetype_check_extension({ls.name}, '.nev')) || any(filetype_check_extension({ls.name}, '.Nev')))
  % a regular Neuralynx dataset directory that contains an event file
  type = 'neuralynx_ds';
  manufacturer = 'Neuralynx';
  content = 'dataset';
elseif isdir(filename) && most(filetype_check_extension({ls.name}, '.ncs'))
  % a directory containing continuously sampled channels in Neuralynx format
  type = 'neuralynx_ds';
  manufacturer = 'Neuralynx';
  content = 'continuously sampled channels';
elseif isdir(filename) && most(filetype_check_extension({ls.name}, '.nse'))
  % a directory containing spike waveforms in Neuralynx format
  type = 'neuralynx_ds';
  manufacturer = 'Neuralynx';
  content = 'spike waveforms';
elseif isdir(filename) && most(filetype_check_extension({ls.name}, '.nte'))
  % a directory containing spike timestamps in Neuralynx format
  type = 'neuralynx_ds';
  manufacturer = 'Neuralynx';
  content = 'spike timestamps';

  % these are formally not Neuralynx file formats, but at the FCDC we use them together with Neuralynx
elseif isdir(filename) && any(filetype({ls.name}, 'neuralynx_ds'))
  % a downsampled Neuralynx DMA file can be split into three seperate lfp/mua/spike directories
  % treat them as one combined dataset
  type = 'neuralynx_cds';
  manufacturer = 'F.C. Donders Centre';
  content = 'dataset containing seperate lfp/mua/spike directories';
elseif filetype_check_extension(filename, '.tsl') && filetype_check_header(filename, 'tsl')
  type = 'neuralynx_tsl';
  manufacturer = 'F.C. Donders Centre';
  content = 'timestamps from DMA log file';
elseif filetype_check_extension(filename, '.tsh') && filetype_check_header(filename, 'tsh')
  type = 'neuralynx_tsh';
  manufacturer = 'F.C. Donders Centre';
  content = 'timestamps from DMA log file';
elseif filetype_check_extension(filename, '.ttl') && filetype_check_header(filename, 'ttl')
  type = 'neuralynx_ttl';
  manufacturer = 'F.C. Donders Centre';
  content = 'Parallel_in from DMA log file';
elseif filetype_check_extension(filename, '.bin') && filetype_check_header(filename, {'uint8', 'uint16', 'uint32', 'int8', 'int16', 'int32', 'int64', 'float32', 'float64'})
  type = 'neuralynx_bin';
  manufacturer = 'F.C. Donders Centre';
  content = 'single channel continuous data';
elseif isdir(filename) && any(filetype_check_extension({ls.name}, '.ttl')) && any(filetype_check_extension({ls.name}, '.tsl')) && any(filetype_check_extension({ls.name}, '.tsh'))
  % a directory containing the split channels from a DMA logfile
  type = 'neuralynx_sdma';
  manufacturer = 'F.C. Donders Centre';
  content = 'split DMA log file';
elseif isdir(filename) && filetype_check_extension(filename, '.sdma')
  % a directory containing the split channels from a DMA logfile
  type = 'neuralynx_sdma';
  manufacturer = 'F.C. Donders Centre';
  content = 'split DMA log file';

  % known Plexon file types
elseif filetype_check_extension(filename, '.nex')  && filetype_check_header(filename, 'NEX1')
  type = 'plexon_nex';
  manufacturer = 'Plexon';
  content = 'electrophysiological data';
elseif filetype_check_extension(filename, '.plx')  && filetype_check_header(filename, 'PLEX')
  type = 'plexon_plx';
  manufacturer = 'Plexon';
  content = 'electrophysiological data';
elseif filetype_check_extension(filename, '.ddt')
  type = 'plexon_ddt';
  manufacturer = 'Plexon';
elseif isdir(filename) && most(filetype_check_extension({ls.name}, '.nex')) && most(filetype_check_header({ls.name}, 'NEX1'))
  % a directory containing multiple plexon NEX files
  type = 'plexon_ds';
  manufacturer = 'Plexon';
  content = 'electrophysiological data';

  % known Cambridge Electronic Design file types
elseif filetype_check_extension(filename, '.smr')
  type = 'ced_son';
  manufacturer = 'Cambridge Electronic Design';
  content = 'Spike2 SON filing system';

  % known BESA file types
elseif filetype_check_extension(filename, '.avr') && strcmp(type, 'unknown')
  type = 'besa_avr';  % FIXME, can also be EEProbe average EEG
  manufacturer = 'BESA';
  content = 'average EEG';
elseif filetype_check_extension(filename, '.elp')
  type = 'besa_elp';
  manufacturer = 'BESA';
  content = 'electrode positions';
elseif filetype_check_extension(filename, '.eps')
  type = 'besa_eps';
  manufacturer = 'BESA';
  content = 'digitizer information';
elseif filetype_check_extension(filename, '.sfp')
  type = 'besa_sfp';
  manufacturer = 'BESA';
  content = 'sensor positions';
elseif filetype_check_extension(filename, '.ela')
  type = 'besa_ela';
  manufacturer = 'BESA';
  content = 'sensor information';
elseif filetype_check_extension(filename, '.pdg')
  type = 'besa_pdg';
  manufacturer = 'BESA';
  content = 'paradigm file';
elseif filetype_check_extension(filename, '.tfc')
  type = 'besa_tfc';
  manufacturer = 'BESA';
  content = 'time frequency coherence';
elseif filetype_check_extension(filename, '.mul')
  type = 'besa_mul';
  manufacturer = 'BESA';
  content = 'multiplexed ascii format';
elseif filetype_check_extension(filename, '.dat') && filetype_check_header(filename, 'BESA_SA')  % header can start with BESA_SA_IMAGE or BESA_SA_MN_IMAGE
  type = 'besa_src';
  manufacturer = 'BESA';
  content = 'beamformer source reconstruction';
elseif filetype_check_extension(filename, '.swf') && filetype_check_header(filename, 'Npts=')
  type = 'besa_swf';
  manufacturer = 'BESA';
  content = 'beamformer source waveform';
elseif filetype_check_extension(filename, '.bsa')
  type = 'besa_bsa';
  manufacturer = 'BESA';
  content = 'beamformer source locations and orientations';
elseif exist(fullfile(p, [f '.dat']), 'file') && (exist(fullfile(p, [f '.gen']), 'file') || exist(fullfile(p, [f '.generic']), 'file'))
  type = 'besa_sb';
  manufacturer = 'BESA';
  content = 'simple binary channel data with a seperate generic ascii header';

  % old files from Pascal Fries' PhD research at the MPI
elseif filetype_check_extension(filename, '.dap') && filetype_check_header(filename, char(1))
  type = 'mpi_dap';
  manufacturer = 'MPI Frankfurt';
  content = 'electrophysiological data';
elseif isdir(filename) && ~isempty(cell2mat(regexp({ls.name}, '.dap$')))
  type = 'mpi_ds';
  manufacturer = 'MPI Frankfurt';
  content = 'electrophysiological data';

  % Frankfurt SPASS format, which uses the Labview Datalog (DTLG) format
elseif  filetype_check_header(filename, 'DTLG') && filetype_check_extension(filename, '.ana')
  type = 'spass_ana';
  manufacturer = 'MPI Frankfurt';
  content = 'electrophysiological data';
elseif  filetype_check_header(filename, 'DTLG') && filetype_check_extension(filename, '.swa')
  type = 'spass_swa';
  manufacturer = 'MPI Frankfurt';
  content = 'electrophysiological data';
elseif  filetype_check_header(filename, 'DTLG') && filetype_check_extension(filename, '.spi')
  type = 'spass_spi';
  manufacturer = 'MPI Frankfurt';
  content = 'electrophysiological data';
elseif  filetype_check_header(filename, 'DTLG') && filetype_check_extension(filename, '.stm')
  type = 'spass_stm';
  manufacturer = 'MPI Frankfurt';
  content = 'electrophysiological data';
elseif  filetype_check_header(filename, 'DTLG') && filetype_check_extension(filename, '.bhv')
  type = 'spass_bhv';
  manufacturer = 'MPI Frankfurt';
  content = 'electrophysiological data';

  % known Chieti ITAB file types
elseif filetype_check_extension(filename, '.raw') && filetype_check_header(filename, 'FORMAT: ATB-BIOMAGDATA')
  type = 'itab_raw';
  manufacturer = 'Chieti ITAB';
  content = 'MEG data, including sensor positions';

  % known Nexstim file types
elseif filetype_check_extension(filename, '.nxe')
  type = 'nexstim_nxe';
  manufacturer = 'Nexstim';
  content = 'electrophysiological data';

  % known Curry V4 file types
elseif filetype_check_extension(filename, '.dap')
  type = 'curry_dap';   % FIXME, can also be MPI Frankfurt electrophysiological data
  manufacturer = 'Curry';
  content = 'data parameter file';
elseif filetype_check_extension(filename, '.dat')
  type = 'curry_dat';
  manufacturer = 'Curry';
  content = 'raw data file';
elseif filetype_check_extension(filename, '.rs4')
  type = 'curry_rs4';
  manufacturer = 'Curry';
  content = 'sensor geometry file';
elseif filetype_check_extension(filename, '.par')
  type = 'curry_par';
  manufacturer = 'Curry';
  content = 'data or image parameter file';
elseif filetype_check_extension(filename, '.bd0') || filetype_check_extension(filename, '.bd1') || filetype_check_extension(filename, '.bd2') || filetype_check_extension(filename, '.bd3') || filetype_check_extension(filename, '.bd4') || filetype_check_extension(filename, '.bd5') || filetype_check_extension(filename, '.bd6') || filetype_check_extension(filename, '.bd7') || filetype_check_extension(filename, '.bd8') || filetype_check_extension(filename, '.bd9')
  type = 'curry_bd';
  manufacturer = 'Curry';
  content = 'BEM description file';
elseif filetype_check_extension(filename, '.bt0') || filetype_check_extension(filename, '.bt1') || filetype_check_extension(filename, '.bt2') || filetype_check_extension(filename, '.bt3') || filetype_check_extension(filename, '.bt4') || filetype_check_extension(filename, '.bt5') || filetype_check_extension(filename, '.bt6') || filetype_check_extension(filename, '.bt7') || filetype_check_extension(filename, '.bt8') || filetype_check_extension(filename, '.bt9')
  type = 'curry_bt';
  manufacturer = 'Curry';
  content = 'BEM transfer matrix file';
elseif filetype_check_extension(filename, '.bm0') || filetype_check_extension(filename, '.bm1') || filetype_check_extension(filename, '.bm2') || filetype_check_extension(filename, '.bm3') || filetype_check_extension(filename, '.bm4') || filetype_check_extension(filename, '.bm5') || filetype_check_extension(filename, '.bm6') || filetype_check_extension(filename, '.bm7') || filetype_check_extension(filename, '.bm8') || filetype_check_extension(filename, '.bm9')
  type = 'curry_bm';
  manufacturer = 'Curry';
  content = 'BEM full matrix file';
elseif filetype_check_extension(filename, '.dig')
  type = 'curry_dig';
  manufacturer = 'Curry';
  content = 'digitizer file';

  % known Curry V2 file types
elseif filetype_check_extension(filename, '.sp0') || filetype_check_extension(filename, '.sp1') || filetype_check_extension(filename, '.sp2') || filetype_check_extension(filename, '.sp3') || filetype_check_extension(filename, '.sp4') || filetype_check_extension(filename, '.sp5') || filetype_check_extension(filename, '.sp6') || filetype_check_extension(filename, '.sp7') || filetype_check_extension(filename, '.sp8') || filetype_check_extension(filename, '.sp9')
  type = 'curry_sp';
  manufacturer = 'Curry';
  content = 'point list';
elseif filetype_check_extension(filename, '.s10') || filetype_check_extension(filename, '.s11') || filetype_check_extension(filename, '.s12') || filetype_check_extension(filename, '.s13') || filetype_check_extension(filename, '.s14') || filetype_check_extension(filename, '.s15') || filetype_check_extension(filename, '.s16') || filetype_check_extension(filename, '.s17') || filetype_check_extension(filename, '.s18') || filetype_check_extension(filename, '.s19') || filetype_check_extension(filename, '.s20') || filetype_check_extension(filename, '.s21') || filetype_check_extension(filename, '.s22') || filetype_check_extension(filename, '.s23') || filetype_check_extension(filename, '.s24') || filetype_check_extension(filename, '.s25') || filetype_check_extension(filename, '.s26') || filetype_check_extension(filename, '.s27') || filetype_check_extension(filename, '.s28') || filetype_check_extension(filename, '.s29') || filetype_check_extension(filename, '.s30') || filetype_check_extension(filename, '.s31') || filetype_check_extension(filename, '.s32') || filetype_check_extension(filename, '.s33') || filetype_check_extension(filename, '.s34') || filetype_check_extension(filename, '.s35') || filetype_check_extension(filename, '.s36') || filetype_check_extension(filename, '.s37') || filetype_check_extension(filename, '.s38') || filetype_check_extension(filename, '.s39')
  type = 'curry_s';
  manufacturer = 'Curry';
  content = 'triangle or tetraedra list';
elseif filetype_check_extension(filename, '.pom')
  type = 'curry_pom';
  manufacturer = 'Curry';
  content = 'anatomical localization file';
elseif filetype_check_extension(filename, '.res')
  type = 'curry_res';
  manufacturer = 'Curry';
  content = 'functional localization file';

  % known MBFYS file types
elseif filetype_check_extension(filename, '.tri')
  type = 'mbfys_tri';
  manufacturer = 'MBFYS';
  content = 'triangulated surface';
elseif filetype_check_extension(filename, '.ama') && filetype_check_header(filename, [10 0 0 0])
  type = 'mbfys_ama';
  manufacturer = 'MBFYS';
  content = 'BEM volume conduction model';

  % Electrical Geodesics Incorporated format
elseif (filetype_check_extension(filename, '.egis') || filetype_check_extension(filename, '.ave') || filetype_check_extension(filename, '.gave') || filetype_check_extension(filename, '.raw')) && (filetype_check_header(filename, [char(1) char(2) char(3) char(4) char(255) char(255)]) || filetype_check_header(filename, [char(3) char(4) char(1) char(2) char(255) char(255)]))
  type = 'egi_egia';
  manufacturer = 'Electrical Geodesics Incorporated';
  content = 'averaged EEG data';
elseif (filetype_check_extension(filename, '.egis') || filetype_check_extension(filename, '.ses') || filetype_check_extension(filename, '.raw')) && (filetype_check_header(filename, [char(1) char(2) char(3) char(4) char(0) char(3)]) || filetype_check_header(filename, [char(3) char(4) char(1) char(2) char(0) char(3)]))
  type = 'egi_egis';
  manufacturer = 'Electrical Geodesics Incorporated';
  content = 'raw EEG data';
elseif (filetype_check_extension(filename, '.sbin') || filetype_check_extension(filename, '.raw'))
  % note that the Chieti MEG data format also has the extension *.raw
  % but that can be detected by looking at the file header
  type = 'egi_sbin';
  manufacturer = 'Electrical Geodesics Incorporated';
  content = 'averaged EEG data';


  % some other known file types
elseif length(filename)>4 && exist([filename(1:(end-4)) '.mat'], 'file') && exist([filename(1:(end-4)) '.bin'], 'file')
  % this is a self-defined FCDC data format, consisting of two files
  % there is a matlab V6 file with the header and a binary file with the data (multiplexed, ieee-le, double)
  type = 'fcdc_matbin';
  manufacturer = 'F.C. Donders Centre';
  content = 'multiplexed electrophysiology data';
elseif filetype_check_extension(filename, '.lay')
  type = 'layout';
  manufacturer = 'Ole Jensen';
  content = 'layout of channels for plotting';
elseif filetype_check_extension(filename, '.dcm') || filetype_check_extension(filename, '.ima') || filetype_check_header(filename, 'DICM', 128)
  type = 'dicom';
  manufacturer = 'Dicom';
  content = 'image data';
elseif filetype_check_extension(filename, '.trl')
  type = 'fcdc_trl';
  manufacturer = 'F.C.Donders';
  content = 'trial definitions';
elseif filetype_check_header(filename, [255 'BIOSEMI']) % filetype_check_extension(filename, '.bdf')
  type = 'biosemi_bdf';
  %   type = 'bham_bdf';
  manufacturer = 'Biosemi Data Format';
  content = 'electrophysiological data';
elseif filetype_check_extension(filename, '.edf')
  type = 'edf';
  manufacturer = 'European Data Format';
  content = 'electrophysiological data';
elseif filetype_check_extension(filename, '.mat') && exist(filename, 'file') && exist([filename(1:(end-4)) '.dat'], 'file') && numel(whos('-file', filename))==1 && strcmp('D', getfield(whos('-file', filename), {1}, 'name')) && strcmp('struct', getfield(whos('-file', filename), {1}, 'class'))
  type = 'spmeeg_mat';
  manufacturer = 'Wellcome Trust Centre for Neuroimaging, UCL, UK';
  content = 'electrophysiological data';
elseif filetype_check_extension(filename, '.mat') && exist(filename, 'file') && filetype_check_ced_spike6mat(filename)
  type = 'ced_spike6mat';
  manufacturer = 'Cambridge Electronic Design Limited';
  content = 'electrophysiological data';
elseif filetype_check_extension(filename, '.mat') && filetype_check_header(filename, 'MATLAB')
  type = 'matlab';
  manufacturer = 'Matlab';
  content = 'Matlab binary data';
elseif filetype_check_header(filename, 'RIFF', 0) && filetype_check_header(filename, 'WAVE', 8)
  type = 'riff_wave';
  manufacturer = 'Microsoft';
  content = 'audio';
elseif filetype_check_extension(filename, '.txt')
  type = 'ascii_txt';
  manufacturer = '';
  content = '';
elseif filetype_check_extension(filename, '.pol')
  type = 'polhemus_fil';
  manufacturer = 'Functional Imaging Lab, London, UK';
  content = 'headshape points';
elseif filetype_check_extension(filename, '.set')
  type = 'eeglab_set';
  manufacturer = 'Swartz Center for Computational Neuroscience, San Diego, USA';
  content = 'electrophysiological data';
elseif filetype_check_extension(filename, '.t') && filetype_check_header(filename, '%%BEGINHEADER')
  type = 'mclust_t';
  manufacturer = 'MClust';
  content = 'sorted spikes';
elseif filetype_check_header(filename, 26)
  type = 'nimh_cortex';
  manufacturer = 'NIMH Laboratory of Neuropsychology, http://www.cortex.salk.edu';
  content = 'events and eye channels';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finished determining the filetype
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(type, 'unknown')
  warning('could not determine filetype of %s', filename);
end

if ~isempty(desired)
  % return a boolean value instead of a descriptive string
  type = strcmp(type, desired);
end

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
current_argout = {type};
if isempty(previous_argin) && ~strcmp(type, 'unknown')
  previous_argin  = current_argin;
  previous_argout = current_argout;
  previous_pwd    = current_pwd;
else
  % don't remember in case unknown
  previous_argin  = [];
  previous_argout = [];
  previous_pwd    = [];
end

return % filetype main()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that helps in deciding whether a directory with files should
% be treated as a "dataset". This function returns a logical 1 (TRUE) if more
% than half of the element of a vector are nonzero number or are logical 1
% (TRUE).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = most(x)
x = x(find(~isnan(x(:))));
y = sum(x==0)<ceil(length(x)/2);

% SUBFUNCTION that always returns a true value
function y = filetype_true(varargin)
y = 1;

% SUBFUNCTION that checks for CED spike6 mat file
function res = filetype_check_ced_spike6mat(filename)
res = 1;
var = whos('-file', filename);

% Check whether all the variables in the file are structs (representing channels)
if ~all(strcmp('struct', unique({var(:).class})) == 1)
  res = 0;
  return;
end

var = load(filename, var(1).name);
var = struct2cell(var);

% Check whether the fields of the first struct have some particular names
fnames = {
  'title'
  'comment'
  'interval'
  'scale'
  'offset'
  'units'
  'start'
  'length'
  'values'
  'times'};

res = (numel(intersect(fieldnames(var{1}), fnames)) == 10);
