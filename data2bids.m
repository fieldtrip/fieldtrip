function cfg = data2bids(cfg, varargin)

% DATA2BIDS is a helper function to convert MEG, EEG, iEEG or MRI data to the Brain
% Imaging Data Structure. The overall idea is that you write a MATLAB script in which
% you call this function multiple times, once for each individually recorded data
% file (or data set). It will write the corresponding sidecar JSON and TSV files for
% each data file.
%
% Use as
%   data2bids(cfg)
% or as
%   data2bids(cfg, data)
%
% The first input argument 'cfg' is the configuration structure, which contains the
% details for the (meta)data and which specifies the sidecar files you want to write.
% The optional 'data' argument corresponds to preprocessed raw data according to
% FT_DATAYPE_RAW or an anatomical MRI according to FT_DATAYPE_VOLUME. The optional
% data input argument allows you to write preprocessed electrophysiological data
% and/or realigned and defaced anatomical MRI to disk.
%
% The configuration structure should contains
%   cfg.method       = string, can be 'decorate', 'convert' or 'copy', see below (default is automatic)
%   cfg.dataset      = string, filename of the input data
%   cfg.outputfile   = string, optional filename for the output data (default is automatic)
%   cfg.writejson    = string, 'yes', 'replace', 'merge' or 'no' (default = 'yes')
%   cfg.writetsv     = string, 'yes', 'replace', 'merge' or 'no' (default = 'yes')
%
% This function starts from existing data file on disk or from a FieldTrip compatible
% data structure in MATLAB memory that is passed as the second input argument.
% Depending on cfg.method it will add the sidecar files, copy the dataset and add
% sidecar files, or convert the dataset and add the sidecar files. Each of the
% methods is discussed here.
%
% DECORATE - data2bids will read the header details and events from the data and write
% the appropriate sidecar files alongside the existing dataset. You would use this to
% obtain the sidecar files for data files that are already in the BIDS organization.
%
% CONVERT - data2bids will read the input data (or use the specified input data) and
% write it to a new output file that is BIDS compliant. The output format is NIfTI
% for MRI data, and BrainVision for EEG and iEEG. Note that MEG data files are stored
% in BIDS in their native format and this function will NOT convert them for you.
%
% COPY - data2bids will copy the data from the input data file to the output data
% file, which renames it, but does not change its content. Furthermore, it will read
% the header details and events from the data and construct the appropriate sidecar
% files.
%
% Although you can explicitly specify cfg.outputfile yourself, it is recommended to
% use the following configuration options. This results in a BIDS compliant output
% directory and file name. With these options data2bids will also write, or if
% already present update the participants.tsv and scans.tsv files.
%   cfg.bidsroot                = string, top level directory for the BIDS output
%   cfg.sub                     = string, subject name
%   cfg.ses                     = string, optional session name
%   cfg.run                     = number, optional
%   cfg.task                    = string, task name is required for functional data
%   cfg.datatype                = string, can be any of 'FLAIR', 'FLASH', 'PD', 'PDT2', 'PDmap', 'T1map', 'T1rho', 'T1w', 'T2map', 'T2star', 'T2w', 'angio', 'bold', 'bval', 'bvec', 'channels', 'coordsystem', 'defacemask', 'dwi', 'eeg', 'epi', 'events', 'fieldmap', 'headshape', 'ieeg', 'inplaneT1', 'inplaneT2', 'magnitude', 'magnitude1', 'magnitude2', 'meg', 'phase1', 'phase2', 'phasediff', 'photo', 'physio', 'sbref', 'stim'
%   cfg.acq                     = string
%   cfg.ce                      = string
%   cfg.rec                     = string
%   cfg.dir                     = string
%   cfg.mod                     = string
%   cfg.echo                    = string
%   cfg.proc                    = string
%
% When specifying the output directory in cfg.bidsroot, you can also specify
% additional information to be added as extra columns in the participants.tsv and
% scans.tsv files. For example:
%   cfg.participant.age         = scalar
%   cfg.participant.sex         = string, 'm' or 'f'
%   cfg.scans.acq_time          = string, should be formatted according to  RFC3339 as '2019-05-22T15:13:38'
%   cfg.dataset_description     = structure with additional fields, see below
% In case any of these values is specified as empty (i.e. []) or as nan, it will be
% written to the tsv file as 'n/a'.
%
% For anatomical and functional MRI data you can specify cfg.dicomfile to read the
% detailed MRI scanner and sequence details from the header of that DICOM file. This
% will be used to fill in the details of the corresponding JSON file.
%   cfg.dicomfile               = string, filename of a matching DICOM file for header details (default = [])
%   cfg.deface                  = string, 'yes' or 'no' (default = 'no')
%
% You can specify cfg.events as a Nx3 matrix with the trial definition (see
% FT_DEFINETRIAL) or as a MATLAB table. When specified as table, the first three
% columns containing integer values corresponding to the begsample, endsample and
% offset, the additional colums can be of another type and have any name. If you do
% not specify the trial definition, the events will be read from the MEG/EEG/iEEG
% dataset. Events from the trial definition or from the data will be written to
% events.tsv.
%   cfg.events                  = trial definition (see FT_DEFINETRIAL) or event structure (see FT_READ_EVENT)
%
% If NBS Presentation was used in combination with another functional data type, you
% can specify cfg.presentationfile with the name of the presentation log file, which
% will be aligned with the data based on triggers (MEG/EEG/iEEG) or based on the
% volumes (fMRI). Events from the presentation log file will also be written to
% events.tsv. To indicate how triggers (in MEG/EEG/iEEG) or volumes (in fMRI) match
% the presentation events, you should specify the mapping between them.
%   cfg.presentationfile        = string, optional filename for the presentation log file
%   cfg.trigger.eventtype       = string (default = [])
%   cfg.trigger.eventvalue      = string or number
%   cfg.trigger.skip            = 'last'/'first'/'none'
%   cfg.presentation.eventtype  = string (default = [])
%   cfg.presentation.eventvalue = string or number
%   cfg.presentation.skip       = 'last'/'first'/'none'
%
% For EEG and iEEG data you can specify an electrode definition according to
% FT_DATATYPE_SENS as an "elec" field in the input data, or you can specify it as
% cfg.elec or you can specify a filename with electrode information.
%   cfg.elec                    = structure with electrode positions or filename, see FT_READ_SENS
%
% General BIDS options that apply to all data types are
%   cfg.InstitutionName             = string
%   cfg.InstitutionAddress          = string
%   cfg.InstitutionalDepartmentName = string
%   cfg.Manufacturer                = string
%   cfg.ManufacturersModelName      = string
%   cfg.DeviceSerialNumber          = string
%   cfg.SoftwareVersions            = string
%
% If you specify cfg.bidsroot, this function will also write the dataset_description.json
% file. You can specify the following fields
%   cfg.dataset_description.writesidecar        = string
%   cfg.dataset_description.Name                = string
%   cfg.dataset_description.BIDSVersion         = string
%   cfg.dataset_description.License             = string
%   cfg.dataset_description.Authors             = string or cell-array of strings
%   cfg.dataset_description.Acknowledgements    = string
%   cfg.dataset_description.HowToAcknowledge    = string
%   cfg.dataset_description.Funding             = string or cell-array of strings
%   cfg.dataset_description.ReferencesAndLinks  = string or cell-array of strings
%   cfg.dataset_description.DatasetDOI          = string
%
% General BIDS options that apply to all functional data types are
%   cfg.TaskName                    = string
%   cfg.TaskDescription             = string
%   cfg.Instructions                = string
%   cfg.CogAtlasID                  = string
%   cfg.CogPOID                     = string
%
% There are more BIDS options for the mri/meg/eeg/ieeg data type specific sidecars.
% Rather than listing them all here, please open this function in the MATLAB editor,
% and scroll down a bit to see what those are. In general the information in the JSON
% files is specified in CamelCase, whereas the information for TSV files is in
% lowercase.
%   cfg.mri.SomeOption              = string in CamelCase, please check the MATLAB code
%   cfg.meg.SomeOption              = string in CamelCase, please check the MATLAB code
%   cfg.eeg.SomeOption              = string in CamelCase, please check the MATLAB code
%   cfg.ieeg.SomeOption             = string in CamelCase, please check the MATLAB code
%   cfg.channels.someoption         = string in lowercase, please check the MATLAB code
%   cfg.events.someoption           = string in lowercase, please check the MATLAB code
%   cfg.coordsystem.someoption      = string in lowercase, please check the MATLAB code
%
% The implementation in this function corresponds to BIDS version 1.2.0. See
% https://bids-specification.readthedocs.io/ for the full specification and
% http://bids.neuroimaging.io/ for further details.
%
% See also FT_DATAYPE_RAW, FT_DATAYPE_VOLUME, FT_DATATYPE_SENS, FT_DEFINETRIAL,
% FT_PREPROCESSING, FT_READ_MRI, FT_READ_EVENT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Undocumented options exist for converting some other data types to BIDS:
% - EMG
% - EXG
% - audio
% - video
% - eyetracking
% - motion
% - physio
% - stim
%
% Most of these data types are currently (Aug 2019) not supported in the BIDS
% specification, but this function converts them in a very similar way as the
% officially supported data types.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2018-2019, Robert Oostenveld
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug

if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% ensure backward compatibility
cfg = ft_checkconfig(cfg, 'renamed', {'anat', 'mri'});
cfg = ft_checkconfig(cfg, 'renamedval', {'native', 'no', 'convert'});
cfg = ft_checkconfig(cfg, 'renamedval', {'native', 'yes', 'copy'});
cfg = ft_checkconfig(cfg, 'renamed', {'native', 'method'});
cfg = ft_checkconfig(cfg, 'renamed', {'mri.deface', 'deface'});
cfg = ft_checkconfig(cfg, 'renamed', {'mri.writesidecar', 'writejson'});
cfg = ft_checkconfig(cfg, 'renamed', {'meg.writesidecar', 'writejson'});
cfg = ft_checkconfig(cfg, 'renamed', {'eeg.writesidecar', 'writejson'});
cfg = ft_checkconfig(cfg, 'renamed', {'ieeg.writesidecar', 'writejson'});
cfg = ft_checkconfig(cfg, 'renamed', {'events.writesidecar', 'writejson'});
cfg = ft_checkconfig(cfg, 'renamed', {'channels.writesidecar', 'writejson'});
cfg = ft_checkconfig(cfg, 'renamed', {'electrodes.writesidecar', 'writejson'});
cfg = ft_checkconfig(cfg, 'renamed', {'coordsystem.writesidecar', 'writejson'});

% get the options and set the defaults
cfg.method                  = ft_getopt(cfg, 'method');                     % default is handled below
cfg.dataset                 = ft_getopt(cfg, 'dataset');
cfg.feedback                = ft_getopt(cfg, 'feedback', 'yes');
cfg.outputfile              = ft_getopt(cfg, 'outputfile');                 % default is handled below
cfg.presentationfile        = ft_getopt(cfg, 'presentationfile');           % full path to the NBS presentation log file, it will be read and parsed using FT_READ_EVENT
cfg.presentation            = ft_getopt(cfg, 'presentation');
cfg.presentation.eventtype  = ft_getopt(cfg.presentation, 'eventtype');
cfg.presentation.eventvalue = ft_getopt(cfg.presentation, 'eventvalue');
cfg.presentation.skip       = ft_getopt(cfg.presentation, 'skip',  'last');  % this is a sensible default for fMRI, for MEG one should probably do 'none'
cfg.trigger                 = ft_getopt(cfg, 'trigger');
cfg.trigger.eventtype       = ft_getopt(cfg.trigger, 'eventtype');
cfg.trigger.eventvalue      = ft_getopt(cfg.trigger, 'eventvalue');
cfg.trigger.skip            = ft_getopt(cfg.trigger, 'skip',  'none');

% these are used to construct the directory and file name
cfg.bidsroot  = ft_getopt(cfg, 'bidsroot');
cfg.sub       = ft_getopt(cfg, 'sub');
cfg.ses       = ft_getopt(cfg, 'ses');
cfg.task      = ft_getopt(cfg, 'task');
cfg.acq       = ft_getopt(cfg, 'acq');
cfg.ce        = ft_getopt(cfg, 'ce');
cfg.rec       = ft_getopt(cfg, 'rec');
cfg.dir       = ft_getopt(cfg, 'dir');
cfg.run       = ft_getopt(cfg, 'run');
cfg.mod       = ft_getopt(cfg, 'mod');
cfg.echo      = ft_getopt(cfg, 'echo');
cfg.proc      = ft_getopt(cfg, 'proc');
cfg.datatype  = ft_getopt(cfg, 'datatype');

if isempty(cfg.datatype)
  modality = {'meg', 'eeg', 'ieeg', 'emg', 'exg', 'audio', 'video', 'eyetracker', 'physio', 'stim', 'motion'};
  for i=1:numel(modality)
    if isfield(cfg, modality{i}) && ~isempty(cfg.(modality{i}))
      % the user specified modality-specific options, assume that the datatype matches
      cfg.datatype = modality{i};
      ft_notice('assuming that the datatype is %s', cfg.datatype);
      continue
    end
  end % for each modality
end

cfg.dicomfile     = ft_getopt(cfg, 'dicomfile');                      % get header details from the specified DICOM files
cfg.deface        = ft_getopt(cfg, 'deface', 'no');                   % whether to deface the anatomical MRI
cfg.writejson     = ft_getopt(cfg, 'writejson', 'merge');             % whether to write the json file
cfg.writetsv      = ft_getopt(cfg, 'writetsv', 'merge');              % whether to write the tsv file

cfg.mri           = ft_getopt(cfg, 'mri');
cfg.meg           = ft_getopt(cfg, 'meg');
cfg.eeg           = ft_getopt(cfg, 'eeg');
cfg.ieeg          = ft_getopt(cfg, 'ieeg');
cfg.emg           = ft_getopt(cfg, 'emg');
cfg.exg           = ft_getopt(cfg, 'exg');
cfg.audio         = ft_getopt(cfg, 'audio');
cfg.video         = ft_getopt(cfg, 'video');
cfg.eyetracker    = ft_getopt(cfg, 'eyetracker');
cfg.physio        = ft_getopt(cfg, 'physio');
cfg.stim          = ft_getopt(cfg, 'stim');
cfg.motion        = ft_getopt(cfg, 'motion');
cfg.channels      = ft_getopt(cfg, 'channels');
cfg.electrodes    = ft_getopt(cfg, 'electrodes');
cfg.events        = ft_getopt(cfg, 'events');     % this can contain the trial definition as Nx3 array, as table, or an event structure
cfg.coordsystem   = ft_getopt(cfg, 'coordsystem');

%% Dataset description

cfg.dataset_description                     = ft_getopt(cfg, 'dataset_description'                       );
cfg.dataset_description.writesidecar        = ft_getopt(cfg.dataset_description, 'writesidecar', 'yes'   );
cfg.dataset_description.Name                = ft_getopt(cfg.dataset_description, 'Name'                  ); % REQUIRED. Name of the dataset.
cfg.dataset_description.BIDSVersion         = ft_getopt(cfg.dataset_description, 'BIDSVersion', 1.2      ); % REQUIRED. The version of the BIDS standard that was used.
cfg.dataset_description.License             = ft_getopt(cfg.dataset_description, 'License'               ); % RECOMMENDED. What license is this dataset distributed under? The use of license name abbreviations is suggested for specifying a license. A list of common licenses with suggested abbreviations can be found in Appendix II.
cfg.dataset_description.Authors             = ft_getopt(cfg.dataset_description, 'Authors'               ); % OPTIONAL. List of individuals who contributed to the creation/curation of the dataset.
cfg.dataset_description.Acknowledgements    = ft_getopt(cfg.dataset_description, 'Acknowledgements'      ); % OPTIONAL. Text acknowledging contributions of individuals or institutions beyond those listed in Authors or Funding.
cfg.dataset_description.HowToAcknowledge    = ft_getopt(cfg.dataset_description, 'HowToAcknowledge'      ); % OPTIONAL. Instructions how researchers using this dataset should acknowledge the original authors. This field can also be used to define a publication that should be cited in publications that use the dataset.
cfg.dataset_description.Funding             = ft_getopt(cfg.dataset_description, 'Funding'               ); % OPTIONAL. List of sources of funding (grant numbers)
cfg.dataset_description.ReferencesAndLinks  = ft_getopt(cfg.dataset_description, 'ReferencesAndLinks'    ); % OPTIONAL. List of references to publication that contain information on the dataset, or links.
cfg.dataset_description.DatasetDOI          = ft_getopt(cfg.dataset_description, 'DatasetDOI'            ); % OPTIONAL. The Document Object Identifier of the dataset (not the corresponding paper).

%% Generic fields for all data types
cfg.TaskName                          = ft_getopt(cfg, 'TaskName'                    ); % REQUIRED. Name of the task (for resting state use the "rest" prefix). Different Tasks SHOULD NOT have the same name. The Task label is derived from this field by removing all non alphanumeric ([a-zA-Z0-9]) characters.
cfg.TaskDescription                   = ft_getopt(cfg, 'TaskDescription'             ); % OPTIONAL. Description of the task.
cfg.Instructions                      = ft_getopt(cfg, 'Instructions'                ); % OPTIONAL. Text of the instructions given to participants before the scan. This is not only important for behavioral or cognitive tasks but also in resting state paradigms (e.g. to distinguish between eyes open and eyes closed).
cfg.CogAtlasID                        = ft_getopt(cfg, 'CogAtlasID'                  ); % OPTIONAL. URL of the corresponding "Cognitive Atlas term that describes the task (e.g. Resting State with eyes closed ""http://www.cognitiveatlas.org/term/id/trm_54e69c642d89b""")
cfg.CogPOID                           = ft_getopt(cfg, 'CogPOID'                     ); % OPTIONAL. URL of the corresponding "CogPO term that describes the task (e.g. Rest "http://wiki.cogpo.org/index.php?title=Rest")
cfg.Manufacturer                      = ft_getopt(cfg, 'Manufacturer'                ); % OPTIONAL. Manufacturer of the recording system ("CTF", ""Neuromag/Elekta"", ""4D/BTi"", ""KIT/Yokogawa"", ""ITAB"", "KRISS", "Other")
cfg.ManufacturersModelName            = ft_getopt(cfg, 'ManufacturersModelName'      ); % OPTIONAL. Manufacturer's designation of the model (e.g. "CTF-275"). See "Appendix VII" with preferred names
cfg.DeviceSerialNumber                = ft_getopt(cfg, 'DeviceSerialNumber'          ); % OPTIONAL. The serial number of the equipment that produced the composite instances. A pseudonym can also be used to prevent the equipment from being identifiable, as long as each pseudonym is unique within the dataset.
cfg.SoftwareVersions                  = ft_getopt(cfg, 'SoftwareVersions'            ); % OPTIONAL. Manufacturer's designation of the acquisition software.
cfg.InstitutionName                   = ft_getopt(cfg, 'InstitutionName'             ); % OPTIONAL. The name of the institution in charge of the equipment that produced the composite instances.
cfg.InstitutionAddress                = ft_getopt(cfg, 'InstitutionAddress'          ); % OPTIONAL. The address of the institution in charge of the equipment that produced the composite instances.
cfg.InstitutionalDepartmentName       = ft_getopt(cfg, 'InstitutionalDepartmentName' ); % The department in the institution in charge of the equipment that produced the composite instances. Corresponds to DICOM Tag 0008, 1040 "Institutional Department Name".

%% MR Scanner Hardware
cfg.mri.StationName                   = ft_getopt(cfg.mri, 'StationName'                    ); % Institution defined name of the machine that produced the composite instances. Corresponds to DICOM Tag 0008, 1010 "Station Name"
cfg.mri.HardcopyDeviceSoftwareVersion = ft_getopt(cfg.mri, 'HardcopyDeviceSoftwareVersion'  ); % (Deprecated) Manufacturer's designation of the software of the device that created this Hardcopy Image (the printer). Corresponds to DICOM Tag 0018, 101A "Hardcopy Device Software Version".
cfg.mri.MagneticFieldStrength         = ft_getopt(cfg.mri, 'MagneticFieldStrength'          ); % Nominal field strength of MR magnet in Tesla. Corresponds to DICOM Tag 0018,0087 "Magnetic Field Strength" .
cfg.mri.ReceiveCoilName               = ft_getopt(cfg.mri, 'ReceiveCoilName'                ); % Information describing the receiver coil. Corresponds to DICOM Tag 0018, 1250 "Receive Coil Name", although not all vendors populate that DICOM Tag, in which case this field can be derived from an appropriate private DICOM field.
cfg.mri.ReceiveCoilActiveElements     = ft_getopt(cfg.mri, 'ReceiveCoilActiveElements'      ); % Information describing the active/selected elements of the receiver coil.  This doesn't correspond to a tag in the DICOM ontology. The vendor-defined terminology for active coil elements can go in this field. As an example, for Siemens, coil channels are typically not activated/selected individually, but rather  in pre-defined selectable "groups" of individual channels, and the list of the  "groups" of elements that are active/selected in any given scan populates  the "Coil String" entry in Siemen's private DICOM fields (e.g., "HEA;HEP" for the Siemens standard 32 ch coil when both the anterior and posterior groups are activated). This is a flexible field that can be used as most appropriate for a given vendor and coil to define the "active" coil elements. Since individual scans can sometimes not have the intended coil elements selected, it is preferable for this field to be populated directly from the DICOM for each individual scan, so that it can be used as a mechanism for checking that a given scan was collected with the intended coil elements selected.
cfg.mri.GradientSetType               = ft_getopt(cfg.mri, 'GradientSetType'                ); % It should be possible to infer the gradient coil from the scanner model. If not,e.g. because of a custom upgrade or use of a gradient insert set, then the specifications of the actual gradient coil should be reported independently.
cfg.mri.MRTransmitCoilSequence        = ft_getopt(cfg.mri, 'MRTransmitCoilSequence'         ); % This is a relevant field if a non-standard transmit coil is used. Corresponds to DICOM Tag 0018, 9049 "MR Transmit Coil Sequence".
cfg.mri.MatrixCoilMode                = ft_getopt(cfg.mri, 'MatrixCoilMode'                 ); % (If used) A method for reducing the number of independent channels by combining in analog the signals from multiple coil elements. There are typically different default modes when using un-accelerated or accelerated (e.g. GRAPPA, SENSE) imaging.
cfg.mri.CoilCombinationMethod         = ft_getopt(cfg.mri, 'CoilCombinationMethod'          ); % Almost all fMRI studies using phased-array coils use root-sum-of-squares (rSOS) combination, but other methods exist. The image reconstruction is changed by the coil combination method (as for the matrix coil mode above), so anything non-standard should be reported.

%% MR Sequence Specifics
cfg.mri.PulseSequenceType             = ft_getopt(cfg.mri, 'PulseSequenceType'              ); % A general description of the pulse sequence used for the scan (i.e. MPRAGE, Gradient Echo EPI, Spin Echo EPI, Multiband gradient echo EPI).
cfg.mri.ScanningSequence              = ft_getopt(cfg.mri, 'ScanningSequence'               ); % Description of the type of data acquired. Corresponds to DICOM Tag 0018, 0020 "Sequence Sequence".
cfg.mri.SequenceVariant               = ft_getopt(cfg.mri, 'SequenceVariant'                ); % Variant of the ScanningSequence. Corresponds to DICOM Tag 0018, 0021 "Sequence Variant".
cfg.mri.ScanOptions                   = ft_getopt(cfg.mri, 'ScanOptions'                    ); % Parameters of ScanningSequence. Corresponds to DICOM Tag 0018, 0022 "Scan Options".
cfg.mri.SequenceName                  = ft_getopt(cfg.mri, 'SequenceName'                   ); % Manufacturer's designation of the sequence name. Corresponds to DICOM Tag 0018, 0024 "Sequence Name".
cfg.mri.PulseSequenceDetails          = ft_getopt(cfg.mri, 'PulseSequenceDetails'           ); % Information beyond pulse sequence type that identifies the specific pulse sequence used (i.e. "Standard Siemens Sequence distributed with the VB17 software," "Siemens WIP ### version #.##," or "Sequence written by X using a version compiled on MM/DD/YYYY").
cfg.mri.NonlinearGradientCorrection   = ft_getopt(cfg.mri, 'NonlinearGradientCorrection'    ); % Boolean stating if the image saved  has been corrected for gradient nonlinearities by the scanner sequence.

%% MR In-Plane Spatial Encoding
cfg.mri.NumberShots                   = ft_getopt(cfg.mri, 'NumberShots'                    ); % The number of RF excitations need to reconstruct a slice or volume. Please mind that  this is not the same as Echo Train Length which denotes the number of lines of k-space collected after an excitation.
cfg.mri.ParallelReductionFactorInPlan = ft_getopt(cfg.mri, 'ParallelReductionFactorInPlane' ); % The parallel imaging (e.g, GRAPPA) factor. Use the denominator of the fraction of k-space encoded for each slice. For example, 2 means half of k-space is encoded. Corresponds to DICOM Tag 0018, 9069 "Parallel Reduction Factor In-plane".
cfg.mri.ParallelAcquisitionTechnique  = ft_getopt(cfg.mri, 'ParallelAcquisitionTechnique'   ); % The type of parallel imaging used (e.g. GRAPPA, SENSE). Corresponds to DICOM Tag 0018, 9078 "Parallel Acquisition Technique".
cfg.mri.PartialFourier                = ft_getopt(cfg.mri, 'PartialFourier'                 ); % The fraction of partial Fourier information collected. Corresponds to DICOM Tag 0018, 9081 "Partial Fourier".
cfg.mri.PartialFourierDirection       = ft_getopt(cfg.mri, 'PartialFourierDirection'        ); % The direction where only partial Fourier information was collected. Corresponds to DICOM Tag 0018, 9036 "Partial Fourier Direction".
cfg.mri.PhaseEncodingDirection        = ft_getopt(cfg.mri, 'PhaseEncodingDirection'         ); % Possible values = [];                     % "i", "j", "k", "i-", "j-", "k-". The letters "i", "j", "k" correspond to the first, second and third axis of the data in the NIfTI file. The polarity of the phase encoding is assumed to go from zero index to maximum index unless "-" sign is present (then the order is reversed - starting from the highest index instead of zero). PhaseEncodingDirection is defined as the direction along which phase is was modulated which may result in visible distortions. Note that this is not the same as the DICOM term InPlanePhaseEncodingDirection which can have "ROW" or "COL" values. This parameter is REQUIRED if corresponding fieldmap data is present or when using multiple runs with different phase encoding directions (which can be later used for field inhomogeneity correction).
cfg.mri.EffectiveEchoSpacing          = ft_getopt(cfg.mri, 'EffectiveEchoSpacing'           ); % The "effective" sampling interval, specified in seconds, between lines in the phase-encoding direction, defined based on the size of the reconstructed image in the phase direction.  It is frequently, but incorrectly, referred to as  "dwell time" (see DwellTime parameter below for actual dwell time).  It is  required for unwarping distortions using field maps. Note that beyond just in-plane acceleration, a variety of other manipulations to the phase encoding need to be accounted for properly, including partial fourier, phase oversampling, phase resolution, phase field-of-view and interpolation. This parameter is REQUIRED if corresponding fieldmap data is present.
cfg.mri.TotalReadoutTime              = ft_getopt(cfg.mri, 'TotalReadoutTime'               ); % This is actually the "effective" total readout time , defined as the readout duration, specified in seconds, that would have generated data with the given level of distortion.  It is NOT the actual, physical duration of the readout train.  If EffectiveEchoSpacing has been properly computed, it is just EffectiveEchoSpacing * (ReconMatrixPE - 1). . This parameter is REQUIRED if corresponding "field/distortion" maps acquired with opposing phase encoding directions are present  (see 8.9.4).

%% MR Timing Parameters
cfg.mri.EchoTime                      = ft_getopt(cfg.mri, 'EchoTime'                       ); % The echo time (TE) for the acquisition, specified in seconds. This parameter is REQUIRED if corresponding fieldmap data is present or the data comes from a multi echo sequence. Corresponds to DICOM Tag 0018, 0081 "Echo Time"  (please note that the DICOM term is in milliseconds not seconds).
cfg.mri.InversionTime                 = ft_getopt(cfg.mri, 'InversionTime'                  ); % The inversion time (TI) for the acquisition, specified in seconds. Inversion time is the time after the middle of inverting RF pulse to middle of excitation pulse to detect the amount of longitudinal magnetization. Corresponds to DICOM Tag 0018, 0082 "Inversion Time"  (please note that the DICOM term is in milliseconds not seconds).
cfg.mri.SliceTiming                   = ft_getopt(cfg.mri, 'SliceTiming'                    ); % The time at which each slice was acquired within each volume (frame) of  the acquisition.  Slice timing is not slice order -- rather, it  is a list of times (in JSON format) containing the time (in seconds) of each slice acquisition in relation to the beginning of volume acquisition.  The list goes through the slices along the slice axis in the slice encoding dimension (see below). Note that to ensure the proper interpretation of the SliceTiming field, it is important to check if the (optional) SliceEncodingDirection exists. In particular,  if SliceEncodingDirection is negative, the entries in SliceTiming are defined in reverse order with respect to the slice axis (i.e., the final entry in the SliceTiming list is the time of acquisition of slice 0). This parameter is REQUIRED for sparse sequences that do not have the DelayTime field set. In addition without this parameter slice time correction will not be possible.
cfg.mri.SliceEncodingDirection        = ft_getopt(cfg.mri, 'SliceEncodingDirection'         ); % Possible values = [];                     % "i", "j", "k", "i-", "j-", "k-" (the axis of the NIfTI data along which slices were acquired, and the direction in which SliceTiming is  defined with respect to). "i", "j", "k" identifiers correspond to the first, second and third axis of the data in the NIfTI file. A "-" sign indicates that the contents of SliceTiming are defined in reverse order -- that is, the first entry corresponds to the slice with the largest index, and the final entry corresponds to slice index zero. When present ,the axis defined by SliceEncodingDirection  needs to be consistent with the "slice_dim" field in the NIfTI header. When absent, the entries in SliceTiming must be in the order of increasing slice index as defined by the NIfTI header.
cfg.mri.DwellTime                     = ft_getopt(cfg.mri, 'DwellTime'                      ); % Actual dwell time (in seconds) of the receiver per point in the readout direction, including any oversampling.  For Siemens, this corresponds to DICOM field (0019,1018) (in ns).   This value is necessary for the (optional) readout distortion correction of anatomicals in the HCP Pipelines.  It also usefully provides a handle on the readout bandwidth, which isn't captured in the other metadata tags.  Not to be confused with "EffectiveEchoSpacing", and the frequent mislabeling of echo spacing (which is spacing in the phase encoding direction) as "dwell time" (which is spacing in the readout direction).

%% MR RF & Contrast
cfg.mri.FlipAngle                     = ft_getopt(cfg.mri, 'FlipAngle'                      ); % Flip angle for the acquisition, specified in degrees. Corresponds to = [];                     % DICOM Tag 0018, 1314 "Flip Angle".
cfg.mri.MultibandAccelerationFactor   = ft_getopt(cfg.mri, 'MultibandAccelerationFactor'    ); % RECOMMENDED. The multiband factor, for multiband acquisitions.
cfg.mri.NegativeContrast              = ft_getopt(cfg.mri, 'NegativeContrast'               ); % OPTIONAL. Boolean (true or false) value specifying whether increasing voxel intensity (within sample voxels) denotes a decreased value with respect to the contrast suffix. This is commonly the case when Cerebral Blood Volume is estimated via usage of a contrast agent in conjunction with a T2* weighted acquisition protocol.

%% MR Slice Acceleration
cfg.mri.MultibandAccelerationFactor   = ft_getopt(cfg.mri, 'MultibandAccelerationFactor'    ); % The multiband factor, for multiband acquisitions.

%% Anatomical landmarks, useful for multimodaltimodal co-registration with MEG, (S)HeadCoil, TMS,etc
cfg.mri.AnatomicalLandmarkCoordinates = ft_getopt(cfg.mri, 'AnatomicalLandmarkCoordinates'  ); % Key:value pairs of any number of additional anatomical landmarks and their coordinates in voxel units (where first voxel has index 0,0,0) relative to the associated anatomical MRI, (e.g. {"AC" = []; % [127,119,149], "PC" = []; % [128,93,141], "IH" = []; % [131,114,206]}, or {"NAS" = []; % [127,213,139], "LPA" = []; % [52,113,96], "RPA" = []; % [202,113,91]}).

%% MR Anatomical scan information
cfg.mri.ContrastBolusIngredient       = ft_getopt(cfg.mri, 'ContrastBolusIngredient'        ); % OPTIONAL. Active ingredient of agent. Values MUST be one of: IODINE, GADOLINIUM, CARBON DIOXIDE, BARIUM, XENON Corresponds to DICOM Tag 0018,1048.

%% MR Functional scan information
cfg.mri.RepetitionTime                = ft_getopt(cfg.mri, 'RepetitionTime'                 ); % REQUIRED. The time in seconds between the beginning of an acquisition of one volume and the beginning of acquisition of the volume following it (TR). Please note that this definition includes time between scans (when no data has been acquired) in case of sparse acquisition schemes. This value needs to be consistent with the pixdim[4] field (after accounting for units stored in xyzt_units field) in the NIfTI header. This field is mutually exclusive with VolumeTiming and is derived from DICOM Tag 0018, 0080 and converted to seconds.
cfg.mri.VolumeTiming                  = ft_getopt(cfg.mri, 'VolumeTiming'                   ); % REQUIRED. The time at which each volume was acquired during the acquisition. It is described using a list of times (in JSON format) referring to the onset of each volume in the BOLD series. The list must have the same length as the BOLD series, and the values must be non-negative and monotonically increasing. This field is mutually exclusive with RepetitionTime and DelayTime. If defined, this requires acquisition time (TA) be defined via either SliceTiming or AcquisitionDuration be defined.

%% MEG specific fields
% Manufacturer and ManufacturersModelName are general
cfg.meg.SamplingFrequency             = ft_getopt(cfg.meg, 'SamplingFrequency'           ); % REQUIRED. Sampling frequency (in Hz) of all the data in the recording, regardless of their type (e.g., 2400)
cfg.meg.PowerLineFrequency            = ft_getopt(cfg.meg, 'PowerLineFrequency'          ); % REQUIRED. Frequency (in Hz) of the power grid at the geographical location of the MEG instrument (i.e. 50 or 60)
cfg.meg.DewarPosition                 = ft_getopt(cfg.meg, 'DewarPosition'               ); % REQUIRED. Position of the dewar during the MEG scan: "upright", "supine" or "degrees" of angle from vertical: for example on CTF systems, upright=15??, supine = 90??.
cfg.meg.SoftwareFilters               = ft_getopt(cfg.meg, 'SoftwareFilters'             ); % REQUIRED. List of temporal and/or spatial software filters applied, orideally key:valuepairsofpre-appliedsoftwarefiltersandtheir parameter values: e.g., {"SSS": {"frame": "head", "badlimit": 7}}, {"SpatialCompensation": {"GradientOrder": Order of the gradient compensation}}. Write "n/a" if no software filters applied.
cfg.meg.DigitizedLandmarks            = ft_getopt(cfg.meg, 'DigitizedLandmarks'          ); % REQUIRED. Boolean ("true" or "false") value indicating whether anatomical landmark points (i.e. fiducials) are contained within this recording.
cfg.meg.DigitizedHeadPoints           = ft_getopt(cfg.meg, 'DigitizedHeadPoints'         ); % REQUIRED. Boolean ("true" or "false") value indicating whether head points outlining the scalp/face surface are contained within this recording.
cfg.meg.MEGChannelCount               = ft_getopt(cfg.meg, 'MEGChannelCount'             ); % OPTIONAL. Number of MEG channels (e.g. 275)
cfg.meg.MEGREFChannelCount            = ft_getopt(cfg.meg, 'MEGREFChannelCount'          ); % OPTIONAL. Number of MEG reference channels (e.g. 23). For systems without such channels (e.g. Neuromag Vectorview), MEGREFChannelCount"=0
cfg.meg.EEGChannelCount               = ft_getopt(cfg.meg, 'EEGChannelCount'             ); % OPTIONAL. Number of EEG channels recorded simultaneously (e.g. 21)
cfg.meg.ECOGChannelCount              = ft_getopt(cfg.meg, 'ECOGChannelCount'            ); % OPTIONAL. Number of ECoG channels
cfg.meg.SEEGChannelCount              = ft_getopt(cfg.meg, 'SEEGChannelCount'            ); % OPTIONAL. Number of SEEG channels
cfg.meg.EOGChannelCount               = ft_getopt(cfg.meg, 'EOGChannelCount'             ); % OPTIONAL. Number of EOG channels
cfg.meg.ECGChannelCount               = ft_getopt(cfg.meg, 'ECGChannelCount'             ); % OPTIONAL. Number of ECG channels
cfg.meg.EMGChannelCount               = ft_getopt(cfg.meg, 'EMGChannelCount'             ); % OPTIONAL. Number of EMG channels
cfg.meg.MiscChannelCount              = ft_getopt(cfg.meg, 'MiscChannelCount'            ); % OPTIONAL. Number of miscellaneous analog channels for auxiliary signals
cfg.meg.TriggerChannelCount           = ft_getopt(cfg.meg, 'TriggerChannelCount'         ); % OPTIONAL. Number of channels for digital (TTL bit level) triggers
cfg.meg.RecordingDuration             = ft_getopt(cfg.meg, 'RecordingDuration'           ); % OPTIONAL. Length of the recording in seconds (e.g. 3600)
cfg.meg.RecordingType                 = ft_getopt(cfg.meg, 'RecordingType'               ); % OPTIONAL. Defines whether the recording is "continuous" or "epoched"; this latter limited to time windows about events of interest (e.g., stimulus presentations, subject responses etc.)
cfg.meg.EpochLength                   = ft_getopt(cfg.meg, 'EpochLength'                 ); % OPTIONAL. Duration of individual epochs in seconds (e.g. 1) in case of epoched data
cfg.meg.ContinuousHeadLocalization    = ft_getopt(cfg.meg, 'ContinuousHeadLocalization'  ); % OPTIONAL. Boolean ("true" or "false") value indicating whether continuous head localisation was performed.
cfg.meg.HeadCoilFrequency             = ft_getopt(cfg.meg, 'HeadCoilFrequency'           ); % OPTIONAL. List of frequencies (in Hz) used by the head localisation coils ("HLC" in CTF systems, "HPI" in Neuromag/Elekta, "COH" in 4D/BTi) that track the subject's head position in the MEG helmet (e.g. [293, 307, 314, 321])
cfg.meg.MaxMovement                   = ft_getopt(cfg.meg, 'MaxMovement'                 ); % OPTIONAL. Maximum head movement (in mm) detected during the recording, as measured by the head localisation coils (e.g., 4.8)
cfg.meg.SubjectArtefactDescription    = ft_getopt(cfg.meg, 'SubjectArtefactDescription'  ); % OPTIONAL. Freeform description of the observed subject artefact and its possible cause (e.g. "Vagus Nerve Stimulator", "non-removable implant"). If this field is set to "n/a", it will be interpreted as absence of major source of artifacts except cardiac and blinks.
cfg.meg.AssociatedEmptyRoom           = ft_getopt(cfg.meg, 'AssociatedEmptyRoom'         ); % OPTIONAL. Relative path in BIDS folder structure to empty-room file associated with the subject's MEG recording. The path needs to use forward slashes instead of backward slashes (e.g. "sub-emptyroom/ses-<label>/meg/sub-emptyroom_ses-<label>_ta sk-noise_run-<label>_meg.ds").
cfg.meg.HardwareFilters               = ft_getopt(cfg.meg, 'HardwareFilters'             ); % RECOMMENDED. List of temporal hardware filters applied. Ideally key:value pairs of pre-applied hardware filters and their parameter values: e.g., {"HardwareFilters": {"Highpass RC filter": {"Half amplitude cutoff (Hz)": 0.0159, "Roll-off": "6dB/Octave"}}}. Write n/a if no hardware filters applied.

%% Specific EEG fields - if recorded with the MEG system
cfg.meg.EEGPlacementScheme            = ft_getopt(cfg.meg, 'EEGPlacementScheme'          ); % OPTIONAL. Placement scheme of EEG electrodes. Either the name of a standardised placement system (e.g., "10-20") or a list of standardised electrode names (e.g. ["Cz", "Pz"]).
cfg.meg.CapManufacturer               = ft_getopt(cfg.meg, 'CapManufacturer'             ); % OPTIONAL. Manufacturer of the EEG cap (e.g. EasyCap)
cfg.meg.CapManufacturersModelName     = ft_getopt(cfg.meg, 'CapManufacturersModelName'   ); % OPTIONAL. Manufacturer’s designation of the EEG cap model (e.g., M10)
cfg.meg.EEGReference                  = ft_getopt(cfg.meg, 'EEGReference'                ); % OPTIONAL. Description of the type of EEG reference used (e.g., M1 for left mastoid, average, or longitudinal bipolar).

%% EEG specific fields
cfg.eeg.EEGReference                  = ft_getopt(cfg.eeg, 'EEGReference'                ); % Description of the type of reference used (common", "average", "DRL", "bipolar" ).  Any specific electrode used as reference should be indicated as such in the channels.tsv file
cfg.eeg.SamplingFrequency             = ft_getopt(cfg.eeg, 'SamplingFrequency'           ); % Sampling frequency (in Hz) of the EEG recording (e.g. 2400)
cfg.eeg.PowerLineFrequency            = ft_getopt(cfg.eeg, 'PowerLineFrequency'          ); % Frequency (in Hz) of the power grid where the EEG is installed (i.e. 50 or 60).
cfg.eeg.SoftwareFilters               = ft_getopt(cfg.eeg, 'SoftwareFilters'             ); % List of temporal software filters applied or ideally  key:value pairs of pre-applied filters and their parameter values
cfg.eeg.CapManufacturer               = ft_getopt(cfg.eeg, 'CapManufacturer'             ); % name of the cap manufacturer
cfg.eeg.CapModelName                  = ft_getopt(cfg.eeg, 'CapModelName'                ); % Manufacturer's designation of the EEG cap model (e.g. "CAPML128", "actiCAP 64Ch Standard-2")
% Manufacturer and ManufacturersModelName are general
cfg.eeg.EEGChannelCount               = ft_getopt(cfg.eeg, 'EEGChannelCount'             ); % Number of EEG channels included in the recording (e.g. 128).
cfg.eeg.ECGChannelCount               = ft_getopt(cfg.eeg, 'ECGChannelCount'             ); % Number of ECG channels included in the recording (e.g. 1).
cfg.eeg.EMGChannelCount               = ft_getopt(cfg.eeg, 'EMGChannelCount'             ); % Number of EMG channels included in the recording (e.g. 2).
cfg.eeg.EOGChannelCount               = ft_getopt(cfg.eeg, 'EOGChannelCount'             ); % Number of EOG channels included in the recording (e.g. 2).
cfg.eeg.MiscChannelCount              = ft_getopt(cfg.eeg, 'MiscChannelCount'            ); % Number of miscellaneous analog channels for auxiliary  signals
cfg.eeg.TriggerChannelCount           = ft_getopt(cfg.eeg, 'TriggerChannelCount'         ); % Number of channels for digital and analog triggers.
cfg.eeg.RecordingDuration             = ft_getopt(cfg.eeg, 'RecordingDuration'           ); % Length of the recording in seconds (e.g. 3600)
cfg.eeg.RecordingType                 = ft_getopt(cfg.eeg, 'RecordingType'               ); % "continuous", "epoched"
cfg.eeg.EpochLength                   = ft_getopt(cfg.eeg, 'EpochLength'                 ); % Duration of individual epochs in seconds (e.g. 1). If recording was continuous, set value to Inf or leave out the field.
cfg.eeg.HeadCircumference             = ft_getopt(cfg.eeg, 'HeadCircumference'           ); % RECOMMENDED. Circumference of the participants head, expressed in cm (e.g., 58).
cfg.eeg.EEGPlacementScheme            = ft_getopt(cfg.eeg, 'EEGPlacementScheme'          ); % Placement scheme of the EEG electrodes. Either the name of a placement system (e.g. "10-20", "equidistant", "geodesic") or a list of electrode positions (e.g. "Cz", "Pz").
cfg.eeg.EEGGround                     = ft_getopt(cfg.eeg, 'EEGGround'                   ); % RECOMMENDED. Description of the location of the ground electrode (e.g., "placed on right mastoid (M2)").
cfg.eeg.HardwareFilters               = ft_getopt(cfg.eeg, 'HardwareFilters'             ); % List of hardware (amplifier) filters applied or ideally  key:value pairs of pre-applied filters and their parameter values
cfg.eeg.SubjectArtefactDescription    = ft_getopt(cfg.eeg, 'SubjectArtefactDescription'  ); % Freeform description of the observed subject artefact and its possible cause (e.g. "Vagus Nerve Stimulator", "non-removable implant"). If this field is left empty, it will be interpreted as absence of  a source of (constantly present) artifacts.

%% iEEG specific fields
cfg.ieeg.iEEGReference                   = ft_getopt(cfg, 'iEEGReference'                  ); % REQUIRED. General description of the reference scheme used and (when applicable) of location of the reference electrode in the raw recordings (e.g. "left mastoid”, “bipolar”, “T01” for electrode with name T01, “intracranial electrode on top of a grid, not included with data”, “upside down electrode”). If different channels have a different reference, this field should have a general description and the channel specific reference should be defined in the _channels.tsv file.
cfg.ieeg.SamplingFrequency               = ft_getopt(cfg, 'SamplingFrequency'              ); % REQUIRED. Sampling frequency (in Hz) of all the iEEG channels in the recording (e.g., 2400). All other channels should have frequency specified as well in the channels.tsv file.
cfg.ieeg.PowerLineFrequency              = ft_getopt(cfg, 'PowerLineFrequency'             ); % REQUIRED. Frequency (in Hz) of the power grid where the iEEG recording was done (i.e. 50 or 60)
cfg.ieeg.SoftwareFilters                 = ft_getopt(cfg, 'SoftwareFilters'                ); % REQUIRED. List of temporal software filters applied or ideally  key:value pairs of pre-applied filters and their parameter values. (n/a if none).
% cfg.ieeg.DCOffsetCorrection see https://github.com/bids-standard/bids-specification/issues/237
cfg.ieeg.HardwareFilters                 = ft_getopt(cfg, 'HardwareFilters'                ); % REQUIRED. List of hardware (amplifier) filters applied with  key:value pairs of filter parameters and their values.
cfg.ieeg.ElectrodeManufacturer           = ft_getopt(cfg, 'ElectrodeManufacturer'          ); % RECOMMENDED. can be used if all electrodes are of the same manufacturer (e.g. AD-TECH, DIXI). If electrodes of different manufacturers are used, please use the corresponding table in the _electrodes.tsv file.
cfg.ieeg.ElectrodeManufacturersModelName = ft_getopt(cfg, 'ElectrodeManufacturersModelName'); % RECOMMENDED. If different electrode types are used, please use the corresponding table in the _electrodes.tsv file.
% Manufacturer and ManufacturersModelName are general
cfg.ieeg.ECOGChannelCount                = ft_getopt(cfg, 'ECOGChannelCount'               ); % RECOMMENDED. Number of iEEG surface channels included in the recording (e.g. 120)
cfg.ieeg.SEEGChannelCount                = ft_getopt(cfg, 'SEEGChannelCount'               ); % RECOMMENDED. Number of iEEG depth channels included in the recording (e.g. 8)
cfg.ieeg.EEGChannelCount                 = ft_getopt(cfg, 'EEGChannelCount'                ); % RECOMMENDED. Number of scalp EEG channels recorded simultaneously (e.g. 21)
cfg.ieeg.EOGChannelCount                 = ft_getopt(cfg, 'EOGChannelCount'                ); % RECOMMENDED. Number of EOG channels
cfg.ieeg.ECGChannelCount                 = ft_getopt(cfg, 'ECGChannelCount'                ); % RECOMMENDED. Number of ECG channels
cfg.ieeg.EMGChannelCount                 = ft_getopt(cfg, 'EMGChannelCount'                ); % RECOMMENDED. Number of EMG channels
cfg.ieeg.MiscChannelCount                = ft_getopt(cfg, 'MiscChannelCount'               ); % RECOMMENDED. Number of miscellaneous analog channels for auxiliary  signals
cfg.ieeg.TriggerChannelCount             = ft_getopt(cfg, 'TriggerChannelCount'            ); % RECOMMENDED. Number of channels for digital (TTL bit level) triggers
cfg.ieeg.RecordingDuration               = ft_getopt(cfg, 'RecordingDuration'              ); % RECOMMENDED. Length of the recording in seconds (e.g. 3600)
cfg.ieeg.RecordingType                   = ft_getopt(cfg, 'RecordingType'                  ); % RECOMMENDED. Defines whether the recording is  “continuous” or  “epoched”; this latter limited to time windows about events of interest (e.g., stimulus presentations, subject responses etc.)
cfg.ieeg.EpochLength                     = ft_getopt(cfg, 'EpochLength'                    ); % RECOMMENDED. Duration of individual epochs in seconds (e.g. 1) in case of epoched data
cfg.ieeg.iEEGGround                      = ft_getopt(cfg, 'iEEGGround'                     ); % RECOMMENDED. Description  of the location of the ground electrode (“placed on right mastoid (M2)”).
cfg.ieeg.iEEGPlacementScheme             = ft_getopt(cfg, 'iEEGPlacementScheme'            ); % RECOMMENDED. Freeform description of the placement of the iEEG electrodes. Left/right/bilateral/depth/surface (e.g. “left frontal grid and bilateral hippocampal depth” or “surface strip and STN depth” or “clinical indication bitemporal, bilateral temporal strips and left grid”).
cfg.ieeg.iEEGElectrodeGroups             = ft_getopt(cfg, 'iEEGElectrodeGroups'            ); % RECOMMENDED. Field to describe the way electrodes are grouped into strips, grids or depth probes e.g. {'grid1': "10x8 grid on left temporal pole", 'strip2': "1x8 electrode strip on xxx"}.
cfg.ieeg.SubjectArtefactDescription      = ft_getopt(cfg, 'SubjectArtefactDescription'     ); % RECOMMENDED. Freeform description of the observed subject artefact and its possible cause (e.g. “door open”, ”nurse walked into room at 2 min”, ”seizure at 10 min”). If this field is left empty, it will be interpreted as absence of artifacts.
cfg.ieeg.ElectricalStimulation           = ft_getopt(cfg, 'ElectricalStimulation'          ); % OPTIONAL. Boolean field to specify if electrical stimulation was done during the recording (options are “true” or “false”). Parameters for event-like stimulation should be specified in the _events.tsv file (see example underneath).
cfg.ieeg.ElectricalStimulationParameters = ft_getopt(cfg, 'ElectricalStimulationParameters'); % OPTIONAL. Free form description of stimulation parameters, such as frequency, shape etc. Specific onsets can be specified in the _events.tsv file. Specific shapes can be described here in freeform text.

%% EMG is not part of the official BIDS specification
cfg.emg.SamplingFrequency                 = ft_getopt(cfg.emg, 'SamplingFrequency'                 );
cfg.emg.RecordingDuration                 = ft_getopt(cfg.emg, 'RecordingDuration'                 );
cfg.emg.RecordingType                     = ft_getopt(cfg.emg, 'RecordingType'                     );
cfg.emg.PowerLineFrequency                = ft_getopt(cfg.emg, 'PowerLineFrequency'                );
cfg.emg.HardwareFilters                   = ft_getopt(cfg.emg, 'HardwareFilters'                   );
cfg.emg.SoftwareFilters                   = ft_getopt(cfg.emg, 'SoftwareFilters'                   );
cfg.emg.EMGChannelCount                   = ft_getopt(cfg.emg, 'EMGChannelCount'                   );
cfg.emg.EOGChannelCount                   = ft_getopt(cfg.emg, 'EOGChannelCount'                   );
cfg.emg.ECGChannelCount                   = ft_getopt(cfg.emg, 'ECGChannelCount'                   );
cfg.emg.ElectrodeManufacturer             = ft_getopt(cfg.emg, 'ElectrodeManufacturer'             );
cfg.emg.ElectrodeManufacturersModelName   = ft_getopt(cfg.emg, 'ElectrodeManufacturersModelName'   );
cfg.emg.EMGPlacementScheme                = ft_getopt(cfg.emg, 'EMGPlacementScheme'                );
cfg.emg.EMGReference                      = ft_getopt(cfg.emg, 'EMGReference'                      );
cfg.emg.EMGGround                         = ft_getopt(cfg.emg, 'EMGGround'                         );

%% EXG is not part of the official BIDS specification
cfg.exg.SamplingFrequency                 = ft_getopt(cfg.exg, 'SamplingFrequency'                 );
cfg.exg.RecordingDuration                 = ft_getopt(cfg.exg, 'RecordingDuration'                 );
cfg.exg.RecordingType                     = ft_getopt(cfg.exg, 'RecordingType'                     );

%% audio is not part of the official BIDS specification
cfg.audio.SampleRate                      = ft_getopt(cfg.audio, 'SampleRate'     );
cfg.audio.ChannelCount                    = ft_getopt(cfg.audio, 'ChannelCount'   );

%% video is not part of the official BIDS specification
cfg.video.FrameRate                       = ft_getopt(cfg.video, 'FrameRate'           );
cfg.video.Width                           = ft_getopt(cfg.video, 'Width'               );
cfg.video.Height                          = ft_getopt(cfg.video, 'Height'              );
cfg.video.BitsPerPixel                    = ft_getopt(cfg.video, 'BitsPerPixel'        );
cfg.video.AudioSampleRate                 = ft_getopt(cfg.video, 'AudioSampleRate'     );
cfg.video.AudioChannelCount               = ft_getopt(cfg.video, 'AudioChannelCount'   );

%% physio is part of the official BIDS specification and goes into 'beh' or in 'func'
cfg.physio.Columns                        = ft_getopt(cfg.physio, 'Columns'              );
cfg.physio.StartTime                      = ft_getopt(cfg.physio, 'StartTime'            );
cfg.physio.SamplingFrequency              = ft_getopt(cfg.physio, 'SamplingFrequency'    );

%% stim is part of the official BIDS specification and goes into 'beh' or in 'func'
cfg.stim.Columns                          = ft_getopt(cfg.stim, 'Columns'              );
cfg.stim.StartTime                        = ft_getopt(cfg.stim, 'StartTime'            );
cfg.stim.SamplingFrequency                = ft_getopt(cfg.stim, 'SamplingFrequency'    );

%% eyetracker is not part of the official BIDS specification
% this follows https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/06-physiological-and-other-continuous-recordings.html
cfg.eyetracker.Columns                    = ft_getopt(cfg.eyetracker, 'Columns'              );
cfg.eyetracker.StartTime                  = ft_getopt(cfg.eyetracker, 'StartTime'            );
cfg.eyetracker.SamplingFrequency          = ft_getopt(cfg.eyetracker, 'SamplingFrequency'    );

%% motion is not part of the official BIDS specification
% this follows https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/06-physiological-and-other-continuous-recordings.html
cfg.motion.Columns                        = ft_getopt(cfg.motion, 'Columns'              );
cfg.motion.StartTime                      = ft_getopt(cfg.motion, 'StartTime'            );
cfg.motion.SamplingFrequency              = ft_getopt(cfg.motion, 'SamplingFrequency'    );

%% information for the coordsystem.json file for MEG, EEG and iEEG
cfg.coordsystem.MEGCoordinateSystem                             = ft_getopt(cfg.coordsystem, 'MEGCoordinateSystem'                            ); % REQUIRED. Defines the coordinate system for the MEG sensors. See Appendix VIII: preferred names of Coordinate systems. If "Other", provide definition of the coordinate system in [MEGCoordinateSystemDescription].
cfg.coordsystem.MEGCoordinateUnits                              = ft_getopt(cfg.coordsystem, 'MEGCoordinateUnits'                             ); % REQUIRED. Units of the coordinates of MEGCoordinateSystem. MUST be ???m???, ???cm???, or ???mm???.
cfg.coordsystem.MEGCoordinateSystemDescription                  = ft_getopt(cfg.coordsystem, 'MEGCoordinateSystemDescription'                 ); % OPTIONAL. Freeform text description or link to document describing the MEG coordinate system system in detail.
cfg.coordsystem.HeadCoilCoordinates                             = ft_getopt(cfg.coordsystem, 'HeadCoilCoordinateSystem'                       ); % OPTIONAL. Key:value pairs describing head localization coil labels and their coordinates, interpreted following the HeadCoilCoordinateSystem, e.g., {"NAS": [12.7,21.3,13.9], "LPA": [5.2,11.3,9.6], "RPA": [20.2,11.3,9.1]}. Note that coils are not always placed at locations that have a known anatomical name (e.g. for Neuromag/Elekta, Yokogawa systems); in that case generic labels can be used (e.g. {"coil1": [122,213,123], "coil2": [67,123,86], "coil3": [219,110,81]} ).
cfg.coordsystem.HeadCoilCoordinateSystem                        = ft_getopt(cfg.coordsystem, 'HeadCoilCoordinates'                            ); % OPTIONAL. Defines the coordinate system for the coils. See Appendix VIII: preferred names of Coordinate systems. If "Other", provide definition of the coordinate system in HeadCoilCoordinateSystemDescription.
cfg.coordsystem.HeadCoilCoordinateUnits                         = ft_getopt(cfg.coordsystem, 'HeadCoilCoordinateUnits'                        ); % OPTIONAL. Units of the coordinates of HeadCoilCoordinateSystem. MUST be ???m???, ???cm???, or ???mm???.
cfg.coordsystem.HeadCoilCoordinateSystemDescription             = ft_getopt(cfg.coordsystem, 'HeadCoilCoordinateSystemDescription'            ); % OPTIONAL. Freeform text description or link to document describing the Head Coil coordinate system system in detail.
cfg.coordsystem.DigitizedHeadPoints                             = ft_getopt(cfg.coordsystem, 'DigitizedHeadPoints'                            ); % OPTIONAL. Relative path to the file containing the locations of digitized head points collected during the session (e.g., "sub-01_headshape.pos"). RECOMMENDED for all MEG systems, especially for CTF and 4D/BTi. For Neuromag/Elekta the head points will be stored in the fif file.
cfg.coordsystem.DigitizedHeadPointsCoordinateSystem             = ft_getopt(cfg.coordsystem, 'DigitizedHeadPointsCoordinateSystem'            ); % OPTIONAL. Defines the coordinate system for the digitized head points. See Appendix VIII: preferred names of Coordinate systems. If "Other", provide definition of the coordinate system in DigitizedHeadPointsCoordinateSystemDescription.
cfg.coordsystem.DigitizedHeadPointsCoordinateUnits              = ft_getopt(cfg.coordsystem, 'DigitizedHeadPointsCoordinateUnits'             ); % OPTIONAL. Units of the coordinates of DigitizedHeadPointsCoordinateSystem. MUST be ???m???, ???cm???, or ???mm???.
cfg.coordsystem.DigitizedHeadPointsCoordinateSystemDescription  = ft_getopt(cfg.coordsystem, 'DigitizedHeadPointsCoordinateSystemDescription' ); % OPTIONAL. Freeform text description or link to document describing the Digitized head Points coordinate system system in detail.
cfg.coordsystem.EEGCoordinateSystem                             = ft_getopt(cfg.coordsystem, 'EEGCoordinateSystem'                            ); % OPTIONAL. Describes how the coordinates of the EEG sensors are to be interpreted.
cfg.coordsystem.EEGCoordinateUnits                              = ft_getopt(cfg.coordsystem, 'EEGCoordinateUnits'                             ); % OPTIONAL. Units of the coordinates of EEGCoordinateSystem. MUST be ???m???, ???cm???, or ???mm???.
cfg.coordsystem.EEGCoordinateSystemDescription                  = ft_getopt(cfg.coordsystem, 'EEGCoordinateSystemDescription'                 ); % OPTIONAL. Freeform text description or link to document describing the EEG coordinate system system in detail.
cfg.coordsystem.iEEGCoordinateSystem                            = ft_getopt(cfg.coordsystem, 'iEEGCoordinateSystem'                           ); % REQUIRED. Defines the coordinate system for the iEEG electrodes. See Appendix VIII for a list of restricted keywords. If positions correspond to pixel indices in a 2D image (of either a volume-rendering, surface-rendering, operative photo, or operative drawing), this must be "Pixels". For more information, see the section on 2D coordinate systems
cfg.coordsystem.iEEGCoordinateUnits	                            = ft_getopt(cfg.coordsystem, 'iEEGCoordinateUnits'                            ); % REQUIRED. Units of the _electrodes.tsv, MUST be "m", "mm", "cm" or "pixels".
cfg.coordsystem.iEEGCoordinateSystemDescription	                = ft_getopt(cfg.coordsystem, 'iEEGCoordinateSystemDescription'                ); % RECOMMENDED. Freeform text description or link to document describing the iEEG coordinate system system in detail (e.g., "Coordinate system with the origin at anterior commissure (AC), negative y-axis going through the posterior commissure (PC), z-axis going to a mid-hemisperic point which lies superior to the AC-PC line, x-axis going to the right").
cfg.coordsystem.iEEGCoordinateProcessingDescription             = ft_getopt(cfg.coordsystem, 'iEEGCoordinateProcessingDescription'            ); % RECOMMENDED. Has any post-processing (such as projection) been done on the electrode positions (e.g., "surface_projection", "none").
cfg.coordsystem.iEEGCoordinateProcessingReference	              = ft_getopt(cfg.coordsystem, 'iEEGCoordinateProcessingReference'              ); % RECOMMENDED. A reference to a paper that defines in more detail the method used to localize the electrodes and to post-process the electrode positions. .
cfg.coordsystem.IntendedFor                                     = ft_getopt(cfg.coordsystem, 'IntendedFor'                                    ); % OPTIONAL. Path or list of path relative to the subject subfolder pointing to the structural MRI, possibly of different types if a list is specified, to be used with the MEG recording. The path(s) need(s) to use forward slashes instead of backward slashes (e.g. "ses-<label>/anat/sub-01_T1w.nii.gz").
cfg.coordsystem.AnatomicalLandmarkCoordinates                   = ft_getopt(cfg.coordsystem, 'AnatomicalLandmarkCoordinates'                  ); % OPTIONAL. Key:value pairs of the labels and 3-D digitized locations of anatomical landmarks, interpreted following the AnatomicalLandmarkCoordinateSystem, e.g., {"NAS": [12.7,21.3,13.9], "LPA": [5.2,11.3,9.6], "RPA": [20.2,11.3,9.1]}.
cfg.coordsystem.AnatomicalLandmarkCoordinateSystem              = ft_getopt(cfg.coordsystem, 'AnatomicalLandmarkCoordinateSystem'             ); % OPTIONAL. Defines the coordinate system for the anatomical landmarks. See Appendix VIII: preferred names of Coordinate systems. If "Other", provide definition of the coordinate system in AnatomicalLandmarkCoordinateSystemDescripti on.
cfg.coordsystem.AnatomicalLandmarkCoordinateUnits               = ft_getopt(cfg.coordsystem, 'AnatomicalLandmarkCoordinateUnits'              ); % OPTIONAL. Units of the coordinates of AnatomicalLandmarkCoordinateSystem. MUST be ???m???, ???cm???, or ???mm???.
cfg.coordsystem.AnatomicalLandmarkCoordinateSystemDescription   = ft_getopt(cfg.coordsystem, 'AnatomicalLandmarkCoordinateSystemDescription'  ); % OPTIONAL. Freeform text description or link to document describing the Head Coil coordinate system system in detail.
cfg.coordsystem.FiducialsDescription                            = ft_getopt(cfg.coordsystem, 'FiducialsDescription'                           ); % OPTIONAL. A freeform text field documenting the anatomical landmarks that were used and how the head localization coils were placed relative to these. This field can describe, for instance, whether the true anatomical locations of the left and right pre-auricular points were used and digitized, or rather whether they were defined as the intersection between the tragus and the helix (the entry of the ear canal), or any other anatomical description of selected points in the vicinity of the ears.

%% columns in the channels.tsv
cfg.channels.name               = ft_getopt(cfg.channels, 'name'               , nan);  % REQUIRED. Channel name (e.g., MRT012, MEG023)
cfg.channels.type               = ft_getopt(cfg.channels, 'type'               , nan);  % REQUIRED. Type of channel; MUST use the channel types listed below.
cfg.channels.units              = ft_getopt(cfg.channels, 'units'              , nan);  % REQUIRED. Physical unit of the data values recorded by this channel in SI (see Appendix V: Units for allowed symbols).
cfg.channels.sampling_frequency = ft_getopt(cfg.channels, 'sampling_frequency' , nan);  % OPTIONAL. Sampling rate of the channel in Hz.
cfg.channels.description        = ft_getopt(cfg.channels, 'description'        , nan);  % OPTIONAL. Brief free-text description of the channel, or other information of interest. See examples below.
cfg.channels.low_cutoff         = ft_getopt(cfg.channels, 'low_cutoff'         , nan);  % OPTIONAL. Frequencies used for the high-pass filter applied to the channel in Hz. If no high-pass filter applied, use n/a.
cfg.channels.high_cutoff        = ft_getopt(cfg.channels, 'high_cutoff'        , nan);  % OPTIONAL. Frequencies used for the low-pass filter applied to the channel in Hz. If no low-pass filter applied, use n/a. Note that hardware anti-aliasing in A/D conversion of all MEG/EEG electronics applies a low-pass filter; specify its frequency here if applicable.
cfg.channels.notch              = ft_getopt(cfg.channels, 'notch'              , nan);  % OPTIONAL. Frequencies used for the notch filter applied to the channel, in Hz. If no notch filter applied, use n/a.
cfg.channels.software_filters   = ft_getopt(cfg.channels, 'software_filters'   , nan);  % OPTIONAL. List of temporal and/or spatial software filters applied (e.g. "SSS", "SpatialCompensation"). Note that parameters should be defined in the general MEG sidecar .json file. Indicate n/a in the absence of software filters applied.
cfg.channels.status             = ft_getopt(cfg.channels, 'status'             , nan);  % OPTIONAL. Data quality observed on the channel (good/bad). A channel is considered bad if its data quality is compromised by excessive noise. Description of noise type SHOULD be provided in [status_description].
cfg.channels.status_description = ft_getopt(cfg.channels, 'status_description' , nan);  % OPTIONAL. Freeform text description of noise or artifact affecting data quality on the channel. It is meant to explain why the channel was declared bad in [status].

%% columns in the electrodes.tsv
cfg.electrodes.name             = ft_getopt(cfg.channels, 'name'               , nan);  % REQUIRED. Name of the electrode
cfg.electrodes.x                = ft_getopt(cfg.channels, 'x'                  , nan);  % REQUIRED. Recorded position along the x-axis
cfg.electrodes.y                = ft_getopt(cfg.channels, 'y'                  , nan);  % REQUIRED. Recorded position along the y-axis
cfg.electrodes.z                = ft_getopt(cfg.channels, 'z'                  , nan);  % REQUIRED. Recorded position along the z-axis
cfg.electrodes.type             = ft_getopt(cfg.channels, 'type'               , nan);  % RECOMMENDED. Type of the electrode (e.g., cup, ring, clip-on, wire, needle)
cfg.electrodes.material         = ft_getopt(cfg.channels, 'material'           , nan);  % RECOMMENDED. Material of the electrode, e.g., Tin, Ag/AgCl, Gold
cfg.electrodes.impedance        = ft_getopt(cfg.channels, 'impedance'          , nan);  % RECOMMENDED. Impedance of the electrode in kOhm

%% information for the participants.tsv
cfg.participants = ft_getopt(cfg, 'participants', struct());

%% information for the scans.tsv
cfg.scans = ft_getopt(cfg, 'scans', struct());

%% sanity checks and determine the default method/outputfile

% the task is both part of the file name (cfg.task) and is also one of the general JSON metadata fields (cfg.TaskName)
if isempty(cfg.task) && isempty(cfg.TaskName)
  % this is fine
elseif ~isempty(cfg.task) && isempty(cfg.TaskName)
  cfg.TaskName = cfg.task;
elseif isempty(cfg.task) && ~isempty(cfg.TaskName)
  cfg.task = cfg.TaskName;
elseif ~isempty(cfg.task) && ~isempty(cfg.TaskName)
  if ~strcmp(cfg.task, cfg.TaskName)
    ft_error('cfg.task and cfg.TaskName should be identical');
  end
end

% construct the output directory and file name
if isempty(cfg.outputfile)
  if isempty(cfg.method) && isempty(cfg.bidsroot) && ~isempty(cfg.dataset)
    cfg.outputfile = cfg.dataset;
    cfg.method = 'decorate';
    ft_notice('using cfg.outputfile=''%s'' and cfg.method=''%s''', cfg.outputfile, cfg.method);
  elseif isempty(cfg.bidsroot)
    ft_error('cfg.bidsroot is required to construct BIDS output directory and file');
  elseif isempty(cfg.sub)
    ft_error('cfg.sub is required to construct BIDS output directory and file');
  elseif isempty(cfg.datatype)
    ft_error('cfg.datatype is required to construct BIDS output directory and file');
  else
    dirname = datatype2dirname(cfg.datatype);
    filename = ['sub-' cfg.sub];
    filename = add_entity(filename, 'ses',  cfg.ses);
    filename = add_entity(filename, 'task', cfg.task);
    filename = add_entity(filename, 'acq',  cfg.acq);
    filename = add_entity(filename, 'ce',   cfg.ce);
    filename = add_entity(filename, 'rec',  cfg.rec);
    filename = add_entity(filename, 'dir',  cfg.dir);
    filename = add_entity(filename, 'run',  cfg.run);
    filename = add_entity(filename, 'mod',  cfg.mod);
    filename = add_entity(filename, 'echo', cfg.echo);
    filename = add_entity(filename, 'proc', cfg.proc);
    filename = add_datatype(filename, cfg.datatype);
    if ~isempty(cfg.ses)
      % construct the output filename, with session directory
      cfg.outputfile = fullfile(cfg.bidsroot, ['sub-' cfg.sub], ['ses-' cfg.ses], dirname, filename);
    else
      % construct the output filename, without session directory
      cfg.outputfile = fullfile(cfg.bidsroot, ['sub-' cfg.sub], dirname, filename);
    end
    if strcmp(cfg.method, 'copy') && ~isempty(cfg.dataset)
      % copy the file extension from the input dataset
      [~, ~, x] = fileparts(cfg.dataset);
      cfg.outputfile = [cfg.outputfile x];
    end
  end
end

% set the default method
if isempty(cfg.method)
  if ~isempty(cfg.dataset) && ~isequal(cfg.dataset, cfg.outputfile)
    cfg.method = 'convert';
  elseif isempty(cfg.dataset) && ~isempty(varargin)
    cfg.method = 'convert';
  else
    cfg.method = 'decorate';
  end
  ft_notice('using cfg.method=''%s''', cfg.method);
end

% do some sanity checks on the input and the method
if istrue(cfg.deface) && ~strcmp(cfg.method, 'convert')
  ft_error('defacing only works in combination with cfg.method=''convert''');
elseif ft_nargin>1 && ~strcmp(cfg.method, 'convert')
  ft_error('input data only works in combination with cfg.method=''convert''');
end

% do some more sanity checks on the input and the method
switch cfg.method
  case 'decorate'
    if ~isempty(cfg.outputfile) && ~isempty(cfg.dataset) && ~isequal(cfg.dataset, cfg.outputfile)
      ft_error('cfg.dataset and cfg.outputfile should be the same');
    end
  case 'convert'
    if ~isempty(cfg.outputfile) && isequal(cfg.dataset, cfg.outputfile)
      ft_error('cfg.dataset and cfg.outputfile should not be the same');
    end
  case 'copy'
    if ~isempty(cfg.outputfile) && isequal(cfg.dataset, cfg.outputfile)
      ft_error('cfg.dataset and cfg.outputfile should not be the same');
    end
  otherwise
    ft_error('unsupported value for cfg.method')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read the information from the dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert dataset to headerfile and datafile
if nargin>1
  % input data was specified
  varargin{1} = ft_checkdata(varargin{1}, 'datatype', {'raw', 'volume'});
  typ = ft_datatype(varargin{1});
elseif ~isempty(cfg.dataset)
  % data should be read from disk
  cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');
  if isfield(cfg, 'headerfile')
    typ = ft_filetype(cfg.headerfile);
  else
    typ = ft_filetype(cfg.presentationfile);
  end
else
  typ = 'unknown';
end

% determine the json files that are required
need_mri_json           = false;
need_meg_json           = false;
need_eeg_json           = false;
need_ieeg_json          = false;
need_emg_json           = false;
need_exg_json           = false;
need_audio_json         = false;
need_video_json         = false;
need_physio_json        = false;
need_stim_json          = false;
need_eyetracker_json    = false;
need_motion_json        = false;
need_coordsystem_json   = false;
% determine the tsv files that are required
need_events_tsv         = false; % for functional and behavioral experiments
need_channels_tsv       = false; % only needed for MEG/EEG/iEEG/EMG
need_electrodes_tsv     = false; % only needed when actually present as cfg.electrodes, data.elec or as cfg.elec

switch typ
  case {'nifti', 'nifti2', 'nifti_fsl'}
    mri = ft_read_mri(cfg.dataset);
    if ~isempty(cfg.dicomfile)
      % read the header details from the matching DICOM file specified by the user
      dcm = dicominfo(cfg.dicomfile);
    else
      dcm = [];
    end
    need_mri_json = true;

  case 'dicom'
    mri = ft_read_mri(cfg.dataset);
    dcm = dicominfo(cfg.dataset);
    need_mri_json = true;

  case 'volume'
    % the input data structure represents imaging data
    mri = varargin{1};
    if ~isempty(cfg.dicomfile)
      % read the header details from the dicom matching file that was specified by the user
      dcm = dicominfo(cfg.dicomfile);
    elseif isfield(mri, 'hdr') && numel(mri.hdr)>1
      % it looks like an MRI read in using FT_READ_MRI using the FreeSurfer code
      % take the DICOM details from the first slice
      dcm = mri.hdr(1);
    else
      dcm = [];
    end
    need_mri_json = true;

  case {'ctf_ds', 'ctf_meg4', 'ctf_res4', 'ctf151', 'ctf275', 'neuromag_fif', 'neuromag122', 'neuromag306'}
    % it is MEG data from disk and in a supported format
    hdr = ft_read_header(cfg.headerfile, 'checkmaxfilter', false);
    if strcmp(cfg.method, 'convert')
      % the data should be converted and written to disk
      dat = ft_read_data(cfg.datafile, 'header', hdr, 'checkboundary', false, 'begsample', 1, 'endsample', hdr.nSamples*hdr.nTrials);
    end
    % read the triggers from disk
    trigger = ft_read_event(cfg.datafile, 'header', hdr);
    need_meg_json = true;

  case {'brainvision_vhdr', 'edf', 'eeglab_set'}
    % the file on disk contains ExG data in a BIDS compiant format
    hdr = ft_read_header(cfg.headerfile);
    if strcmp(cfg.method, 'convert')
      % the data should be converted and written to disk
      dat = ft_read_data(cfg.datafile, 'header', hdr, 'checkboundary', false, 'begsample', 1, 'endsample', hdr.nSamples*hdr.nTrials);
    end
    % read the triggers from disk
    trigger = ft_read_event(cfg.datafile, 'header', hdr);
    if isequal(cfg.datatype, 'eeg')
      need_eeg_json = true;
    elseif isequal(cfg.datatype, 'ieeg')
      need_ieeg_json = true;
    elseif isequal(cfg.datatype, 'emg')
      need_emg_json = true;
    elseif isequal(cfg.datatype, 'exg')
      need_exg_json = true;
    else
      ft_warning('assuming that the dataset represents EEG');
      need_eeg_json = true;
    end

  case 'presentation_log'
    trigger = ft_read_event(cfg.dataset);
    need_events_tsv = true;

  case {'audio_wav', 'audio_ogg', 'audio_flac', 'audio_au', 'audio_aiff', 'audio_aif', 'audio_aifc', 'audio_mp3', 'audio_m4a', 'audio_mp4'}
    % the file on disk contains audio
    need_audio_json = true;
    try
      audio = audioinfo(cfg.dataset);
    catch
      ft_warning('audio format is unsupported on this MATLAB version and/or operating system');
      audio = struct('SampleRate', nan, 'Duration', nan, 'NumChannels', nan);
    end

  case 'video'
    % the file on disk contains video
    need_video_json = true;
    try
      video = VideoReader(cfg.dataset);
    catch
      ft_warning('video format is unsupported on this MATLAB version and/or operating system');
      video = struct('FrameRate', nan, 'Width', nan, 'Height', nan, 'Duration', nan);
    end
    try
      audio = audioinfo(cfg.dataset);
    catch
      ft_warning('audio format is unsupported on this MATLAB version and/or operating system');
      audio = struct('SampleRate', nan, 'Duration', nan, 'NumChannels', nan);
    end

  case 'raw'
    % the input data structure contains raw timeseries data
    if isequal(cfg.datatype, 'meg')
      need_meg_json = true;
    elseif isequal(cfg.datatype, 'eeg')
      need_eeg_json = true;
    elseif isequal(cfg.datatype, 'ieeg')
      need_ieeg_json = true;
    elseif isequal(cfg.datatype, 'emg')
      need_emg_json = true;
    elseif isequal(cfg.datatype, 'exg')
      need_exg_json = true;
    elseif isequal(cfg.datatype, 'physio')
      need_physio_json = true;
    elseif isequal(cfg.datatype, 'stim')
      need_stim_json = true;
    elseif isequal(cfg.datatype, 'eyetracker')
      need_eyetracker_json = true;
    elseif isequal(cfg.datatype, 'motion')
      need_motion_json = true;
    else
      ft_error('cannot determine the type of the data, please specify cfg.dataype');
    end

    hdr = ft_fetch_header(varargin{1});
    if strcmp(cfg.method, 'convert')
      % the data should be written to disk
      dat     = ft_fetch_data(varargin{1}, 'checkboundary', false, 'begsample', 1, 'endsample', hdr.nSamples*hdr.nTrials);
      trigger = ft_fetch_event(varargin{1});
    end
    try
      % try to get the electrode definition, either from the data or from the configuration
      tmpcfg = keepfields(cfg, {'elec'});
      tmpcfg.senstype = 'eeg';
      elec = ft_fetch_sens(tmpcfg, varargin{1});
      need_electrodes_tsv = true;
    catch
      need_electrodes_tsv = false;
    end

    if ft_senstype(varargin{1}, 'ctf') || ft_senstype(varargin{1}, 'neuromag')
      % use the subsequent MEG-specific metadata handling for the JSON and TSV sidecar files
      typ = ft_senstype(varargin{1});
    end

  otherwise
    % assume that the file on disk contains raw timeseries data that can be read by FieldTrip
    if isequal(cfg.datatype, 'meg')
      need_meg_json = true;
    elseif isequal(cfg.datatype, 'eeg')
      need_eeg_json = true;
    elseif isequal(cfg.datatype, 'ieeg')
      need_ieeg_json = true;
    elseif isequal(cfg.datatype, 'emg')
      need_emg_json = true;
    elseif isequal(cfg.datatype, 'exg')
      need_exg_json = true;
    elseif isequal(cfg.datatype, 'physio')
      need_physio_json = true;
    elseif isequal(cfg.datatype, 'stim')
      need_stim_json = true;
    elseif isequal(cfg.datatype, 'eyetracker')
      need_eyetracker_json = true;
    elseif isequal(cfg.datatype, 'motion')
      need_motion_json = true;
    elseif isequal(cfg.datatype, 'events')
      need_events_tsv = true;
    else
      ft_error('cannot determine the type of the data, please specify cfg.dataype');
    end

    if ~isempty(cfg.dataset)
      hdr = ft_read_header(cfg.headerfile, 'checkmaxfilter', false);
      if strcmp(cfg.method, 'convert')
        % the data should be converted and written to disk
        dat = ft_read_data(cfg.datafile, 'header', hdr, 'checkboundary', false, 'begsample', 1, 'endsample', hdr.nSamples*hdr.nTrials);
        trigger = ft_read_event(cfg.datafile, 'header', hdr);
      end
      % FIXME try to get the electrode definition, either from the data or from the configuration
    end

end % switch typ

if need_meg_json || need_eeg_json || need_ieeg_json
  try
    % try to get the electrode definition, either from data.elec or from cfg.elec
    tmpcfg = keepfields(cfg, {'elec'});
    tmpcfg.senstype = 'eeg';
    if ~isempty(varargin)
      elec = ft_fetch_sens(tmpcfg, varargin{1});
    else
      elec = ft_fetch_sens(tmpcfg);
    end
    need_electrodes_tsv = true;
  catch
    % electrodes can also be specified as cfg.electrodes
    need_electrodes_tsv = ~isnan(cfg.electrodes.name);
  end
end

need_events_tsv       = need_events_tsv       || need_meg_json || need_eeg_json || need_ieeg_json || need_emg_json || need_exg_json || need_eyetracker_json || need_motion_json || (contains(cfg.outputfile, 'task') || ~isempty(cfg.TaskName) || ~isempty(cfg.task));
need_channels_tsv     = need_channels_tsv     || need_meg_json || need_eeg_json || need_ieeg_json || need_emg_json || need_exg_json;
need_coordsystem_json = need_coordsystem_json || need_meg_json || need_electrodes_tsv;

if need_emg_json
  ft_warning('EMG data is not yet part of the official BIDS specification');
  cfg.dataset_description.BIDSVersion = 'n/a';
elseif need_exg_json
  ft_warning('EXG data is not yet part of the official BIDS specification');
  cfg.dataset_description.BIDSVersion = 'n/a';
elseif need_audio_json
  ft_warning('audio data is not yet part of the official BIDS specification');
  cfg.dataset_description.BIDSVersion = 'n/a';
elseif need_video_json
  ft_warning('video data is not yet part of the official BIDS specification');
  cfg.dataset_description.BIDSVersion = 'n/a';
elseif need_eyetracker_json
  ft_warning('eyetracker data is not yet part of the official BIDS specification');
  cfg.dataset_description.BIDSVersion = 'n/a';
elseif need_motion_json
  ft_warning('motion data is not yet part of the official BIDS specification');
  cfg.dataset_description.BIDSVersion = 'n/a';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get the defaults and user-specified settings for each possible sidecar file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make the relevant selection, all json fields start with a capital letter
fn = fieldnames(cfg.dataset_description);
fn = fn(~cellfun(@isempty, regexp(fn, '^[A-Z].*')));
dataset_description_settings = keepfields(cfg.dataset_description, fn);

% make the relevant selection, all json fields start with a capital letter
fn = fieldnames(cfg);
fn = fn(~cellfun(@isempty, regexp(fn, '^[A-Z].*')));
generic_settings = keepfields(cfg, fn);

% make the relevant selection, all json fields start with a capital letter
fn = fieldnames(cfg.mri);
fn = fn(~cellfun(@isempty, regexp(fn, '^[A-Z].*')));
mri_settings = keepfields(cfg.mri, fn);

% make the relevant selection, all json fields start with a capital letter
fn = fieldnames(cfg.meg);
fn = fn(~cellfun(@isempty, regexp(fn, '^[A-Z].*')));
meg_settings = keepfields(cfg.meg, fn);

% make the relevant selection, all json fields start with a capital letter
fn = fieldnames(cfg.eeg);
fn = fn(~cellfun(@isempty, regexp(fn, '^[A-Z].*')));
eeg_settings = keepfields(cfg.eeg, fn);

% make the relevant selection, most json fields start with a capital letter, but some start with iEEG
fn = fieldnames(cfg.ieeg);
fn = fn(~cellfun(@isempty, regexp(fn, '^[A-Z].*|^iEEG')));
ieeg_settings = keepfields(cfg.ieeg, fn);

% make the relevant selection, all json fields start with a capital letter
fn = fieldnames(cfg.emg);
fn = fn(~cellfun(@isempty, regexp(fn, '^[A-Z].*')));
emg_settings = keepfields(cfg.emg, fn);

% make the relevant selection, all json fields start with a capital letter
fn = fieldnames(cfg.exg);
fn = fn(~cellfun(@isempty, regexp(fn, '^[A-Z].*')));
exg_settings = keepfields(cfg.exg, fn);

% make the relevant selection, all json fields start with a capital letter
fn = fieldnames(cfg.audio);
fn = fn(~cellfun(@isempty, regexp(fn, '^[A-Z].*')));
audio_settings = keepfields(cfg.audio, fn);

% make the relevant selection, all json fields start with a capital letter
fn = fieldnames(cfg.video);
fn = fn(~cellfun(@isempty, regexp(fn, '^[A-Z].*')));
video_settings = keepfields(cfg.video, fn);

% make the relevant selection, all json fields start with a capital letter
fn = fieldnames(cfg.physio);
fn = fn(~cellfun(@isempty, regexp(fn, '^[A-Z].*')));
physio_settings = keepfields(cfg.physio, fn);

% make the relevant selection, all json fields start with a capital letter
fn = fieldnames(cfg.stim);
fn = fn(~cellfun(@isempty, regexp(fn, '^[A-Z].*')));
stim_settings = keepfields(cfg.stim, fn);

% make the relevant selection, all json fields start with a capital letter
fn = fieldnames(cfg.eyetracker);
fn = fn(~cellfun(@isempty, regexp(fn, '^[A-Z].*')));
eyetracker_settings = keepfields(cfg.eyetracker, fn);

% make the relevant selection, all json fields start with a capital letter
fn = fieldnames(cfg.motion);
fn = fn(~cellfun(@isempty, regexp(fn, '^[A-Z].*')));
motion_settings = keepfields(cfg.motion, fn);

% make the relevant selection, most json fields start with a capital letter, but some start with iEEG
fn = fieldnames(cfg.coordsystem);
fn = fn(~cellfun(@isempty, regexp(fn, '^[A-Z].*|^iEEG')));
coordsystem_settings = keepfields(cfg.coordsystem, fn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% construct the content for the json and tsv files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% need_mri_json
if need_mri_json
  % start with the information from the DICOM header
  mri_json = keepfields(dcm, [fieldnames(mri_settings); fieldnames(generic_settings)]);
  % merge the information specified by the user with that from the data
  % in case fields appear in both, the first input overrules the second
  mri_json = mergeconfig(mri_settings,     mri_json, false);
  mri_json = mergeconfig(generic_settings, mri_json, false);
end % if need_mri_json

%% need_meg_json
if need_meg_json
  meg_json.SamplingFrequency          = hdr.Fs;
  meg_json.MEGChannelCount            = sum(strcmp(hdr.chantype, 'megmag') | strcmp(hdr.chantype, 'meggrad') | strcmp(hdr.chantype, 'megplanar') | strcmp(hdr.chantype, 'megaxial'));
  meg_json.MEGREFChannelCount         = sum(strcmp(hdr.chantype, 'refmag') | strcmp(hdr.chantype, 'refgrad') | strcmp(hdr.chantype, 'refplanar') | strcmp(hdr.chantype, 'ref'));
  meg_json.EEGChannelCount            = sum(strcmp(hdr.chantype, 'eeg'));
  meg_json.ECOGChannelCount           = sum(strcmp(hdr.chantype, 'ecog'));
  meg_json.SEEGChannelCount           = sum(strcmp(hdr.chantype, 'seeg'));
  meg_json.EOGChannelCount            = sum(strcmp(hdr.chantype, 'eog'));
  meg_json.ECGChannelCount            = sum(strcmp(hdr.chantype, 'ecg'));
  meg_json.EMGChannelCount            = sum(strcmp(hdr.chantype, 'emg'));
  meg_json.MiscChannelCount           = sum(strcmp(hdr.chantype, 'misc') | strcmp(hdr.chantype, 'unknown'));
  meg_json.TriggerChannelCount        = sum(strcmp(hdr.chantype, 'trigger'));
  meg_json.RecordingDuration          = (hdr.nTrials*hdr.nSamples)/hdr.Fs;
  meg_json.EpochLength                = hdr.nSamples/hdr.Fs;
  if ft_senstype(hdr.grad, 'ctf151')
    meg_json.ContinuousHeadLocalization = any(strcmp(hdr.chantype, 'headloc')); % CTF specific
    meg_json.Manufacturer             = 'CTF';
    meg_json.ManufacturersModelName   = 'CTF-151';
  elseif ft_senstype(hdr.grad, 'ctf275')
    meg_json.ContinuousHeadLocalization = any(strcmp(hdr.chantype, 'headloc')); % CTF specific
    meg_json.Manufacturer             = 'CTF';
    meg_json.ManufacturersModelName   = 'CTF-275';
  elseif ft_senstype(hdr.grad, 'neuromag122')
    meg_json.Manufacturer             = 'Neuromag/Elekta/MEGIN';
    meg_json.ManufacturersModelName   = 'Neuromag-122';
  elseif ft_senstype(hdr.grad, 'neuromag306')
    meg_json.Manufacturer             = 'Neuromag/Elekta/MEGIN';
    % the ManufacturersModelName could be either Vectorview or Triux
  end

  % merge the information specified by the user with that from the data
  % in case fields appear in both, the first input overrules the second
  meg_json = mergeconfig(meg_settings,     meg_json, false);
  meg_json = mergeconfig(generic_settings, meg_json, false);
end % if need_meg_json

%% need_eeg_json
if need_eeg_json
  eeg_json.SamplingFrequency          = hdr.Fs;
  eeg_json.EEGChannelCount            = sum(strcmp(hdr.chantype, 'eeg'));
  eeg_json.EOGChannelCount            = sum(strcmp(hdr.chantype, 'eog'));
  eeg_json.ECGChannelCount            = sum(strcmp(hdr.chantype, 'ecg'));
  eeg_json.EMGChannelCount            = sum(strcmp(hdr.chantype, 'emg'));
  eeg_json.TriggerChannelCount        = sum(strcmp(hdr.chantype, 'trigger'));
  eeg_json.MiscChannelCount           = sum(strcmp(hdr.chantype, 'misc') | strcmp(hdr.chantype, 'unknown'));
  eeg_json.RecordingDuration          = (hdr.nTrials*hdr.nSamples)/hdr.Fs;
  eeg_json.EpochLength                = hdr.nSamples/hdr.Fs;

  % merge the information specified by the user with that from the data
  % in case fields appear in both, the first input overrules the second
  eeg_json = mergeconfig(eeg_settings,     eeg_json, false);
  eeg_json = mergeconfig(generic_settings, eeg_json, false);
end % if need_eeg_json

%% need_ieeg_json
if need_ieeg_json
  ieeg_json.SamplingFrequency          = hdr.Fs;
  ieeg_json.ECOGChannelCount           = sum(strcmp(hdr.chantype, 'ecog'));
  ieeg_json.SEEGChannelCount           = sum(strcmp(hdr.chantype, 'seeg'));
  ieeg_json.EEGChannelCount            = sum(strcmp(hdr.chantype, 'eeg'));
  ieeg_json.EOGChannelCount            = sum(strcmp(hdr.chantype, 'eog'));
  ieeg_json.ECGChannelCount            = sum(strcmp(hdr.chantype, 'ecg'));
  ieeg_json.EMGChannelCount            = sum(strcmp(hdr.chantype, 'emg'));
  ieeg_json.TriggerChannelCount        = sum(strcmp(hdr.chantype, 'trigger'));
  ieeg_json.MiscChannelCount           = sum(strcmp(hdr.chantype, 'misc') | strcmp(hdr.chantype, 'unknown'));
  ieeg_json.RecordingDuration          = (hdr.nTrials*hdr.nSamples)/hdr.Fs;
  ieeg_json.EpochLength                = hdr.nSamples/hdr.Fs;

  % merge the information specified by the user with that from the data
  % in case fields appear in both, the first input overrules the second
  ieeg_json = mergeconfig(ieeg_settings,    ieeg_json, false);
  ieeg_json = mergeconfig(generic_settings, ieeg_json, false);
end

%% need_emg_json
if need_emg_json
  emg_json.SamplingFrequency          = hdr.Fs;
  emg_json.EOGChannelCount            = sum(strcmp(hdr.chantype, 'eog'));
  emg_json.ECGChannelCount            = sum(strcmp(hdr.chantype, 'ecg'));
  emg_json.EMGChannelCount            = sum(strcmp(hdr.chantype, 'emg'));
  emg_json.TriggerChannelCount        = sum(strcmp(hdr.chantype, 'trigger'));
  emg_json.MiscChannelCount           = sum(strcmp(hdr.chantype, 'misc') | strcmp(hdr.chantype, 'unknown'));
  emg_json.RecordingDuration          = (hdr.nTrials*hdr.nSamples)/hdr.Fs;
  emg_json.EpochLength                = hdr.nSamples/hdr.Fs;

  % merge the information specified by the user with that from the data
  % in case fields appear in both, the first input overrules the second
  emg_json = mergeconfig(emg_settings,     emg_json, false);
  emg_json = mergeconfig(generic_settings, emg_json, false);
end

%% need_exg_json
if need_exg_json
  exg_json.SamplingFrequency          = hdr.Fs;
  exg_json.RecordingDuration          = (hdr.nTrials*hdr.nSamples)/hdr.Fs;
  exg_json.EpochLength                = hdr.nSamples/hdr.Fs;

  % merge the information specified by the user with that from the data
  % in case fields appear in both, the first input overrules the second
  exg_json = mergeconfig(exg_settings,     exg_json, false);
  exg_json = mergeconfig(generic_settings, exg_json, false);
end

%% need_audio_json
if need_audio_json
  audio_json.SampleRate    = audio.SampleRate;
  audio_json.Duration      = audio.Duration;
  audio_json.ChannelCount  = audio.NumChannels;

  % merge the information specified by the user with that from the data
  % in case fields appear in both, the first input overrules the second
  audio_json = mergeconfig(audio_settings,   audio_json, false);
  audio_json = mergeconfig(generic_settings, audio_json, false);
end

%% need_video_json
if need_video_json
  ws = warning('off', 'MATLAB:structOnObject');
  video_json = keepfields(struct(video), {'FrameRate', 'Width', 'Height', 'Duration'});
  warning(ws);
  video_json.AudioSampleRate    = audio.SampleRate;
  video_json.AudioDuration      = audio.Duration;
  video_json.AudioChannelCount  = audio.NumChannels;

  % merge the information specified by the user with that from the data
  % in case fields appear in both, the first input overrules the second
  video_json = mergeconfig(video_settings,   video_json, false);
  video_json = mergeconfig(generic_settings, video_json, false);
end

%% need_physio_json
if need_physio_json
  physio_json.SamplingFrequency = hdr.Fs;
  physio_json.StartTime = nan;
  physio_json.Columns = hdr.label;

  % merge the information specified by the user with that from the data
  % in case fields appear in both, the first input overrules the second
  physio_json = mergeconfig(physio_settings,  physio_json, false);
  physio_json = mergeconfig(generic_settings, physio_json, false);
end

%% need_stim_json
if need_stim_json
  stim_json.SamplingFrequency = hdr.Fs;
  stim_json.StartTime = nan;
  stim_json.Columns = hdr.label;

  % merge the information specified by the user with that from the data
  % in case fields appear in both, the first input overrules the second
  stim_json = mergeconfig(stim_settings,    stim_json, false);
  stim_json = mergeconfig(generic_settings, stim_json, false);
end

%% need_eyetracker_json
if need_eyetracker_json
  eyetracker_json.SamplingFrequency = hdr.Fs;
  eyetracker_json.StartTime = nan;
  eyetracker_json.Columns = hdr.label;

  % merge the information specified by the user with that from the data
  % in case fields appear in both, the first input overrules the second
  eyetracker_json = mergeconfig(eyetracker_settings,  eyetracker_json, false);
  eyetracker_json = mergeconfig(generic_settings,     eyetracker_json, false);
end

%% need_motion_json
if need_motion_json
  motion_json.SamplingFrequency = hdr.Fs;
  motion_json.StartTime = nan;
  motion_json.Columns = hdr.label;

  % merge the information specified by the user with that from the data
  % in case fields appear in both, the first input overrules the second
  motion_json = mergeconfig(motion_settings,  motion_json, false);
  motion_json = mergeconfig(generic_settings, motion_json, false);
end

%% need_channels_tsv
if need_channels_tsv

  if isstruct(cfg.channels)
    try
      cfg.channels = struct2table(cfg.channels);
    catch
      ft_error('incorrect specification of cfg.channels');
    end
  end

  % channel information can come from the header and from cfg.channels
  channels_tsv = hdr2table(hdr);
  channels_tsv = merge_table(channels_tsv, cfg.channels, 'name');

  % the default for cfg.channels consists of one row where all values are nan, this needs to be removed
  keep = false(size(channels_tsv.name));
  for i=1:numel(channels_tsv.name)
    keep(i) = ischar(channels_tsv.name{i});
  end
  channels_tsv = channels_tsv(keep,:);

  % do a sanity check on the number of channels for the electrophysiology data types
  if need_meg_json
    type_json = meg_json;
  elseif need_eeg_json
    type_json = eeg_json;
  elseif need_ieeg_json
    type_json = ieeg_json;
  elseif need_emg_json
    type_json = emg_json;
  elseif need_exg_json
    type_json = exg_json;
  end
  fn = fieldnames(type_json);
  fn = fn(endsWith(fn, 'ChannelCount'));
  jsoncount = 0;
  for i=1:numel(fn)
    if ~isempty(type_json.(fn{i}))
      jsoncount = jsoncount + type_json.(fn{i});
    end
  end
  if size(channels_tsv,1)~=jsoncount
    ft_warning('incorrect specification of the channel count: %d in the json, %d in the tsv', jsoncount, size(channels_tsv,1));
  end
end % if need_channels_tsv

%% need_electrodes_tsv
if need_electrodes_tsv

  if isstruct(cfg.electrodes)
    try
      cfg.electrodes = struct2table(cfg.electrodes);
    catch
      ft_error('incorrect specification of cfg.electrodes.%s', fn{i});
    end
  end

  % electrode details can be specified in cfg.elec, data.elec or in cfg.electrodes
  electrodes_tsv = elec2table(elec);
  electrodes_tsv = merge_table(electrodes_tsv, cfg.electrodes, 'name');

  % the default for cfg.electrodes consists of one row where all values are nan, this needs to be removed
  keep = false(size(electrodes_tsv.name));
  for i=1:numel(electrodes_tsv.name)
    keep(i) = ischar(electrodes_tsv.name{i});
  end
  electrodes_tsv = electrodes_tsv(keep,:);

end % need_electrodes_tsv

%% need_coordsystem_json
if need_coordsystem_json
  if isfield(hdr, 'grad') && ft_senstype(hdr.grad, 'ctf')
    % coordinate system for MEG sensors
    coordsystem_json.MEGCoordinateSystem            = 'CTF';
    coordsystem_json.MEGCoordinateUnits             = 'cm';
    coordsystem_json.MEGCoordinateSystemDescription = 'CTF head coordinates, orientation ALS, origin between the ears';
    % coordinate system for head localization coils
    coordsystem_json.HeadCoilCoordinates                 = []; % see below
    coordsystem_json.HeadCoilCoordinateSystem            = 'CTF';
    coordsystem_json.HeadCoilCoordinateUnits             = 'cm';
    coordsystem_json.HeadCoilCoordinateSystemDescription = 'CTF head coordinates, orientation ALS, origin between the ears';
    if isempty(coordsystem_json.HeadCoilCoordinates)
      % get the positions from the dataset header
      label = cellstr(hdr.orig.hc.names);
      position = hdr.orig.hc.head;
      for i=1:numel(label)
        coordsystem_json.HeadCoilCoordinates.(fixname(label{i})) = position(:,i)';
      end
    end
  elseif isfield(hdr, 'grad') && ft_senstype(hdr.grad, 'neuromag')
    % coordinate system for MEG sensors
    coordsystem_json.MEGCoordinateSystem            = 'Neuromag';
    coordsystem_json.MEGCoordinateUnits             = 'm';
    coordsystem_json.MEGCoordinateSystemDescription = 'Neuromag head coordinates, orientation RAS, origin between the ears';
    % coordinate system for head localization coils
    coordsystem_json.HeadCoilCoordinates                 = []; % FIXME it might be possible to get these from the dataset header
    coordsystem_json.HeadCoilCoordinateSystem            = 'Neuromag';
    coordsystem_json.HeadCoilCoordinateUnits             = 'm';
    coordsystem_json.HeadCoilCoordinateSystemDescription = 'Neuromag head coordinates, orientation RAS, origin between the ears';
  else
    ft_warning('coordsystem handling not yet supported for this data, you MUST specify cfg.coordsystem');
    coordsystem_json = table();
  end

  % merge the information specified by the user with that from the data
  % in case fields appear in both, the first input overrules the second
  coordsystem_json = mergeconfig(coordsystem_settings, coordsystem_json, false); % FIXME the order of precedence is different here
end % if need_coordsystem_json

%% need_events_tsv
if need_events_tsv

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% align the presentation events with BOLD volumes or MEG/EEG triggers
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if ~isempty(cfg.presentationfile) && need_mri_json
    % the code can either align presentation and BOLD volumes, or use the user-supplied events that must be aligned
    assert(isempty(cfg.events) || (isstruct(cfg.events) && numel(fieldnames(cfg.events))==0), 'this is mutually exclusive with specifying cfg.events')

    % read the events from the presentation file
    presentation = ft_read_event(cfg.presentationfile);

    % merge the information with the json sidecar file
    % in case fields appear in both, the first input overrules the second
    tmp = mergeconfig(mri_json, read_json(corresponding_json(cfg.outputfile)), false);
    assert(~isempty(tmp.RepetitionTime), 'you must specify cfg.mri.RepetitionTime');

    % create a header structure that represents the fMRI timeseries
    hdr.Fs = 1/tmp.RepetitionTime;
    hdr.nSamples = mri.dim(4);

    % create a event structure with one trigger for each BOLD volume
    trigger = [];
    for i=1:hdr.nSamples
      trigger(i).type   = 'volume';
      trigger(i).sample = i;
    end

    % align the presentation events with the triggers
    cfg.events = align_presentation(presentation, cfg.presentation, trigger, cfg.trigger, hdr, istrue(cfg.feedback));

  elseif ~isempty(cfg.presentationfile) && (need_meg_json || need_eeg_json || need_ieeg_json || need_emg_json || need_exg_json)
    % the code can either align presentation and trigger channel, or use the user-supplied events that must be aligned
    assert(isempty(cfg.events) || (isstruct(cfg.events) && numel(fieldnames(cfg.events))==0), 'this is mutually exclusive with specifying cfg.events')

    % read the events from the presentation file
    presentation = ft_read_event(cfg.presentationfile);

    % align the presentation events with the triggers
    cfg.events = align_presentation(presentation, cfg.presentation, trigger, cfg.trigger, hdr, istrue(cfg.feedback));

    % also include all triggers from the MEG/EEG dataset as events
    cfg.events = appendstruct(cfg.events, trigger);

  end % if presentationfile

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% convert the events into a table
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if istable(cfg.events) && all(ismember({'onset', 'duration'}, fieldnames(cfg.events)))
    % use the events table as it is
    events_tsv = cfg.events;
  elseif istable(cfg.events) && all(ismember({'begsample', 'endsample', 'offset'}, fieldnames(cfg.events)))
    % it is a "trl" matrix formatted as table, use it as it is
    events_tsv = cfg.events;
  elseif isstruct(cfg.events) && ~isempty(cfg.events) && numel(fieldnames(cfg.events))>0
    % it is the output from FT_READ_EVENT
    if exist('hdr', 'var')
      events_tsv = event2table(hdr, cfg.events);
    else
      events_tsv = event2table([], cfg.events);
    end
  elseif ismatrix(cfg.events) && ~isempty(cfg.events) && numel(fieldnames(cfg.events))>0
    % it is a "trl" matrix formatted as numeric array, convert it to an events table
    begsample = cfg.events(:,1);
    endsample = cfg.events(:,2);
    offset    = cfg.events(:,3); % this is not used for the events.tsv
    if any(offset~=0)
      ft_warning('the offset in the trl matrix is ignored');
    end
    if size(trl,2)>3
      ft_warning('additional columns in the trl matrix are ignored');
    end
    % convert to the required fields
    onset     = (begsample-1)/hdr.Fs;
    duration  = (endsample-begsample+1)/hdr.Fs;
    events_tsv = table(onset, duration);
  elseif exist('trigger', 'var')
    % convert the triggers from FT_READ_EVENT into a table
    if exist('hdr', 'var')
      events_tsv = event2table(hdr, trigger);
    else
      events_tsv = event2table([], trigger);
    end
  elseif ~isempty(cfg.presentationfile)
    % read the presentation file and convert into a table
    events_tsv = event2table([], ft_read_event(cfg.presentationfile));
  else
    ft_warning('no events were specified');
    % make an empty table with columns for onset and duration
    onset    = [];
    duration = [];
    events_tsv = table(onset, duration);
  end

  if isempty(events_tsv)
    ft_warning('no events found');
    need_events_tsv = false;

  else
    % ensure that column names are in lower case
    events_tsv.Properties.VariableNames = lower(events_tsv.Properties.VariableNames);

    % ensure that the onset and duration appear as the first two columns
    order = nan(1, size(events_tsv,2));
    order(1) = find(strcmp(events_tsv.Properties.VariableNames, 'onset'));
    order(2) = find(strcmp(events_tsv.Properties.VariableNames, 'duration'));
    order(3:end) = setdiff(1:size(events_tsv,2), order([1 2]));
    events_tsv = events_tsv(:,order);

    % sort the events ascending on the onset
    events_tsv = sortrows(events_tsv, 'onset');
  end

  if ~isempty(cfg.presentationfile) && need_mri_json
    % rename the column 'sample' into 'volume'
    events_tsv.Properties.VariableNames('sample') = {'volume'};
  end

end % if need_events_tsv

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% write or copy the data to the output file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[p, f, x] = fileparts(cfg.outputfile);
isdir_or_mkdir(p);

switch cfg.method
  case 'decorate'
    % look at the user's specification of cfg.datatype
    if strcmp(cfg.datatype, 'events')
      % add the TSV file extension, this is needed for behavioral data represented in scans.tsv
      % the events.tsv file will be written further down
      [p, f, x] = fileparts(cfg.outputfile);
      cfg.outputfile = fullfile(p, [f '.tsv']);
    end

  case 'convert'
    % the output depends on the type of input data
    switch typ
      case {'nifti', 'nifti2', 'nifti_fsl', 'dicom', 'volume'}
        [p, f, x] = fileparts(cfg.outputfile);
        if ~isequal(x, '.nii')
          cfg.outputfile = fullfile(p, [f '.nii']);
        end
        if istrue(cfg.deface)
          % deface the anatomical MRI
          ft_info('defacing anatomical MRI');
          if istrue(cfg.feedback)
            tmpcfg = [];
            tmpcfg.figurename = sprintf('ORIGINAL - %s', f);
            tmpcfg.showcallinfo = 'no';
            tmpcfg.colorbar = 'no';
            ft_sourceplot(tmpcfg, mri);
          end
          tmpcfg = [];
          tmpcfg.method = 'spm';
          if ~isfield(mri, 'coordsys')
            % try to determine it from the BIDS filename
            coordsys = get_token(cfg.outputfile, 'space');
            if ~isempty(coordsys)
              mri.coordsys = coordsys;
            else
              mri.coordsys = 'unknown';
            end
          end
          mri = ft_defacevolume(tmpcfg, mri);
          if istrue(cfg.feedback)
            tmpcfg = [];
            tmpcfg.figurename = sprintf('DEFACED - %s', f);
            tmpcfg.showcallinfo = 'no';
            tmpcfg.colorbar = 'no';
            ft_sourceplot(tmpcfg, mri);
          end
        end
        ft_info('writing %s\n', cfg.outputfile);
        ft_write_mri(cfg.outputfile, mri, 'dataformat', 'nifti');

      case {'ctf_ds', 'ctf_meg4', 'ctf_res4', 'ctf151', 'ctf275', 'neuromag_fif', 'neuromag122', 'neuromag306'}
        ft_error('please use a system specific tool for converting MEG datasets');

      case {'presentation_log'}
        % the events.tsv file will be written further down
        [p, f, x] = fileparts(cfg.outputfile);
        cfg.outputfile = fullfile(p, [f '.tsv']);

      otherwise
        % look at the user's specification of cfg.datatype
        switch cfg.datatype
          case {'eeg', 'ieeg', 'emg', 'exg'}
            % write the data in BrainVision core file format
            [p, f, x] = fileparts(cfg.outputfile);
            cfg.outputfile = fullfile(p, [f '.vhdr']);
            ft_info('writing %s\n', cfg.outputfile);
            ft_write_data(cfg.outputfile, dat, 'dataformat', 'brainvision_eeg', 'header', hdr, 'event', trigger);
          case {'physio', 'stim', 'eyetracker', 'motion'}
            % write the data according to the Stim and Physio format as specified at
            % https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/06-physiological-and-other-continuous-recordings.html
            [p, f, x] = fileparts(cfg.outputfile);
            cfg.outputfile = fullfile(p, [f '.tsv']);
            ft_info('writing %s\n', cfg.outputfile);
            writematrix(dat', cfg.outputfile, 'FileType', 'text', 'Delimiter', '\t'); % without headers, the JSON will be written further down
          case {'events'}
            % add the TSV file extension, this is needed for behavioral data represented in scans.tsv
            % the events.tsv file will be written further down
            [p, f, x] = fileparts(cfg.outputfile);
            cfg.outputfile = fullfile(p, [f '.tsv']);
          otherwise
            ft_error('cannot determine how to write the data')
        end
    end % switch typ

  case 'copy'
    [~, ~, xin] = fileparts(cfg.dataset);
    [p, ~, xout] = fileparts(cfg.outputfile);
    if ~strcmp(xin, xout)
      ft_error('input and output filename extension do not match');
    end

    switch typ
      case {'dicom'}
        ft_error('DICOM files must be converted to NIfTI for BIDS compliance');

      case {'eyelink_edf', 'eyelink_asc'}
        ft_error('Eyelink files must be converted to TSV for BIDS compliance');

      case {'ctf_ds', 'ctf_meg4', 'ctf_res4', 'ctf151', 'ctf275'}
        % the data consists of a directory with multiple files inside
        ft_info('copying %s to %s\n', cfg.dataset, cfg.outputfile);
        copy_ctf_files(cfg.dataset, cfg.outputfile, false);

      case {'brainvision_vhdr', 'brainvision_vmrk', 'brainvision_eeg', 'brainvision_dat', 'brainvision_seg'}
        % the data consists of three files and the header file contains pointers to the markers and data
        ft_info('copying %s to %s\n', cfg.dataset, cfg.outputfile);
        copy_brainvision_files(cfg.dataset, cfg.outputfile, false);

      otherwise
        ft_info('copying %s to %s\n', cfg.dataset, cfg.outputfile);
        copyfile(cfg.dataset, cfg.outputfile);
    end

  otherwise
    ft_error('unsupported value for cfg.method')
end % switch method


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% write the metadata to the json files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% each of these has a corresponding json file
modality = {'mri', 'meg', 'eeg', 'ieeg', 'emg', 'exg', 'audio', 'video', 'eyetracker', 'physio', 'stim', 'motion', 'coordsystem'};
for i=1:numel(modality)
  if eval(sprintf('need_%s_json', modality{i}))
    modality_json = eval(sprintf('%s_json', modality{i}));
    modality_json = remove_empty(modality_json);

    if strcmp(modality{i}, 'coordsystem')
      [p, f] = fileparts(cfg.outputfile);
      f = remove_entity(f, 'task');     % remove _task-something
      f = remove_entity(f, 'acq');      % remove _acq-something
      f = remove_entity(f, 'ce');       % remove _ce-something
      f = remove_entity(f, 'rec');      % remove _rec-something
      f = remove_entity(f, 'dir');      % remove _dir-something
      f = remove_entity(f, 'run');      % remove _run-something
      f = remove_entity(f, 'mod');      % remove _mod-something
      f = remove_entity(f, 'echo');     % remove _echo-something
      f = remove_entity(f, 'proc');     % remove _proc-something
      f = remove_datatype(f);           % remove _meg, _eeg, etc.
      filename = fullfile(p, [f '_coordsystem.json']);
    else
      % just replace the extension with json
      filename = corresponding_json(cfg.outputfile);
    end

    if isfile(filename)
      existing = read_json(filename);
    else
      existing = [];
    end

    switch cfg.writejson
      case 'yes'
        if ~isempty(existing)
          ft_warning('not overwriting the existing and non-empty file ''%s''', filename);
        else
          write_json(filename, modality_json);
        end
      case 'replace'
        write_json(filename, modality_json);
      case 'merge'
        write_json(filename, mergeconfig(modality_json, existing, false))
      case 'no'
        % do nothing
        ft_info('not writing %s\n', filename);
      otherwise
        ft_error('incorrect option for cfg.writejson');
    end % switch writejson
  end % if need_xxx_json
end % for each modality


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% write the metadata to the tsv files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% each of these has a corresponding tsv file
modality = {'channels', 'electrodes', 'events'};
for i=1:numel(modality)
  if eval(sprintf('need_%s_tsv', modality{i}))
    modality_tsv = eval(sprintf('%s_tsv', modality{i}));

    if strcmp(modality{i}, 'electrodes')
      [p, f] = fileparts(cfg.outputfile);
      f = remove_entity(f, 'task');     % remove _task-something
      f = remove_entity(f, 'acq');      % remove _acq-something
      f = remove_entity(f, 'ce');       % remove _ce-something
      f = remove_entity(f, 'rec');      % remove _rec-something
      f = remove_entity(f, 'dir');      % remove _dir-something
      f = remove_entity(f, 'run');      % remove _run-something
      f = remove_entity(f, 'mod');      % remove _mod-something
      f = remove_entity(f, 'echo');     % remove _echo-something
      f = remove_entity(f, 'proc');     % remove _proc-something
      f = remove_datatype(f);           % remove _meg, _eeg, etc.
      filename = fullfile(p, sprintf('%s_%s.tsv', f, modality{i}));
    else
      [p, f] = fileparts(cfg.outputfile);
      f = remove_datatype(f); % remove _bold, _meg, etc.
      filename = fullfile(p, sprintf('%s_%s.tsv', f, modality{i}));
    end

    if isfile(filename)
      existing = read_tsv(filename);
    else
      existing = [];
    end

    switch cfg.writetsv
      case 'yes'
        if ~isempty(existing)
          ft_warning('not overwriting the existing and non-empty file ''%s''', filename);
        else
          write_tsv(filename, modality_tsv);
        end
      case 'replace'
        write_tsv(filename, modality_tsv);
      case 'merge'
        if ~isempty(existing)
          modality_tsv  = output_compatible(modality_tsv);
          existing      = output_compatible(existing);
          if strcmp(modality{i}, 'events')
            % merge complete rows
            modality_tsv = merge_table(modality_tsv, existing);
          else
            % use the channel name as the unique key
            modality_tsv = merge_table(modality_tsv, existing, 'name');
          end
        end
        write_tsv(filename, modality_tsv);
      case 'no'
        % do nothing
        ft_info('not writing %s\n', filename);
      otherwise
        ft_error('incorrect option for cfg.writetsv');
    end % switch
  end % if needd_xxx_tsv
end % for each modality


if ~isempty(cfg.bidsroot)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % update the dataset_description
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  filename = fullfile(cfg.bidsroot, 'dataset_description.json');

  if isfile(filename)
    existing = read_json(filename);
  else
    existing = [];
  end

  switch cfg.writejson
    case 'yes'
      if ~isempty(existing)
        ft_warning('not overwriting the existing and non-empty file ''%s''', filename);
      else
        write_json(filename, dataset_description_settings);
      end
    case 'replace'
      write_json(filename, dataset_description_settings);
    case 'merge'
      write_json(filename, mergeconfig(dataset_description_settings, existing, false));
    case 'no'
      % do nothing
    otherwise
      ft_error('incorrect option for cfg.writejson');
  end % switch

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % update the participants.tsv
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  this = table();
  this.participant_id = ['sub-' cfg.sub];
  fn = fieldnames(cfg.participants);
  for i=1:numel(fn)
    % write [] as 'n/a'
    % write nan as 'n/a'
    % write boolean as 'True' or 'False'
    this.(fn{i}) = output_compatible(cfg.participants.(fn{i}));
  end

  filename = fullfile(cfg.bidsroot, 'participants.tsv');
  if isfile(filename)
    participants = read_tsv(filename);
    % the participants.tsv is always merged
    participants = merge_table(participants, this, 'participant_id');
  else
    participants = this;
  end

  % write the updated file back to disk
  write_tsv(filename, participants);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % update the scans.tsv
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ~isempty(cfg.ses)
    % construct the output filename, with session directory
    filename = fullfile(cfg.bidsroot, ['sub-' cfg.sub], ['ses-' cfg.ses], ['sub-' cfg.sub '_' 'ses-' cfg.ses '_scans.tsv']);
  else
    % construct the output filename, without session directory
    filename = fullfile(cfg.bidsroot, ['sub-' cfg.sub], ['sub-' cfg.sub '_scans.tsv']);
  end

  this = table();
  [~, f, x] = fileparts(cfg.outputfile);
  this.filename = fullfile(datatype2dirname(cfg.datatype), [f x]);
  fn = fieldnames(cfg.scans);
  for i=1:numel(fn)
    % write [] as 'n/a'
    % write nan as 'n/a'
    % write boolean as 'True' or 'False'
    this.(fn{i}) = output_compatible(cfg.scans.(fn{i}));
  end

  if isfile(filename)
    scans_tsv = read_tsv(filename);
    % the scans.tsv is always merged
    scans_tsv = merge_table(scans_tsv, this, 'filename');
  else
    scans_tsv = this;
  end

  % write the updated file back to disk
  write_tsv(filename, scans_tsv);

end % if bidsroot

% do not return an output variable if not requested
if ~ft_nargout
  clear cfg
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = get_token(f, tok)
% ensure that it is only the filename
[p, f] = fileparts(f);
pieces = tokenize(f, '_');
val = '';
for i=1:numel(pieces)
  if startsWith([pieces{i} '-'], tok)
    val = pieces{i}((numel(tok)+2):end);
    return
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = add_entity(f, key, val)
if ~isempty(val)
  if isscalar(val)
    % this applies specifically to run
    val = num2str(val);
  end
  f = [f '_' key '-' val];
end

function f = remove_entity(f, key)
part = regexp(f, sprintf('_%s-[a-zA-Z0-9]*', key), 'split');
if numel(part)>1 && ~isempty(part{2})
  f = [part{1} part{2}];
else
  f = part{1};
end

function f = add_datatype(f, typ)
f = [f '_' typ];

function f = remove_datatype(f)
typ = {'FLAIR', 'FLASH', 'PD', 'PDT2', 'PDmap', 'T1map', 'T1rho', 'T1w', 'T2map', 'T2star', 'T2w', 'angio', 'bold', 'bval', 'bvec', 'channels', 'coordsystem', 'defacemask', 'dwi', 'eeg', 'epi', 'events', 'fieldmap', 'headshape', 'ieeg', 'inplaneT1', 'inplaneT2', 'magnitude', 'magnitude1', 'magnitude2', 'meg', 'phase1', 'phase2', 'phasediff', 'photo', 'sbref', 'physio', 'stim', 'emg', 'exg', 'audio', 'video', 'eyetracker', 'motion'};
for i=1:numel(typ)
  if endsWith(f, ['_' typ{i}])
    f = f(1:end-length(typ{i})-1); % also the '_'
    return
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION convert hdr structure into table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tab = hdr2table(hdr)
if isempty(hdr)
  tab = table();
else
  name = hdr.label(:);
  type = hdr.chantype(:);
  units = hdr.chanunit(:);
  sampling_frequency = repmat(hdr.Fs, hdr.nChans, 1);
  tab = table(name, type, units, sampling_frequency);;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION convert elec structure into table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tab = elec2table(elec)
if isempty(elec)
  tab = table();
else
  name = elec.label(:);
  x = elec.elecpos(:,1);
  y = elec.elecpos(:,2);
  z = elec.elecpos(:,3);
  tab = table(name, x, y, z);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION convert event structure into table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tab = event2table(hdr, event)
if isempty(event)
  tab = table();
else
  for i=1:numel(event)
    % normalize missing values
    if isempty(event(i).sample)
      event(i).sample = NaN;
    end
    if isempty(event(i).duration)
      event(i).duration = NaN;
    end
    if isempty(event(i).value)
      event(i).value = NaN;
    end
  end
  if isempty(hdr) && isfield(event, 'timestamp')
    ft_warning('assuming a timestamp resolution of 0.1 milliseconds');
    onset        = ([event.timestamp]/10000)';     % the resolution is 0.1 milliseconds
    duration     = zeros(size(onset));             % in seconds
    type         = {event.type}';
    value        = {event.value}';
  else
    onset        = (([event.sample]-1)/hdr.Fs)';   % in seconds
    duration     = ([event.duration]/hdr.Fs)';     % in seconds
    sample       = ([event.sample])';              % in samples, the first sample of the file is 1
    type         = {event.type}';
    value        = {event.value}';
  end
  if all(cellfun(@isnumeric, type))
    % this can be an array of strings or values
    type = cell2mat(type);
  end
  if all(cellfun(@isnumeric, value))
    % this can be an array of strings or values
    value = cell2mat(value);
  end
  if exist('sample', 'var')
    tab = table(onset, duration, sample, type, value);
  else
    tab = table(onset, duration, type, value);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tsv = read_tsv(filename)
ft_info('reading %s\n', filename);
tsv = readtable(filename, 'Delimiter', 'tab', 'FileType', 'text', 'TreatAsEmpty', 'n/a', 'ReadVariableNames', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_tsv(filename, tsv)
ft_info('writing %s\n', filename);
fn = tsv.Properties.VariableNames;
for i=1:numel(fn)
  % write [] as 'n/a'
  % write nan as 'n/a'
  % write boolean as 'True' or 'False'
  tsv.(fn{i}) = output_compatible(tsv.(fn{i}));
end
writetable(tsv, filename, 'Delimiter', 'tab', 'FileType', 'text');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function json = read_json(filename)
ft_info('reading %s\n', filename);
ft_hastoolbox('jsonlab', 1);
json = loadjson(filename);
json = ft_struct2char(json); % convert strings into char-arrays

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_json(filename, json)
ft_info('writing %s\n', filename);
json = remove_empty(json);
json = sort_fields(json);
json = ft_struct2char(json); % convert strings into char-arrays
ft_hastoolbox('jsonlab', 1);
% see also the output_compatible helper function
% write nan as 'n/a'
% write boolean as True/False
str = savejson('', json, 'NaN', '"n/a"', 'ParseLogical', true);
fid = fopen_or_error(filename, 'w');
fwrite(fid, str);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filename = corresponding_json(filename)
[p, f, x] = fileparts(filename);
if isequal(x, '.gz') && endsWith(f, '.nii')
  % it is a gzip compressed NIfTI file, remove the .nii from the file name
  f = f(1:end-4);
end
filename = fullfile(p, [f '.json']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = merge_vector(x, y)
assert(isequal(size(x), size(y)));
for i=1:numel(x)
  if isnumeric(x) && isnumeric(y) && isnan(x(i)) && ~isnan(y(i))
    x(i) = y(i);
  end
  if iscell(x) && iscell(y) && isempty(x{i}) && ~isempty(y{i})
    x{i} = y{i};
  end
  if iscell(x) && isnumeric(y) && isempty(x{i}) && ~isnan(y{i})
    x{i} = y(i);
  end
  if iscell(x) && iscell(y) && isequal(x{i}, 'unknown') && ~isequal(y{i}, 'unknown')
    x{i} = y{i};
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sel = select_event(event, eventtype, eventvalue)
% this code is copied from FT_TRIALFUN_GENERAL

% start by selecting all events
sel = true(1, length(event)); % this should be a row vector

% select all events of the specified type
if ~isempty(eventtype)
  for i=1:numel(event)
    sel(i) = sel(i) && ismatch(event(i).type, eventtype);
  end
end

% select all events with the specified value
if ~isempty(eventvalue)
  for i=1:numel(event)
    sel(i) = sel(i) && ismatch(event(i).value, eventvalue);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = remove_empty(s)
if isempty(s)
  return
else
  fn = fieldnames(s);
  fn = fn(structfun(@isempty, s));
  s = removefields(s, fn);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = sort_fields(x)
fn = fieldnames(x);
fn = sort(fn);
y = struct();
for i=1:numel(fn)
  y.(fn{i}) = x.(fn{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dir = datatype2dirname(typ)
% see https://bids-specification.readthedocs.io/en/stable/99-appendices/04-entity-table.html
% emg, exg, eyetracker, motion, audio, and video are not part of the official specification
switch typ
  case {'T1w' 'T2w' 'T1rho' 'T1map' 'T2map' 'T2star' 'FLAIR' 'FLASH' 'PD' 'PDmap' 'PDT2' 'inplaneT1' 'inplaneT2' 'angio' 'defacemask'}
    dir = 'anat';
  case {'bold' 'sbref'} % this could also include 'events' 'physio' 'stim'
    dir = 'func';
  case {'dwi' 'bvec' 'bval'}
    dir = 'dwi';
  case {'phasediff' 'phase1' 'phase2' 'magnitude1' 'magnitude2' 'magnitude' 'fieldmap' 'epi'}
    dir = 'fmap';
  case {'events' 'stim' 'physio' 'eyetracker' 'motion' 'audio' 'video'} % stim and physio could also be stored in 'func'
    dir = 'beh';
  case {'meg'} % this could also include 'events'
    dir = 'meg';
  case {'eeg'} % this could also include 'events'
    dir = 'eeg';
  case {'ieeg'} % this could also include 'events'
    dir = 'ieeg';
  case {'emg'} % this is not part of the official specification
    dir = 'emg';
  case {'exg'} % this is not part of the official specification
    dir = 'exg';
  otherwise
    ft_error('unrecognized data type "%s"', typ);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = output_compatible(val)
if istable(val)
  fn = val.Properties.VariableNames;
  for i=1:numel(fn)
    val.(fn{i}) = output_compatible(val.(fn{i}));
  end
elseif iscell(val)
  % use recursion to make all elements compatible
  val = cellfun(@output_compatible, val, 'UniformOutput', false);
elseif isnumeric(val) && numel(val)>1 && any(isnan(val))
  % convert and use recursion to make all elements compatible
  val = num2cell(val);
  val = cellfun(@output_compatible, val, 'UniformOutput', false);
else
  % write [] as 'n/a'
  % write nan as 'n/a'
  % write boolean as 'True' or 'False'
  if isempty(val)
    val = 'n/a';
  elseif isnan(val)
    val = 'n/a';
  elseif islogical(val)
    if val
      val = 'True';
    else
      val = 'False';
    end
  end
end
