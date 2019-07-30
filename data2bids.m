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
%   cfg.method                  = string, can be 'decorate', 'convert' or 'copy', see below (default is automatic)
%   cfg.dataset                 = string, filename of the input data
%   cfg.outputfile              = string, optional filename for the output data
%   cfg.mri.deface              = string, 'yes' or 'no' (default = 'no')
%   cfg.mri.writesidecar        = string, 'yes', 'replace', 'merge' or 'no' (default = 'yes')
%   cfg.meg.writesidecar        = string, 'yes', 'replace', 'merge' or 'no' (default = 'yes')
%   cfg.eeg.writesidecar        = string, 'yes', 'replace', 'merge' or 'no' (default = 'yes')
%   cfg.ieeg.writesidecar       = string, 'yes', 'replace', 'merge' or 'no' (default = 'yes')
%   cfg.events.writesidecar     = string, 'yes', 'replace', 'merge' or 'no' (default = 'yes')
%   cfg.coordystem.writesidecar = string, 'yes', 'replace', 'merge' or 'no' (default = 'yes')
%   cfg.channels.writesidecar   = string, 'yes', 'replace', 'merge' or 'no' (default = 'yes')
%   cfg.electrodes.writesidecar = string, 'yes', 'replace', 'merge' or 'no' (default = 'yes')
%
% This function starts from existing data file on disk or from a FieldTrip compatible
% data structure in MATLAB memory that is passed as the second input argument.
% Depending on cfg.method it will add the sidecar files, copy the dataset and add
% sidecar files, or convert the dataset and add the sidecar files. Each of the
% methods is discussed here.
%
% DECORATE - data2bids will read the header details and events from the data and write
% the appropriate sidecar files alongside the existing dataset. You would use this to
% obtain the sidecar files for a dataset that already has the correct BIDS name.
%
% CONVERT - data2bids will read the input data (or use the specified input data) and
% write it to a new output file that is BIDS compliant. The output format is NIFTI
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
% In case cfg.dataset points to a NIFTI file, or in case you pass a preprocessed MRI
% as input data structure, you can specify cfg.mri.dicomfile to read the detailed MR
% scanner and sequence details from the DICOM file. This will be used to fill in the
% details of the corresponding JSON file.
%   cfg.mri.dicomfile           = string, filename of a matching DICOM file for header details (default = [])
%
% You can specify cfg.events.trl as a Nx3 matrix with the trial definition (see
% FT_DEFINETRIAL) or as a MATLAB table. When specified as table, the first three
% columns containing integer values corresponding to the begsample, endsample and
% offset, the additional colums can be of another type and have any name. If you do
% not specify the trial definition, the events will be read from the MEG/EEG/iEEG
% dataset. Events from the trial definition or from the data will be written to
% events.tsv.
%   cfg.events.trl              = trial definition, see also FT_DEFINETRIAL
%
% You can specify cfg.presentationfile with the name of a NBS presentation log file,
% which will be aligned with the data based on triggers (MEG/EEG/iEEG) or based on
% the volumes (fMRI). To indicate how triggers (in MEG/EEG/iEEG) or volumes (in fMRI)
% match the presentation events, you should also specify the mapping between them.
% Events from the presentation log file will be written to events.tsv.
%   cfg.presentationfile        = string, optional filename for the presentation log file
%   cfg.trigger.eventtype       = string (default = [])
%   cfg.trigger.eventvalue      = string or number
%   cfg.presentation.eventtype  = string (default = [])
%   cfg.presentation.eventvalue = string or number
%   cfg.presentation.skip       = 'last'/'first'/'none'
%
% For EEG and iEEG data you can specify an electrode definition according to
% FT_DATATYPE_SENS as an "elec" field in the input data, or you can specify it as
% cfg.elec or you can specify a filename with electrode information.
%   cfg.elec                     = structure with electrode positions or filename, see FT_READ_SENS
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
% FT_PREPROCESSING, FT_READ_MRI

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
cfg = ft_checkconfig(cfg, 'forbidden', 'deface');                           % should be cfg.mri.deface
cfg = ft_checkconfig(cfg, 'renamedval', {'native', 'no', 'convert'});
cfg = ft_checkconfig(cfg, 'renamedval', {'native', 'yes', 'copy'});
cfg = ft_checkconfig(cfg, 'renamed', {'native', 'method'});

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
cfg.trigger.event           = ft_getopt(cfg.trigger, 'event');

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
% cfg.space     = ft_getopt(cfg, 'space'); % FIXME
cfg.datatype  = ft_getopt(cfg, 'datatype');

cfg.mri                     = ft_getopt(cfg, 'mri');
cfg.mri.deface              = ft_getopt(cfg.mri, 'deface', 'no');             % deface the anatomical MRI
cfg.mri.dicomfile           = ft_getopt(cfg.mri, 'dicomfile');                % get the details from one of the original DICOM files
cfg.mri.writesidecar        = ft_getopt(cfg.mri, 'writesidecar', 'yes');      % whether to write the sidecar file

cfg.meg                     = ft_getopt(cfg, 'meg');
cfg.meg.writesidecar        = ft_getopt(cfg.meg, 'writesidecar', 'yes');      % whether to write the sidecar file

cfg.eeg                     = ft_getopt(cfg, 'eeg');
cfg.eeg.writesidecar        = ft_getopt(cfg.eeg, 'writesidecar', 'yes');      % whether to write the sidecar file

cfg.ieeg                    = ft_getopt(cfg, 'ieeg');
cfg.ieeg.writesidecar       = ft_getopt(cfg.ieeg, 'writesidecar', 'yes');     % whether to write the sidecar file

cfg.channels                = ft_getopt(cfg, 'channels');
cfg.channels.writesidecar   = ft_getopt(cfg.channels, 'writesidecar', 'yes'); % whether to write the sidecar file

cfg.electrodes              = ft_getopt(cfg, 'electrodes');
cfg.electrodes.writesidecar = ft_getopt(cfg.electrodes, 'writesidecar', 'yes'); % whether to write the sidecar file

cfg.events                  = ft_getopt(cfg, 'events');
cfg.events.trl              = ft_getopt(cfg.events, 'trl');                   % this can contain the trial definition as Nx3 array or as table
cfg.events.writesidecar     = ft_getopt(cfg.events, 'writesidecar', 'yes');   % whether to write the sidecar file

cfg.coordsystem              = ft_getopt(cfg, 'coordsystem');
cfg.coordsystem.writesidecar = ft_getopt(cfg.coordsystem, 'writesidecar', 'yes');

%% Dataset description

cfg.dataset_description                     = ft_getopt(cfg, 'dataset_description'                       );
cfg.dataset_description.writesidecar        = ft_getopt(cfg.dataset_description, 'writesidecar', 'yes'   );
cfg.dataset_description.Name                = ft_getopt(cfg.dataset_description, 'Name'                  ); % REQUIRED. Name of the dataset.
cfg.dataset_description.BIDSVersion         = ft_getopt(cfg.dataset_description, 'BIDSVersion'           ); % REQUIRED. The version of the BIDS standard that was used.
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
cfg.Manufacturer                      = ft_getopt(cfg, 'Manufacturer'                ); % OPTIONAL. Manufacturer of the MEG system ("CTF", ""Neuromag/Elekta"", ""4D/BTi"", ""KIT/Yokogawa"", ""ITAB"", "KRISS", "Other")
cfg.ManufacturersModelName            = ft_getopt(cfg, 'ManufacturersModelName'      ); % OPTIONAL. Manufacturer's designation of the MEG scanner model (e.g. "CTF-275"). See "Appendix VII" with preferred names
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
cfg.mri.PhaseEncodingDirection        = ft_getopt(cfg.mri, 'PhaseEncodingDirection'         ); % Possible values = [];                     % "i", "j", "k", "i-", "j-", "k-". The letters "i", "j", "k" correspond to the first, second and third axis of the data in the NIFTI file. The polarity of the phase encoding is assumed to go from zero index to maximum index unless "-" sign is present (then the order is reversed - starting from the highest index instead of zero). PhaseEncodingDirection is defined as the direction along which phase is was modulated which may result in visible distortions. Note that this is not the same as the DICOM term InPlanePhaseEncodingDirection which can have "ROW" or "COL" values. This parameter is REQUIRED if corresponding fieldmap data is present or when using multiple runs with different phase encoding directions (which can be later used for field inhomogeneity correction).
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
cfg.ieeg.Manufacturer                    = ft_getopt(cfg, 'Manufacturer'                   ); % RECOMMENDED. Manufacturer of the amplifier system  (e.g. "TDT, Blackrock")
cfg.ieeg.ManufacturersModelName          = ft_getopt(cfg, 'ManufacturersModelName'         ); % RECOMMENDED. Manufacturer’s designation of the iEEG amplifier model.
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

%% information for the coordsystem.json file for MEG, EEG and iEEG
cfg.coordsystem.MEGCoordinateSystem                             = ft_getopt(cfg.coordsystem, 'MEGCoordinateSystem'                            ); % REQUIRED. Defines the coordinate system for the MEG sensors. See Appendix VIII: preferred names of Coordinate systems. If "Other", provide definition of the coordinate system in [MEGCoordinateSystemDescription].
cfg.coordsystem.MEGCoordinateUnits                              = ft_getopt(cfg.coordsystem, 'MEGCoordinateUnits'                             ); % REQUIRED. Units of the coordinates of MEGCoordinateSystem. MUST be ???m???, ???cm???, or ???mm???.
cfg.coordsystem.MEGCoordinateSystemDescription                  = ft_getopt(cfg.coordsystem, 'MEGCoordinateSystemDescription'                 ); % OPTIONAL. Freeform text description or link to document describing the MEG coordinate system system in detail.
cfg.coordsystem.EEGCoordinateSystem                             = ft_getopt(cfg.coordsystem, 'EEGCoordinateSystem'                            ); % OPTIONAL. Describes how the coordinates of the EEG sensors are to be interpreted.
cfg.coordsystem.EEGCoordinateUnits                              = ft_getopt(cfg.coordsystem, 'EEGCoordinateUnits'                             ); % OPTIONAL. Units of the coordinates of EEGCoordinateSystem. MUST be ???m???, ???cm???, or ???mm???.
cfg.coordsystem.EEGCoordinateSystemDescription                  = ft_getopt(cfg.coordsystem, 'EEGCoordinateSystemDescription'                 ); % OPTIONAL. Freeform text description or link to document describing the EEG coordinate system system in detail.
cfg.coordsystem.HeadCoilCoordinates                             = ft_getopt(cfg.coordsystem, 'HeadCoilCoordinateSystem'                       ); % OPTIONAL. Key:value pairs describing head localization coil labels and their coordinates, interpreted following the HeadCoilCoordinateSystem, e.g., {"NAS": [12.7,21.3,13.9], "LPA": [5.2,11.3,9.6], "RPA": [20.2,11.3,9.1]}. Note that coils are not always placed at locations that have a known anatomical name (e.g. for Neuromag/Elekta, Yokogawa systems); in that case generic labels can be used (e.g. {"coil1": [122,213,123], "coil2": [67,123,86], "coil3": [219,110,81]} ).
cfg.coordsystem.HeadCoilCoordinateSystem                        = ft_getopt(cfg.coordsystem, 'HeadCoilCoordinates'                            ); % OPTIONAL. Defines the coordinate system for the coils. See Appendix VIII: preferred names of Coordinate systems. If "Other", provide definition of the coordinate system in HeadCoilCoordinateSystemDescription.
cfg.coordsystem.HeadCoilCoordinateUnits                         = ft_getopt(cfg.coordsystem, 'HeadCoilCoordinateUnits'                        ); % OPTIONAL. Units of the coordinates of HeadCoilCoordinateSystem. MUST be ???m???, ???cm???, or ???mm???.
cfg.coordsystem.HeadCoilCoordinateSystemDescription             = ft_getopt(cfg.coordsystem, 'HeadCoilCoordinateSystemDescription'            ); % OPTIONAL. Freeform text description or link to document describing the Head Coil coordinate system system in detail.
cfg.coordsystem.DigitizedHeadPoints                             = ft_getopt(cfg.coordsystem, 'DigitizedHeadPoints'                            ); % OPTIONAL. Relative path to the file containing the locations of digitized head points collected during the session (e.g., "sub-01_headshape.pos"). RECOMMENDED for all MEG systems, especially for CTF and 4D/BTi. For Neuromag/Elekta the head points will be stored in the fif file.
cfg.coordsystem.DigitizedHeadPointsCoordinateSystem             = ft_getopt(cfg.coordsystem, 'DigitizedHeadPointsCoordinateSystem'            ); % OPTIONAL. Defines the coordinate system for the digitized head points. See Appendix VIII: preferred names of Coordinate systems. If "Other", provide definition of the coordinate system in DigitizedHeadPointsCoordinateSystemDescription.
cfg.coordsystem.DigitizedHeadPointsCoordinateUnits              = ft_getopt(cfg.coordsystem, 'DigitizedHeadPointsCoordinateUnits'             ); % OPTIONAL. Units of the coordinates of DigitizedHeadPointsCoordinateSystem. MUST be ???m???, ???cm???, or ???mm???.
cfg.coordsystem.DigitizedHeadPointsCoordinateSystemDescription  = ft_getopt(cfg.coordsystem, 'DigitizedHeadPointsCoordinateSystemDescription' ); % OPTIONAL. Freeform text description or link to document describing the Digitized head Points coordinate system system in detail.
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

% the task is both part of the file name (cfg.task) and one of the general JSON metadata fields (cfg.TaskName)
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
  if ~isequal(cfg.dataset, cfg.outputfile)
    cfg.method = 'convert';
  else
    cfg.method = 'decorate';
  end
  ft_notice('using cfg.method=''%s''', cfg.method);
end

% do some sanity checks on the input and the method
if istrue(cfg.mri.deface) && ~strcmp(cfg.method, 'convert')
  ft_error('defacing only works in combination with cfg.method=''convert''');
elseif ft_nargin>1 && ~strcmp(cfg.method, 'convert')
  ft_error('input data only works in combination with cfg.method=''convert''');
end

% do some more sanity checks on the input and the method
switch cfg.method
  case 'decorate'
    if ~isempty(cfg.outputfile) && ~isequal(cfg.dataset, cfg.outputfile)
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
else
  % data should be read from disk
  cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');
  if isfield(cfg, 'headerfile')
    typ = ft_filetype(cfg.headerfile);
  else
    typ = ft_filetype(cfg.presentationfile);
  end
end

% determine the sidecar files that are required
need_mri_json         = false;
need_meg_json         = false;
need_eeg_json         = false;
need_ieeg_json        = false;
need_events_tsv       = false; % for behavioral experiments
need_electrodes_tsv   = false; % only needed when actually present as data.cfg or as cfg.elec

switch typ
  case {'nifti', 'nifti2', 'nifti_fsl'}
    mri = ft_read_mri(cfg.dataset);
    if ~isempty(cfg.mri.dicomfile)
      % read the header details from the matching DICOM file specified by the user
      dcm = dicominfo(cfg.mri.dicomfile);
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
    if ~isempty(cfg.mri.dicomfile)
      % read the header details from the dicom matching file that was specified by the user
      dcm = dicominfo(cfg.mri.dicomfile);
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
    if isempty(cfg.trigger.event)
      trigger = ft_read_event(cfg.datafile, 'header', hdr);
    else
      % use the triggers as specified in the cfg
      trigger = cfg.trigger.event;
    end
    need_meg_json = true;

  case {'brainvision_vhdr', 'edf', 'eeglab_set'}
    % it is EEG data from disk and in a supported format
    hdr = ft_read_header(cfg.headerfile, 'checkmaxfilter', false);
    if strcmp(cfg.method, 'convert')
      % the data should be converted and written to disk
      dat = ft_read_data(cfg.datafile, 'header', hdr, 'checkboundary', false, 'begsample', 1, 'endsample', hdr.nSamples*hdr.nTrials);
    end
    if isempty(cfg.trigger.event)
      trigger = ft_read_event(cfg.datafile, 'header', hdr);
    else
      % use the triggers as specified in the cfg
      trigger = cfg.trigger.event;
    end
    need_eeg_json = true;

  case 'raw'
    % the input data structure contains raw electrophysiology data
    if isequal(cfg.datatype, 'meg')
      need_meg_json = true;
    elseif isequal(cfg.datatype, 'eeg')
      need_eeg_json = true;
    elseif isequal(cfg.datatype, 'ieeg')
      need_ieeg_json = true;
    else
      ft_warning('assuming that the dataset represents EEG');
      need_eeg_json = true;
    end

    hdr = ft_fetch_header(varargin{1});
    if strcmp(cfg.method, 'convert')
      % the data should be written to disk
      dat = ft_fetch_data(varargin{1}, 'checkboundary', false, 'begsample', 1, 'endsample', hdr.nSamples*hdr.nTrials);
    end
    if isempty(cfg.trigger.event)
      trigger = ft_fetch_event(varargin{1});
    else
      % use the triggers as specified in the cfg
      trigger = cfg.trigger.event;
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

  case 'presentation_log'
    need_events_tsv = true;

  otherwise
    % the file on disk contains raw electrophysiology data
    if isequal(cfg.datatype, 'meg')
      need_meg_json = true;
    elseif isequal(cfg.datatype, 'eeg')
      need_eeg_json = true;
    elseif isequal(cfg.datatype, 'ieeg')
      need_ieeg_json = true;
    else
      ft_warning('assuming that the dataset represents EEG');
      need_eeg_json = true;
    end

    hdr = ft_read_header(cfg.headerfile, 'checkmaxfilter', false);
    if strcmp(cfg.method, 'convert')
      % the data should be converted and written to disk
      dat = ft_read_data(cfg.datafile, 'header', hdr, 'checkboundary', false, 'begsample', 1, 'endsample', hdr.nSamples*hdr.nTrials);
    end
    if isempty(cfg.trigger.event)
      trigger = ft_read_event(cfg.datafile, 'header', hdr);
    else
      % use the triggers as specified in the cfg
      trigger = cfg.trigger.event;
    end

end % switch typ

need_events_tsv       = need_events_tsv || need_meg_json || need_eeg_json || need_ieeg_json || (need_mri_json && (contains(cfg.outputfile, 'task') || ~isempty(cfg.TaskName) || ~isempty(cfg.task)));
need_channels_tsv     = need_meg_json || need_eeg_json || need_ieeg_json;
need_coordsystem_json = need_meg_json; % FIXME this is also needed when EEG and iEEG electrodes are present

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get the defaults and user-specified settings for each possible sidecar file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start with empty metadata descriptions
mri_json         = [];
meg_json         = [];
eeg_json         = [];
ieeg_json        = [];
events_tsv       = [];
channels_tsv     = [];
electrodes_tsv   = [];
coordsystem_json = [];

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

% make the relevant selection, most json fields start with a capital letter, some start with iEEG
fn = fieldnames(cfg.ieeg);
fn = fn(~cellfun(@isempty, regexp(fn, '^[A-Z].*|^iEEG')));
ieeg_settings = keepfields(cfg.ieeg, fn);

% make the relevant selection, all json fields start with a capital letter
fn = fieldnames(cfg.coordsystem);
fn = fn(~cellfun(@isempty, regexp(fn, '^[A-Z].*')));
coordsystem_settings = keepfields(cfg.coordsystem, fn);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% construct the content for the json and tsv files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% need_mri_json
if need_mri_json
  % start with the information from the DICOM header
  mri_json = keepfields(dcm, fieldnames(mri_settings));
  % merge the information specified by the user with that from the data
  % in case fields appear in both, the first input overrules the second
  mri_json = mergeconfig(mri_json, mri_settings, false);
  mri_json = mergeconfig(mri_json, generic_settings, false);
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
  meg_json.MiscChannelCount           = sum(strcmp(hdr.chantype, 'misc'));
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
    meg_json.Manufacturer             = 'Neuromag/Elekta';
    meg_json.ManufacturersModelName   = 'Neuromag-122';
  elseif ft_senstype(hdr.grad, 'neuromag306')
    meg_json.Manufacturer             = 'Neuromag/Elekta';
    % the ManufacturersModelName could be either Vectorview or Triux
  end

  % merge the information specified by the user with that from the data
  % in case fields appear in both, the first input overrules the second
  meg_json = mergeconfig(meg_json, meg_settings, false);
  meg_json = mergeconfig(meg_json, generic_settings, false);
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
  eeg_json = mergeconfig(eeg_json, eeg_settings, false);
  eeg_json = mergeconfig(eeg_json, generic_settings, false);
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
  ieeg_json = mergeconfig(ieeg_json, ieeg_settings, false);
  ieeg_json = mergeconfig(ieeg_json, generic_settings, false);
end

%% need_coordsystem_json
if need_coordsystem_json
  if ft_senstype(hdr.grad, 'ctf')
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
    % merge the information specified by the user with that from the data
    % in case fields appear in both, the first input overrules the second
    coordsystem_json = mergeconfig(coordsystem_settings, coordsystem_json, false); % FIXME the order of precedence is different here
  else
    ft_warning('coordsystem handling not yet supported for %s', ft_senstype(hdr.grad));
  end
end % if need_coordsystem_json

%% need_channels_tsv
if need_channels_tsv
  % start with an empty table
  channels_tsv = table();

  % ensure that all columns have the correct number of elements
  fn = setdiff(fieldnames(cfg.channels), 'writesidecar');
  for i=1:numel(fn)
    if ischar(cfg.channels.(fn{i}))
      cfg.channels.(fn{i}) = repmat({cfg.channels.(fn{i})}, hdr.nChans, 1);
    elseif isnumeric(cfg.channels.(fn{i})) && numel(cfg.channels.(fn{i}))==1
      cfg.channels.(fn{i}) = repmat(cfg.channels.(fn{i}), hdr.nChans, 1);
    elseif iscell(cfg.channels.(fn{i})) && numel(cfg.channels.(fn{i}))==1
      cfg.channels.(fn{i}) = repmat(cfg.channels.(fn{i}), hdr.nChans, 1);
    elseif iscell(cfg.channels.(fn{i})) && numel(cfg.channels.(fn{i}))==hdr.nChans
      cfg.channels.(fn{i}) = cfg.channels.(fn{i})(:); % just make sure it is a column vector
    else
      ft_error('incorrect specification of cfg.channels.%s', fn{i});
    end
  end

  % these columns can be determined from the header
  channels_tsv.name                = merge_vector(hdr.label(:),    cfg.channels.name);
  channels_tsv.type                = merge_vector(hdr.chantype(:), cfg.channels.type);
  channels_tsv.units               = merge_vector(hdr.chanunit(:), cfg.channels.units);
  channels_tsv.sampling_frequency  = merge_vector(repmat(hdr.Fs, hdr.nChans, 1), cfg.channels.sampling_frequency);

  % all other columns have to be specified by the user
  fn = setdiff(fieldnames(cfg.channels), {'writesidecar', 'name', 'type', 'unit', 'sampling_frequency'});
  for i=1:numel(fn)
    if isnumeric(cfg.channels.(fn{i})) && all(isnan(cfg.channels.(fn{i})))
      % this is the default when not specified by the user, do not add it to the table
    else
      channels_tsv.(fn{i}) = cfg.channels.(fn{i});
    end
  end
end % if need_channels_tsv

%% need_electrodes_tsv
if need_electrodes_tsv
  % start with an empty table
  electrodes_tsv = table();

  % ensure that all columns have the correct number of elements
  fn = setdiff(fieldnames(cfg.electrodes), 'writesidecar');
  for i=1:numel(fn)
    if ischar(cfg.electrodes.(fn{i}))
      cfg.electrodes.(fn{i}) = repmat({cfg.electrodes.(fn{i})}, hdr.nChans, 1);
    elseif isnumeric(cfg.electrodes.(fn{i})) && numel(cfg.electrodes.(fn{i}))==1
      cfg.electrodes.(fn{i}) = repmat(cfg.electrodes.(fn{i}), hdr.nChans, 1);
    elseif iscell(cfg.electrodes.(fn{i})) && numel(cfg.electrodes.(fn{i}))==1
      cfg.electrodes.(fn{i}) = repmat(cfg.electrodes.(fn{i}), hdr.nChans, 1);
    elseif iscell(cfg.electrodes.(fn{i})) && numel(cfg.electrodes.(fn{i}))==hdr.nChans
      cfg.electrodes.(fn{i}) = cfg.electrodes.(fn{i})(:); % just make sure it is a column vector
    else
      ft_error('incorrect specification of cfg.electrodes.%s', fn{i});
    end
  end

  % these columns can be determined from the header
  electrodes_tsv.name   = merge_vector(elec.label,        cfg.electrodes.name);
  electrodes_tsv.x      = merge_vector(elec.elecpos(:,1), cfg.electrodes.x);
  electrodes_tsv.y      = merge_vector(elec.elecpos(:,2), cfg.electrodes.y);
  electrodes_tsv.z      = merge_vector(elec.elecpos(:,3), cfg.electrodes.z);

  % all other columns have to be specified by the user
  fn = setdiff(fieldnames(cfg.electrodes), {'writesidecar', 'name', 'x', 'y', 'z'});
  for i=1:numel(fn)
    if isnumeric(cfg.electrodes.(fn{i})) && all(isnan(cfg.electrodes.(fn{i})))
      % this is the default when not specified by the user, do not add it to the table
    else
      electrodes_tsv.(fn{i}) = cfg.electrodes.(fn{i});
    end
  end
end % need_electrodes_tsv

%% need_events_tsv
if need_events_tsv
  % read the presentation file that may accompany the functional data
  if ~isempty(cfg.presentationfile)
    presentation = ft_read_event(cfg.presentationfile);
  else
    presentation = [];
  end

  if need_mri_json

    if isempty(presentation)
      ft_warning('cfg.presentationfile not specified, cannot determine events')

      onset    = [];
      duration = [];
      events_tsv = table(onset, duration);

    else
      % align the events from the presentation log file with the MR volumes
      % this requires one event per volume in the presentation file

      % merge the information with the json sidecar file
      % in case fields appear in both, the first input overrules the second
      tmp = mergeconfig(mri_json, read_json(corresponding_json(cfg.outputfile)), false);
      assert(~isempty(tmp.RepetitionTime), 'you must specify cfg.mri.RepetitionTime');

      % create a header structure that represents the fMRI timeseries
      hdr.Fs = 1/tmp.RepetitionTime;
      hdr.nSamples = mri.dim(4);

      % create a event structure with one event for each fMRI volume
      volume = [];
      for i=1:hdr.nSamples
        volume(i).type   = 'volume';
        volume(i).sample = i;
      end

      % find the presentation events corresponding to each volume
      selpres = select_event(presentation, cfg.presentation.eventtype, cfg.presentation.eventvalue);
      selpres = presentation(selpres);

      ft_info('%d volumes, %d presentation events', length(volume), length(selpres));
      if length(volume)>length(selpres)
        % this happens when the scanner keeps running while presentation has already been stopped
        n = length(volume)-length(selpres);
        ft_warning('discarding last %d volumes for realignment of events', n);
        volume = volume(1:end-n);
      elseif length(selpres)>length(volume)
        % this happens when DICOM volumes that are represented in the
        % presentation log have been deleted from disk
        switch cfg.presentation.skip
          case 'first'
            ft_warning('discarding first %d presentation events for realignment of events', n);
            selpres = selpres((n+1):end);
          case 'last'
            ft_warning('discarding last %d presentation events for realignment of events', n);
            selpres = selpres(1:end-n);
          case 'none'
            ft_error('not enough volumes to match the presentation events');
        end % case
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % the following code is largely shared between the MEG and MRI section
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % predict the sample number from the timestamp
      model     = polyfit([selpres.timestamp], [volume.sample], 1);
      estimated = polyval(model, [selpres.timestamp]);

      if istrue(cfg.feedback)
        [p, f, x] = fileparts(cfg.dataset);
        figure('name', ['PRESENTATION - ' f]);
        subplot(2,1,1)
        hold on
        % presentation timestamps are expressed in units of 0.1 miliseconds
        plot([selpres.timestamp]/1e4, [volume.sample], 'b.')
        plot([selpres.timestamp]/1e4, estimated, 'ro')
        xlabel('presentation time (s)')
        ylabel('MRI volumes')
        legend({'observed', 'predicted'})

        subplot(2,1,2)
        plot([selpres.timestamp]/1e4, ([volume.sample]-estimated)/hdr.Fs, 'g.')
        xlabel('presentation time (s)')
        ylabel('difference (s)')
      end

      % estimate the time in seconds of all presentation events
      estimated = polyval(model, [presentation.timestamp]);
      estimated = round(1000*estimated)/1000; % round to three decimals
      for i=1:numel(estimated)
        presentation(i).sample = estimated(i);
      end
      % convert the event structure to a TSV table
      presentation_tsv = event2table(hdr, presentation);
      % rename the column to "volume" instead of "sample"
      sel = strcmp(presentation_tsv.Properties.VariableNames, 'sample');
      presentation_tsv.Properties.VariableNames{sel} = 'volume';

      % for fMRI the presentation log file is the only source of events
      events_tsv = presentation_tsv;
      clear presentation_tsv selpres volume

      % sort ascending on the onset of each event
      events_tsv = sortrows(events_tsv, 'onset');
    end

  elseif need_meg_json || need_eeg_json || need_ieeg_json
    % merge the events from the trigger channel with those from the (optional) presentation file

    if istable(cfg.events.trl)
      % check that the column names are valid
      assert(stcmp(cfg.events.trl.Properties.VariableNames{1}, 'begsample'));
      assert(stcmp(cfg.events.trl.Properties.VariableNames{2}, 'endsample'));
      assert(stcmp(cfg.events.trl.Properties.VariableNames{3}, 'offset'));
      % use the events table as it is
      events_tsv = cfg.events.trl;
    elseif ~isempty(cfg.events.trl)
      % convert the trl matrix to an events table
      begsample = cfg.events.trl(:,1);
      endsample = cfg.events.trl(:,2);
      offset    = cfg.events.trl(:,3); % FIXME this is not used
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
    else
      % convert the events from the dataset into a table
      events_tsv = event2table(hdr, trigger);
    end

    if ~isempty(presentation) && ~isempty(trigger)
      % align the events from the presentation log file with the triggers

      % select the correspopnding triggers and events in the presentation file
      seltrig = select_event(trigger,      cfg.trigger.eventtype,      cfg.trigger.eventvalue);
      selpres = select_event(presentation, cfg.presentation.eventtype, cfg.presentation.eventvalue);
      seltrig = trigger(seltrig);
      selpres = presentation(selpres);

      %if length(seltrig)~=length(selpres)
      %  ft_error('inconsistent number: %d triggers, %d presentation events', length(seltrig), length(selpres));
      %end

      ft_info('%d triggers, %d presentation events', length(seltrig), length(selpres));
      if length(seltrig)>length(selpres)
        % don't know how to solve this
        ft_error('inconsistent number: %d triggers, %d presentation events', length(seltrig), length(selpres));
      elseif length(selpres)>length(seltrig)
        n = length(selpres)-length(seltrig);
        % This could happen when due to acquisition problems there is more than one
        % *.ds directory. If this is a known case, cfg.presentation.skip can be used.
        % Note that this only works, if there are two ds-datasets (not more).
        switch cfg.presentation.skip
          case 'first'
            ft_warning('discarding first %d presentation events for realignment of events', n);
            selpres = selpres((n+1):end);
          case 'last'
            ft_warning('discarding last %d presentation events for realignment of events', n);
            selpres = selpres(1:end-n);
          case 'none'
            ft_error('not enough triggers to match the presentation events');
        end % case
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % the following code is largely shared between the MEG and MRI section
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % predict the presentation sample number from the presentation timestamp
      model     = polyfit([selpres.timestamp], [seltrig.sample], 1);
      estimated = polyval(model, [selpres.timestamp]);

      if istrue(cfg.feedback)
        [~, f, ~] = fileparts(cfg.dataset);
        figure('name', ['PRESENTATION - ' f]);
        subplot(2,1,1)
        hold on
        % presentation timestamps are expressed in units of 0.1 miliseconds
        plot([selpres.timestamp]/1e4, [seltrig.sample], 'b.')
        plot([selpres.timestamp]/1e4, estimated, 'ro')
        xlabel('presentation time (s)')
        ylabel('data samples')
        legend({'observed', 'predicted'})

        subplot(2,1,2)
        plot([selpres.timestamp]/1e4, ([seltrig.sample]-estimated)/hdr.Fs, 'g.')
        xlabel('presentation time (s)')
        ylabel('difference (s)')
      end

      % estimate the sample number and time in seconds of all presentation events
      estimated = polyval(model, [presentation.timestamp]);
      estimated = round(estimated); % round to the nearest sample
      for i=1:numel(estimated)
        presentation(i).sample = estimated(i);
      end
      % convert the event structure to a TSV table
      presentation_tsv = event2table(hdr, presentation);

      % the events from the the presentation log file should be merged with the triggers
      % trigger values are often numeric, whereas presentation event values are often strings
      if isnumeric(events_tsv.value) && ~isnumeric(presentation_tsv.value)
        % convert them, otherwise the concatenation fails
        events_tsv.value = num2cell(events_tsv.value);
      end
      % concatenate them
      events_tsv = [events_tsv; presentation_tsv];

      clear presentation_tsv selpres seltrig
    end

  elseif need_events_tsv
    % convert the presentation structure to a TSV table
    events_tsv = struct2table(presentation);
    % the sample and offset are specific for events that are aligned with samples in the data
    events_tsv = removevars(events_tsv, 'sample');
    events_tsv = removevars(events_tsv, 'offset');
    events_tsv.onset = events_tsv.timestamp/10000; % the resolution is 0.1 milliseconds
    events_tsv = removevars(events_tsv, 'timestamp');
    % a value of zero implies that the delta function or event is so short as to be effectively modeled as an impulse.
    events_tsv.duration = zeros(size(events_tsv.duration));
    % reorder them so that onset and duration are the first two columns
    order = [];
    order = [order find(strcmp(events_tsv.Properties.VariableNames, 'onset'))];
    order = [order find(strcmp(events_tsv.Properties.VariableNames, 'duration'))];
    order = [order setdiff(1:numel(events_tsv.Properties.VariableNames), order)];
    events_tsv = events_tsv(:,order);

  end

  if ~isempty(events_tsv) && isfield(events_tsv, 'offset')
    % sort the events ascending on the onset
    events_tsv = sortrows(events_tsv, 'onset');
  end

end % if need_events_tsv

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% write the data to the output file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch cfg.method
  case 'decorate'
    % there is nothing to do here

  case 'convert'
    % the output depends on the type of input data
    switch typ
      case {'nifti', 'nifti2', 'nifti_fsl', 'dicom', 'volume'}
        [p, f, x] = fileparts(cfg.outputfile);
        if ~isequal(x, '.nii')
          cfg.outputfile = fullfile(p, [f '.nii']);
        end
        if istrue(cfg.mri.deface)
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

      case {'presentation_log'}
        % do not write data, but only write the events.tsv file

      case {'ctf_ds', 'ctf_meg4', 'ctf_res4', 'ctf151', 'ctf275', 'neuromag_fif', 'neuromag122', 'neuromag306'}
        ft_error('please use a system specific tool for converting MEG datasets');

      otherwise
        [p, f, x] = fileparts(cfg.outputfile);
        if ~isequal(x, '.vhdr')
          cfg.outputfile = fullfile(p, [f '.vhdr']);
        end
        ft_info('writing %s\n', cfg.outputfile);
        ft_write_data(cfg.outputfile, dat, 'dataformat', 'brainvision_eeg', 'header', hdr, 'event', trigger);
    end % switch typ

  case 'copy'
    [~, ~, xin] = fileparts(cfg.dataset);
    [p, ~, xout] = fileparts(cfg.outputfile);
    if ~strcmp(xin, xout)
      ft_error('input and output filename extension do not match');
    end

    switch typ
      case {'ctf_ds', 'ctf_meg4', 'ctf_res4', 'ctf151', 'ctf275'}
        % the data consists of a directory with multiple files inside
        ft_info('copying %s to %s\n', cfg.dataset, cfg.outputfile);
        isdir_or_mkdir(p);
        copy_ctf_files(cfg.dataset, cfg.outputfile, false);

      case {'brainvision_vhdr', 'brainvision_vmrk', 'brainvision_eeg', 'brainvision_dat', 'brainvision_seg'}
        % the data consists of three files and the header file contains pointers to the markers and data
        ft_info('copying %s to %s\n', cfg.dataset, cfg.outputfile);
        isdir_or_mkdir(p);
        copy_brainvision_files(cfg.dataset, cfg.outputfile, false);

      otherwise
        ft_info('copying %s to %s\n', cfg.dataset, cfg.outputfile);
        isdir_or_mkdir(p);
        copyfile(cfg.dataset, cfg.outputfile);
    end

  otherwise
    ft_error('unsupported value for cfg.method')
end % switch method


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% write the metadata to the json and tsv files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remove fields that have an empty value
mri_json  = remove_empty(mri_json);
meg_json  = remove_empty(meg_json);
eeg_json  = remove_empty(eeg_json);
ieeg_json = remove_empty(ieeg_json);

if ~isempty(mri_json)
  filename = corresponding_json(cfg.outputfile);ls
  if isfile(filename)
    existing = read_json(filename);
  else
    existing = [];
  end
  switch cfg.mri.writesidecar
    case 'yes'
      if ~isempty(existing)
        ft_warning('not overwriting the existing and non-empty file ''%s''', filename);
      else
        write_json(filename, mri_json);
      end
    case 'replace'
      write_json(filename, mri_json);
    case 'merge'
      write_json(filename, mergeconfig(mri_json, existing, false))
    case 'no'
      % do nothing
    otherwise
      ft_error('incorrect option for cfg.mri.writesidecar');
  end % switch
end

if ~isempty(meg_json)
  filename = corresponding_json(cfg.outputfile);
  if isfile(filename)
    existing = read_json(filename);
  else
    existing = [];
  end
  switch cfg.meg.writesidecar
    case 'yes'
      if ~isempty(existing)
        ft_warning('not overwriting the existing and non-empty file ''%s''', filename);
      else
        write_json(filename, meg_json);
      end
    case 'replace'
      write_json(filename, meg_json);
    case 'merge'
      write_json(filename, mergeconfig(meg_json, existing, false));
    case 'no'
      % do nothing
    otherwise
      ft_error('incorrect option for cfg.meg.writesidecar');
  end % switch
end

if ~isempty(eeg_json)
  filename = corresponding_json(cfg.outputfile);
  if isfile(filename)
    existing = read_json(filename);
  else
    existing = [];
  end
  switch cfg.eeg.writesidecar
    case 'yes'
      if ~isempty(existing)
        ft_warning('not overwriting the existing and non-empty file ''%s''', filename);
      else
        write_json(filename, eeg_json);
      end
    case 'replace'
      write_json(filename, eeg_json);
    case 'merge'
      write_json(filename, mergeconfig(eeg_json, existing, false));
    case 'no'
      % do nothing
    otherwise
      ft_error('incorrect option for cfg.eeg.writesidecar');
  end % switch
end

if ~isempty(ieeg_json)
  filename = corresponding_json(cfg.outputfile);
  if isfile(filename)
    existing = read_json(filename);
  else
    existing = [];
  end
  switch cfg.ieeg.writesidecar
    case 'yes'
      if ~isempty(existing)
        ft_warning('not overwriting the existing and non-empty file ''%s''', filename);
      else
        write_json(filename, ieeg_json);
      end
    case 'replace'
      write_json(filename, ieeg_json);
    case 'merge'
      write_json(filename, mergeconfig(ieeg_json, existing, false));
    case 'no'
      % do nothing
    otherwise
      ft_error('incorrect option for cfg.ieeg.writesidecar');
  end % switch
end

if ~isempty(coordsystem_json)
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
  f = remove_datatype(f);           % remove _bold, _meg, etc.
  filename = fullfile(p, [f '_coordsystem.json']);
  if isfile(filename)
    existing = read_json(filename);
  else
    existing = [];
  end
  switch cfg.coordsystem.writesidecar
    case 'yes'
      if ~isempty(existing)
        ft_warning('not overwriting the existing and non-empty file ''%s''', filename);
      else
        write_json(filename, coordsystem_json);
      end
    case 'replace'
      write_json(filename, coordsystem_json);
    case 'merge'
      write_json(filename, mergeconfig(coordsystem_json, existing, false))
    case 'no'
      % do nothing
    otherwise
      ft_error('incorrect option for cfg.coordsystem.writesidecar');
  end % switch
end

if ~isempty(channels_tsv)
  [p, f] = fileparts(cfg.outputfile);
  f = remove_datatype(f); % remove _bold, _meg, etc.
  filename = fullfile(p, [f '_channels.tsv']);
  if isfile(filename)
    existing = read_tsv(filename);
  else
    existing = [];
  end
  switch cfg.channels.writesidecar
    case 'yes'
      if ~isempty(existing)
        ft_warning('not overwriting the existing and non-empty file ''%s''', filename);
      else
        write_tsv(filename, channels_tsv);
      end
    case 'replace'
      write_tsv(filename, channels_tsv);
    case 'merge'
      if ~isempty(existing)
        channels_tsv = unique(cat(1, channels_tsv, existing));
      end
      write_tsv(filename, channels_tsv);
    case 'no'
      % do nothing
    otherwise
      ft_error('incorrect option for cfg.channels.writesidecar');
  end % switch
end

if ~isempty(electrodes_tsv)
  [p, f] = fileparts(cfg.outputfile);
  f = remove_datatype(f); % remove _bold, _meg, etc.
  filename = fullfile(p, [f '_electrodes.tsv']);
  if isfile(filename)
    existing = read_tsv(filename);
  else
    existing = [];
  end
  switch cfg.electrodes.writesidecar
    case 'yes'
      if ~isempty(existing)
        ft_warning('not overwriting the existing and non-empty file ''%s''', filename);
      else
        write_tsv(filename, electrodes_tsv);
      end
    case 'replace'
      write_tsv(filename, electrodes_tsv);
    case 'merge'
      if ~isempty(existing)
        electrodes_tsv = unique(cat(1, electrodes_tsv, existing));
      end
      write_tsv(filename, electrodes_tsv);
    case 'no'
      % do nothing
    otherwise
      ft_error('incorrect option for cfg.electrodes.writesidecar');
  end % switch
end

if ~isempty(events_tsv)
  [p, f] = fileparts(cfg.outputfile);
  isdir_or_mkdir(p);
  f = remove_datatype(f); % remove _bold, _meg, etc.
  filename = fullfile(p, [f '_events.tsv']);
  if isfile(filename)
    existing = read_tsv(filename);
  else
    existing = [];
  end
  switch cfg.events.writesidecar
    case 'yes'
      if ~isempty(existing)
        ft_warning('not overwriting the existing and non-empty file ''%s''', filename);
      else
        write_tsv(filename, events_tsv);
      end
    case 'replace'
      write_tsv(filename, events_tsv);
    case 'merge'
      if ~isempty(existing)
        events_tsv = unique(cat(1, events_tsv, existing));
      end
      write_tsv(filename, events_tsv);
    case 'no'
      % do nothing
    otherwise
      ft_error('incorrect option for cfg.events.writesidecar');
  end % switch
end

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
  switch cfg.dataset_description.writesidecar
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
      ft_error('incorrect option for cfg.dataset_description.writesidecar');
  end % switch

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % update the participants.tsv
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  filename = fullfile(cfg.bidsroot, 'participants.tsv');

  this = table();
  this.participant_id = ['sub-' cfg.sub];
  fn = fieldnames(cfg.participants);
  for i=1:numel(fn)
    % write [] as 'n/a'
    % write nan as 'n/a'
    % write boolean as 'True' or 'False'
    this.(fn{i}) = output_compatible(cfg.participants.(fn{i}));
  end

  if isfile(filename)
    participants = read_tsv(filename);
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
  this.filename = fullfile(cfg.datatype, [f x]);
  fn = fieldnames(cfg.scans);
  for i=1:numel(fn)
    % write [] as 'n/a'
    % write nan as 'n/a'
    % write boolean as 'True' or 'False'
    this.(fn{i}) = output_compatible(cfg.scans.(fn{i}));
  end

  if isfile(filename)
    scans = read_tsv(filename);
    scans = merge_table(scans, this, 'filename');
  else
    scans = this;
  end

  % write the updated file back to disk
  write_tsv(filename, scans);

end

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
typ = {'FLAIR', 'FLASH', 'PD', 'PDT2', 'PDmap', 'T1map', 'T1rho', 'T1w', 'T2map', 'T2star', 'T2w', 'angio', 'bold', 'bval', 'bvec', 'channels', 'coordsystem', 'defacemask', 'dwi', 'eeg', 'epi', 'events', 'fieldmap', 'headshape', 'ieeg', 'inplaneT1', 'inplaneT2', 'magnitude', 'magnitude1', 'magnitude2', 'meg', 'phase1', 'phase2', 'phasediff', 'photo', 'physio', 'sbref', 'stim'};
for i=1:numel(typ)
  if endsWith(f, ['_' typ{i}])
    f = f(1:end-length(typ{i})-1); % also the '_'
    return
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION convert event structure into table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tab = event2table(hdr, trigger)
if isempty(trigger)
  tab = [];
else
  for i=1:numel(trigger)
    % normalize missing values
    if isempty(trigger(i).sample)
      trigger(i).sample = NaN;
    end
    if isempty(trigger(i).duration)
      trigger(i).duration = NaN;
    end
    if isempty(trigger(i).value)
      trigger(i).value = NaN;
    end
  end
  onset        = (([trigger.sample]-1)/hdr.Fs)';   % in seconds
  duration     = ([trigger.duration]/hdr.Fs)';     % in seconds
  sample       = ([trigger.sample])';              % in samples, the first sample of the file is 1
  type         = {trigger.type}';
  value        = {trigger.value}';
  if all(cellfun(@isnumeric, type))
    % this can be an array of strings or values
    type = cell2mat(type);
  end
  if all(cellfun(@isnumeric, value))
    % this can be an array of strings or values
    value = cell2mat(value);
  end
  tab = table(onset, duration, sample, type, value);
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
  %   val = [];
  %   if isnumeric(val) && any(isnan(val))
  %     sel = isnan(val);
  %     val = num2cell(val);
  %     val(sel) = {'n/a'};
  %     % replace the numeric value NaN with the string 'n/a'
  %     tsv.(fn{i}) = val;
  %   elseif iscell(val)
  %     for j=1:numel(val)
  %       if isnan(val{j})
  %         val{j} = 'n/a';
  %       end
  %     end
  %     % replace the numeric value NaN with the string 'n/a'
  %     tsv.(fn{i}) = val;
  %   end
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
json = remove_empty(json);
json = sort_fields(json);
json = ft_struct2char(json); % convert strings into char-arrays
ft_info('writing %s\n', filename);
ft_hastoolbox('jsonlab', 1);
% see also the output_compatible helper function
% write nan as 'n/a'
% write boolean as True/False
str = savejson('', json, 'NaN', '"n/a"', 'ParseLogical', true);
fid = fopen(filename, 'w');
fwrite(fid, str);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filename = corresponding_json(filename)
[p, f, x] = fileparts(filename);
if isequal(x, '.gz') && endsWith(f, '.nii')
  % it is a gzip compressed nifti file, remove the .nii from the file name
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
switch typ
  case {'T1w' 'T2w' 'T1rho' 'T1map' 'T2map' 'T2star' 'FLAIR' 'FLASH' 'PD' 'PDmap' 'PDT2' 'inplaneT1' 'inplaneT2' 'angio' 'defacemask'}
    dir = 'anat';
  case {'bold' 'sbref'} % this could also include 'events' 'physio' 'stim'
    dir = 'func';
  case {'dwi' 'bvec' 'bval'}
    dir = 'dwi';
  case {'phasediff' 'phase1' 'phase2' 'magnitude1' 'magnitude2' 'magnitude' 'fieldmap' 'epi'}
    dir = 'fmap';
  case {'events' 'stim' 'physio'}
    dir = 'beh';
  case {'meg'} % this could also include 'channels' 'photo' 'coordsystem' 'headshape'
    dir = 'meg';
  case {'eeg'} % this could also include 'channels' 'photo' 'coordsystem'
    dir = 'eeg';
  case {'ieeg'} % this could also include 'channels' 'photo' 'coordsystem'
    dir = 'ieeg';
  otherwise
    ft_error('unrecognized data type "%s"', typ);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = output_compatible(val)
if iscell(val)
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
