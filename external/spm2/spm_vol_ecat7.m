function V = spm_vol_ecat7(fname,required)
% Get header information etc. for ECAT 7 images.
%  FORMAT V = spm_vol_ecat7(fname,required)
%  P         - an ECAT 7 filename.
%  fname     - a structure containing image volume information.
%  required  - an optional text argument specifying which volumes to
%              use.  The default is '1010001'.  Also, the argument
%              'all' will result in all matrices being extracted.
%
% The elements of V are described by the help for spm_vol, except for
% additional fields (V.mh and V.sh) that contains the main- and sub-
% header information.
%
% _______________________________________________________________________
% @(#)spm_vol_ecat7.m	2.12 John Ashburner and Roger Gunn 03/05/23

V        = [];
if nargin==1,
	required = '1010001';
elseif ischar(required),
else,
	if isfinite(required),
		required = sprintf('%.7x',16842752+required);
	else,
		required = 'all';
	end;
end;

fp       = fopen(fname,'r','ieee-be');
if (fp == -1)
	fname  = [spm_str_manip(fname,'sr') '.v'];
	fp = fopen(fname,'r','ieee-be');
	if (fp == -1)
		return;
	end;
end

mh       = ECAT7_mheader(fp);
if ~strcmp(mh.MAGIC_NUMBER,'MATRIX70v') & ~strcmp(mh.MAGIC_NUMBER,'MATRIX71v') & ~strcmp(mh.MAGIC_NUMBER,'MATRIX72v'),
        fclose(fp);
        return; % Quietly return.
end

if (mh.FILE_TYPE ~= 7)
        error(['"' spm_str_manip(fname,'k20d') '" isn''t an image.']);
        fclose(fp);
        return;
end

list     = s7_matlist(fp);
if strcmp(required(1,:),'all'),
	matches  = find((list(:,4) == 1) | (list(:,4) == 2));
	llist = list(matches,:);
else,
	for i=1:size(required,1),
		matnum   = sscanf(required(i,:),'%x');
		matches  = find( (list(:,1) == matnum) & ((list(:,4) == 1) | (list(:,4) == 2)));
		if (size(matches,1) ~= 1)
			Error(['"' spm_str_manip(fname,'k20d') '" doesn''t have the required image.']);
			fclose(fp);
		end;
		llist(i,:) = list(matches,:);
	end;
end;

frame_num = sscanf(required,'%x')-16842752;
private   = struct('mh',[],'sh',[]);
V         = struct('fname','','dim',[],'pinfo',[],'mat',[], 'descrip','','n',[],'private',private);
V(:)      = [];

for i=1:size(llist,1),
	sh       = ECAT7_sheader(fp,llist(i,2));
	dim      = [sh.X_DIMENSION sh.Y_DIMENSION sh.Z_DIMENSION 4];
	if ~spm_platform('bigend') & dim(4)~=2, dim(4) = dim(4)*256; end;

	pinfo    = [sh.SCALE_FACTOR*mh.ECAT_CALIBRATION_FACTOR ; 0 ; 512*llist(i,2)];

	dircos   = diag([-1 -1 -1]);
	step     = ([sh.X_PIXEL_SIZE sh.Y_PIXEL_SIZE sh.Z_PIXEL_SIZE]*10);
	start    = -(dim(1:3)'/2).*step';
	mat      = [[dircos*diag(step) dircos*start] ; [0 0 0 1]];
	[pth,nam,ext] = fileparts(fname);
	matname  = fullfile(pth,[nam '.mat']);
	if exist(matname) == 2,
		str=load(matname);
		if isfield(str,'mat'),
			mat = str.mat;
		elseif isfield(str,'M'),
			mat = str.M;
		end;
	end;
	V(i).fname      = fname;
	V(i).dim        = dim;
	V(i).mat        = mat;
	V(i).pinfo      = pinfo;
	V(i).n          = llist(i,1)-16842752;
	V(i).descrip    = sh.ANNOTATION;
	V(i).private.mh = mh;
	V(i).private.sh = sh;
end;
fclose(fp);
return;
%_______________________________________________________________________

%_______________________________________________________________________
%S7_MATLIST List the available matrixes in an ECAT 7 file.
%	LIST = S7_MATLIST(FP) lists the available matrixes
%	in the file specified by FP.
%
%	Columns in LIST:
%		1 - Matrix identifier.
%		2 - Matrix subheader record number
%		3 - Last record number of matrix data block.
%		4 - Matrix status:
%			1 - exists - rw
%			2 - exists - ro
%			3 - matrix deleted
%
function list = s7_matlist(fp);

% I believe fp should be opened with:
%	fp = fopen(filename,'r','ieee-be');

fseek(fp,512,'bof');
block = fread(fp,128,'int');
if (size(block,1) ~= 128) list = []; return; end;
block = reshape(block,4,32);

list  = [];

while (block(2,1) ~= 2)
	if (block(1,1)+block(4,1) ~= 31)
		list = []; return;
	end
	list = [list block(:,2:32)];

	fseek(fp,512*(block(2,1)-1),'bof');
	block = fread(fp,128,'int');
	if (size(block,1) ~= 128) list = []; return; end;
	block = reshape(block,4,32);
end

list = [list block(:,2:(block(4,1)+1))];
list = list';
return;
%_______________________________________________________________________

%_______________________________________________________________________
function SHEADER=ECAT7_sheader(fid,record)
%
% Sub header read routine for ECAT 7 image files
%
% Roger Gunn, 260298
off                     = (record-1)*512;
status                  = fseek(fid, off,'bof');
data_type               = fread(fid,1,'uint16',0);
num_dimensions          = fread(fid,1,'uint16',0);
x_dimension             = fread(fid,1,'uint16',0);
y_dimension             = fread(fid,1,'uint16',0);
z_dimension             = fread(fid,1,'uint16',0);
x_offset                = fread(fid,1,'float32',0);
y_offset                = fread(fid,1,'float32',0);
z_offset                = fread(fid,1,'float32',0);
recon_zoom              = fread(fid,1,'float32',0);
scale_factor            = fread(fid,1,'float32',0);
image_min               = fread(fid,1,'int16',0);
image_max               = fread(fid,1,'int16',0);
x_pixel_size            = fread(fid,1,'float32',0);
y_pixel_size            = fread(fid,1,'float32',0);
z_pixel_size            = fread(fid,1,'float32',0);
frame_duration          = fread(fid,1,'uint32',0);
frame_start_time        = fread(fid,1,'uint32',0);
filter_code             = fread(fid,1,'uint16',0);
x_resolution            = fread(fid,1,'float32',0);
y_resolution            = fread(fid,1,'float32',0);
z_resolution            = fread(fid,1,'float32',0);
num_r_elements          = fread(fid,1,'float32',0);
num_angles              = fread(fid,1,'float32',0);
z_rotation_angle        = fread(fid,1,'float32',0);
decay_corr_fctr         = fread(fid,1,'float32',0);
corrections_applied     = fread(fid,1,'uint32',0);
gate_duration           = fread(fid,1,'uint32',0);
r_wave_offset           = fread(fid,1,'uint32',0);
num_accepted_beats      = fread(fid,1,'uint32',0);
filter_cutoff_frequency = fread(fid,1,'float32',0);
filter_resolution       = fread(fid,1,'float32',0);
filter_ramp_slope       = fread(fid,1,'float32',0);
filter_order            = fread(fid,1,'uint16',0);
filter_scatter_fraction = fread(fid,1,'float32',0);
filter_scatter_slope    = fread(fid,1,'float32',0);
annotation              = fread(fid,40,'char',0);
mt_1_1                  = fread(fid,1,'float32',0);
mt_1_2                  = fread(fid,1,'float32',0);
mt_1_3                  = fread(fid,1,'float32',0);
mt_2_1                  = fread(fid,1,'float32',0);
mt_2_2                  = fread(fid,1,'float32',0);
mt_2_3                  = fread(fid,1,'float32',0);
mt_3_1                  = fread(fid,1,'float32',0);
mt_3_2                  = fread(fid,1,'float32',0);
mt_3_3                  = fread(fid,1,'float32',0);
rfilter_cutoff          = fread(fid,1,'float32',0);
rfilter_resolution      = fread(fid,1,'float32',0);
rfilter_code            = fread(fid,1,'uint16',0);
rfilter_order           = fread(fid,1,'uint16',0);
zfilter_cutoff          = fread(fid,1,'float32',0);
zfilter_resolution      = fread(fid,1,'float32',0);
zfilter_code            = fread(fid,1,'uint16',0);
zfilter_order           = fread(fid,1,'uint16',0);
mt_4_1                  = fread(fid,1,'float32',0);
mt_4_2                  = fread(fid,1,'float32',0);
mt_4_3                  = fread(fid,1,'float32',0);
scatter_type            = fread(fid,1,'uint16',0);
recon_type              = fread(fid,1,'uint16',0);
recon_views             = fread(fid,1,'uint16',0);
fill                    = fread(fid,1,'uint16',0);
annotation              = deblank(char(annotation.*(annotation>0))');

SHEADER = struct('DATA_TYPE', data_type, ...
	'NUM_DIMENSIONS', num_dimensions, ...
	'X_DIMENSION', x_dimension, ...
	'Y_DIMENSION', y_dimension, ...
	'Z_DIMENSION', z_dimension, ...
	'X_OFFSET', x_offset, ...
	'Y_OFFSET', y_offset, ...
	'Z_OFFSET', z_offset, ...
	'RECON_ZOOM', recon_zoom, ...
	'SCALE_FACTOR', scale_factor, ...
	'IMAGE_MIN', image_min, ...
	'IMAGE_MAX', image_max, ...
	'X_PIXEL_SIZE', x_pixel_size, ...
	'Y_PIXEL_SIZE', y_pixel_size, ...
	'Z_PIXEL_SIZE', z_pixel_size, ...
	'FRAME_DURATION', frame_duration, ...
	'FRAME_START_TIME', frame_start_time, ...
	'FILTER_CODE', filter_code, ...
	'X_RESOLUTION', x_resolution, ...
	'Y_RESOLUTION', y_resolution, ...
	'Z_RESOLUTION', z_resolution, ...
	'NUM_R_ELEMENTS', num_r_elements, ...
	'NUM_ANGLES', num_angles, ...
	'Z_ROTATION_ANGLE', z_rotation_angle, ...
	'DECAY_CORR_FCTR', decay_corr_fctr, ...
	'CORRECTIONS_APPLIED', corrections_applied, ...
	'GATE_DURATION', gate_duration, ...
	'R_WAVE_OFFSET', r_wave_offset, ...
	'NUM_ACCEPTED_BEATS', num_accepted_beats, ...
	'FILTER_CUTOFF_FREQUENCY', filter_cutoff_frequency, ...
	'FILTER_RESOLUTION', filter_resolution, ...
	'FILTER_RAMP_SLOPE', filter_ramp_slope, ...
	'FILTER_ORDER', filter_order, ...
	'FILTER_SCATTER_CORRECTION', filter_scatter_fraction, ...
	'FILTER_SCATTER_SLOPE', filter_scatter_slope, ...
	'ANNOTATION', annotation, ...
	'MT_1_1', mt_1_1, ...
	'MT_1_2', mt_1_2, ...
	'MT_1_3', mt_1_3, ...
	'MT_2_1', mt_2_1, ...
	'MT_2_2', mt_2_2, ...
	'MT_2_3', mt_2_3, ...
	'MT_3_1', mt_3_1, ...
	'MT_3_2', mt_3_2, ...
	'MT_3_3', mt_3_3, ...
	'RFILTER_CUTOFF', rfilter_cutoff, ...
	'RFILTER_RESOLUTION', rfilter_resolution, ...
	'RFILTER_CODE', rfilter_code, ...
	'RFILTER_ORDER', rfilter_order, ...
	'ZFILTER_CUTOFF', zfilter_cutoff, ...
	'ZFILTER_RESOLUTION', zfilter_resolution, ...
	'ZFILTER_CODE', zfilter_code, ...
	'ZFILTER_ORDER', zfilter_order, ...
	'MT_4_1', mt_4_1, ...
	'MT_4_2', mt_4_2, ...
	'MT_4_3', mt_4_3, ...
	'SCATTER_TYPE', scatter_type, ...
	'RECON_TYPE', recon_type, ...
	'RECON_VIEWS', recon_views, ...
	'FILL', fill);
return;
%_______________________________________________________________________
function [MHEADER]=ECAT7_mheader(fid)
%
% Main header read routine for ECAT 7 image files
%
% Roger Gunn, 260298

status                  = fseek(fid, 0,'bof');
magic_number            = fread(fid,14,'char',0);
original_file_name      = fread(fid,32,'char',0);
sw_version              = fread(fid,1,'uint16',0);
system_type             = fread(fid,1,'uint16',0);
file_type               = fread(fid,1,'uint16',0);
serial_number           = fread(fid,10,'char',0);
scan_start_time         = fread(fid,1,'uint32',0);
isotope_name            = fread(fid,8,'char',0);
isotope_halflife        = fread(fid,1,'float32',0);
radiopharmaceutical     = fread(fid,32,'char',0);
gantry_tilt             = fread(fid,1,'float32',0);
gantry_rotation         = fread(fid,1,'float32',0);
bed_elevation           = fread(fid,1,'float32',0);
intrinsic_tilt           = fread(fid,1,'float32',0);
wobble_speed             = fread(fid,1,'uint16',0);
transm_source_type       = fread(fid,1,'uint16',0);
distance_scanned         = fread(fid,1,'float32',0);
transaxial_fov           = fread(fid,1,'float32',0);
angular_compression      = fread(fid,1,'uint16',0);
coin_samp_mode           = fread(fid,1,'uint16',0);
axial_samp_mode          = fread(fid,1,'uint16',0);
ecat_calibration_factor  = fread(fid,1,'float32',0);
calibration_units        = fread(fid,1,'uint16',0);
calibration_units_type   = fread(fid,1,'uint16',0);
compression_code         = fread(fid,1,'uint16',0);
study_type               = fread(fid,12,'char',0);
patient_id               = fread(fid,16,'char',0);
patient_name             = fread(fid,32,'char',0);
patient_sex              = fread(fid,1,'char',0);
patient_dexterity        = fread(fid,1,'char',0);
patient_age              = fread(fid,1,'float32',0);
patient_height           = fread(fid,1,'float32',0);
patient_weight           = fread(fid,1,'float32',0);
patient_birth_date       = fread(fid,1,'uint32',0);
physician_name           = fread(fid,32,'char',0);
operator_name            = fread(fid,32,'char',0);
study_description        = fread(fid,32,'char',0);
acquisition_type         = fread(fid,1,'uint16',0);
patient_orientation      = fread(fid,1,'uint16',0);
facility_name            = fread(fid,20,'char',0);
num_planes               = fread(fid,1,'uint16',0);
num_frames               = fread(fid,1,'uint16',0);
num_gates                = fread(fid,1,'uint16',0);
num_bed_pos              = fread(fid,1,'uint16',0);
init_bed_position        = fread(fid,1,'float32',0);
bed_position = zeros(15,1);
for bed=1:15,
	bed_position(bed)  = fread(fid,1,'float32',0);
end;
plane_separation         = fread(fid,1,'float32',0);
lwr_sctr_thres           = fread(fid,1,'uint16',0);
lwr_true_thres           = fread(fid,1,'uint16',0);
upr_true_thres           = fread(fid,1,'uint16',0);
user_process_code        = fread(fid,10,'char',0);
acquisition_mode         = fread(fid,1,'uint16',0);
bin_size                 = fread(fid,1,'float32',0);
branching_fraction       = fread(fid,1,'float32',0);
dose_start_time          = fread(fid,1,'uint32',0);
dosage                   = fread(fid,1,'float32',0);
well_counter_corr_factor = fread(fid,1,'float32',0);
data_units               = fread(fid,32,'char',0);
septa_state              = fread(fid,1,'uint16',0);
fill                     = fread(fid,1,'uint16',0);
	
magic_number             = deblank(char(magic_number.*(magic_number>32))');
original_file_name       = deblank(char(original_file_name.*(original_file_name>0))');
serial_number            = deblank(char(serial_number.*(serial_number>0))');
isotope_name             = deblank(char(isotope_name.*(isotope_name>0))');
radiopharmaceutical      = deblank(char(radiopharmaceutical.*(radiopharmaceutical>0))');
study_type               = deblank(char(study_type.*(study_type>0))');
patient_id               = deblank(char(patient_id.*(patient_id>0))');
patient_name             = deblank(char(patient_name.*(patient_name>0))');
patient_sex              = deblank(char(patient_sex.*(patient_sex>0))');
patient_dexterity        = deblank(char(patient_dexterity.*(patient_dexterity>0))');
physician_name           = deblank(char(physician_name.*(physician_name>0))');
operator_name            = deblank(char(operator_name.*(operator_name>0))');
study_description        = deblank(char(study_description.*(study_description>0))');
facility_name            = deblank(char(facility_name.*(facility_name>0))');
user_process_code        = deblank(char(user_process_code.*(user_process_code>0))');
data_units               = deblank(char(data_units.*(data_units>0))');

MHEADER = struct('MAGIC_NUMBER', magic_number, ...
	'ORIGINAL_FILE_NAME', original_file_name, ...
	'SW_VERSION', sw_version, ...
	'SYSTEM_TYPE', system_type, ...
	'FILE_TYPE', file_type, ...
	'SERIAL_NUMBER', serial_number, ...
	'SCAN_START_TIME', scan_start_time, ...
	'ISOTOPE_NAME', isotope_name, ...
	'ISOTOPE_HALFLIFE', isotope_halflife, ...
	'RADIOPHARMACEUTICAL', radiopharmaceutical, ... 
	'GANTRY_TILT', gantry_tilt, ...
	'GANTRY_ROTATION', gantry_rotation, ... 
	'BED_ELEVATION', bed_elevation, ...
	'INTRINSIC_TILT', intrinsic_tilt, ... 
	'WOBBLE_SPEED', wobble_speed, ...
	'TRANSM_SOURCE_TYPE', transm_source_type, ...
	'DISTANCE_SCANNED', distance_scanned, ...
	'TRANSAXIAL_FOV', transaxial_fov, ...
	'ANGULAR_COMPRESSION', angular_compression, ... 
	'COIN_SAMP_MODE', coin_samp_mode, ...
	'AXIAL_SAMP_MODE', axial_samp_mode, ...
	'ECAT_CALIBRATION_FACTOR', ecat_calibration_factor, ...
	'CALIBRATION_UNITS', calibration_units, ...
	'CALIBRATION_UNITS_TYPE', calibration_units_type, ...
	'COMPRESSION_CODE', compression_code, ...
	'STUDY_TYPE', study_type, ...
	'PATIENT_ID', patient_id, ...
	'PATIENT_NAME', patient_name, ...
	'PATIENT_SEX', patient_sex, ...
	'PATIENT_DEXTERITY', patient_dexterity, ...
	'PATIENT_AGE', patient_age, ...
	'PATIENT_HEIGHT', patient_height, ...
	'PATIENT_WEIGHT', patient_weight, ...
	'PATIENT_BIRTH_DATE', patient_birth_date, ...
	'PHYSICIAN_NAME', physician_name, ...
	'OPERATOR_NAME', operator_name, ...
	'STUDY_DESCRIPTION', study_description, ...
	'ACQUISITION_TYPE', acquisition_type, ...
	'PATIENT_ORIENTATION', patient_orientation, ...
	'FACILITY_NAME', facility_name, ...
	'NUM_PLANES', num_planes, ...
	'NUM_FRAMES', num_frames, ...
	'NUM_GATES', num_gates, ...
	'NUM_BED_POS', num_bed_pos, ...
	'INIT_BED_POSITION', init_bed_position, ...
	'BED_POSITION', bed_position, ...
	'PLANE_SEPARATION', plane_separation, ...
	'LWR_SCTR_THRES', lwr_sctr_thres, ...
	'LWR_TRUE_THRES', lwr_true_thres, ...
	'UPR_TRUE_THRES', upr_true_thres, ...
	'USER_PROCESS_CODE', user_process_code, ...
	'ACQUISITION_MODE', acquisition_mode, ...
	'BIN_SIZE', bin_size, ...
	'BRANCHING_FRACTION', branching_fraction, ...
	'DOSE_START_TIME', dose_start_time, ...
	'DOSAGE', dosage, ...
	'WELL_COUNTER_CORR_FACTOR', well_counter_corr_factor, ...
	'DATA_UNITS', data_units, ...
	'SEPTA_STATE', septa_state, ...
	'FILL', fill);
return;
%_______________________________________________________________________
