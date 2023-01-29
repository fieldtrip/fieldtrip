function test_ft_write_mri

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_read_mri ft_write_mri

mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.mri'));
mri = ft_convert_units(mri, 'mm');

%%

dataformat = {
  'analyze'
  'analyze_hdr'
  'analyze_img'
  'analyze_old'
  'freesurfer_mgz'
  'mgh'
  'mgz'
  'nifti'
  'nifti2'
  'nifti_gz'
  'nifti_img'
  'nifti_spm'
  'vmp'
  'vmr'
  };

%%

p = tempname;

for i=1:numel(dataformat)
  filename = fullfile(p, ['Subject01_' dataformat{i}]);
  % pass it as separate fields
  ft_write_mri(filename, mri.anatomy, 'dataformat', dataformat{i}, 'transform', mri.transform, 'unit', mri.unit, 'coordsys', mri.coordsys);
end

% remove the temporary directory with all files inside it
rmdir(p, 's');

%%

p = tempname;

for i=1:numel(dataformat)
  filename = fullfile(p, ['Subject01_' dataformat{i}]);
  % pass it as a structure
  ft_write_mri(filename, mri, 'dataformat', dataformat{i});
end

% remove the temporary directory with all files inside it
rmdir(p, 's');

%%

p = tempname;

% scale between -1 and 1, we will skip the unsigned formats
dat = double(mri.anatomy);
dat = dat./max(abs(dat(:)));

mri.anatomy = double(dat);                         filename = fullfile(p, 'Subject01_double'); ft_write_mri(filename, mri, 'dataformat', 'analyze_old');
mri.anatomy = single(dat);                         filename = fullfile(p, 'Subject01_single'); ft_write_mri(filename, mri, 'dataformat', 'analyze_old');
mri.anatomy = int32(dat.*double(intmax('int32'))); filename = fullfile(p, 'Subject01_int32');  ft_write_mri(filename, mri, 'dataformat', 'analyze_old');
mri.anatomy = int16(dat.*double(intmax('int16'))); filename = fullfile(p, 'Subject01_int16');  ft_write_mri(filename, mri, 'dataformat', 'analyze_old');

mri.anatomy = double(dat);                         filename = fullfile(p, 'Subject01_double'); ft_write_mri(filename, mri, 'dataformat', 'nifti');
mri.anatomy = single(dat);                         filename = fullfile(p, 'Subject01_single'); ft_write_mri(filename, mri, 'dataformat', 'nifti');
mri.anatomy = int32(dat.*double(intmax('int32'))); filename = fullfile(p, 'Subject01_int32');  ft_write_mri(filename, mri, 'dataformat', 'nifti');
mri.anatomy = int16(dat.*double(intmax('int16'))); filename = fullfile(p, 'Subject01_int16');  ft_write_mri(filename, mri, 'dataformat', 'nifti');
mri.anatomy = uint8(dat.*double(intmax('int8')));  filename = fullfile(p, 'Subject01_int8');   ft_write_mri(filename, mri, 'dataformat', 'nifti');

% remove the temporary directory with all files inside it
rmdir(p, 's');
