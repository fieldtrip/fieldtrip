function test_readcompresseddata

% MEM 4000mb
% WALLTIME 00:10:00

% TEST inflate_file ft_read_data ft_read_header test_readcompresseddata

% test these data sets
datasets = {
  dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/ctf151/Subject01.ds.zip')
  dccnpath('/home/common/matlab/fieldtrip/data/test/original/eeg/brainvision/MischaCompressed.zip')
};

for k = 1:numel(datasets)
  cfg = [];
  cfg.continuous = 'yes';
  cfg.dataset = datasets{k};
  data = ft_preprocessing(cfg);
end

% and this mri
mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/mri/dicom/dicomzipped.zip'));

end
