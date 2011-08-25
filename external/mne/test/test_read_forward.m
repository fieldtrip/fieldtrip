function test_read_forward()
% Test IO with FIF forward files

data_path = getenv('MNE_SAMPLE_DATASET_PATH');
if isempty(data_path)
    error('MNE_SAMPLE_DATASET_PATH environment variable no set.')
end

fname = [data_path filesep 'MEG' filesep 'sample' filesep 'sample_audvis-meg-eeg-oct-6-fwd.fif'];
fwd = mne_read_forward_solution(fname);