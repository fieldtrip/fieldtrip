function test_ex_compute_inverse()
% Test MNE dSPM inverse computation

data_path = getenv('MNE_SAMPLE_DATASET_PATH');
if isempty(data_path)
    error('MNE_SAMPLE_DATASET_PATH environment variable no set.')
end

fname_inv = [data_path filesep 'MEG' filesep 'sample' filesep 'sample_audvis-meg-oct-6-meg-inv.fif'];
fname_data = [data_path filesep 'MEG' filesep 'sample' filesep 'sample_audvis-ave.fif'];

setno = 1;
nave = -1;
dSPM = true;
lambda2 = 1/9;

[res] = mne_ex_compute_inverse(fname_data,setno,fname_inv,nave,lambda2,dSPM);

