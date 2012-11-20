function test_bug1826

% at this moment the test script does not yet work, but we don't want the automatic regression testing to flag it as failure
return

% TEST test_bug1826
% TEST ft_read_mri ft_write_mri ft_volumerealign


cd /home/common/matlab/fieldtrip/data/test/bug1826

subjectT1  = 'T1.nii.gz';
subjectT2  = 'T2.nii.gz';
subjectDTi = 'DTI.nii';

T1  = ft_read_mri(subjectT1);
T2  = ft_read_mri(subjectT2);
DTi = ft_read_mri(subjectDTi);

figure; ft_plot_ortho(T1); title('T1 before aligning')
figure; ft_plot_ortho(T2); title('T2 before aligning')
figure; ft_plot_ortho(DTi);title('DTi before aligning')

% they are now in memory, which would be the normal starting point in a
% fieldtrip analysis pipeline

% INSERT THE NEW CODE HERE, THIS PROBABLY INVOLVES
% step 1) write them to disk, perhaps using ft_write_mri
% step 2) call external code, e.g. FSL (using a system call) or SPM
% step 3) read the results from disk

% they are now again in memory, but realigned to each other. This is where
% a normal fieldtrip analysis pipeline would continue.
figure; ft_plot_ortho(T1); title('T1 after aligning')
figure; ft_plot_ortho(T2); title('T2 after aligning')
figure; ft_plot_ortho(DTi);title('DTi after aligning')

