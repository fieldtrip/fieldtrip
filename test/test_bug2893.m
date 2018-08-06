function test_bug2893

% WALLTIME 00:10:00
% MEM 2gb

% TEST ft_read_cifti

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2893/177746_MEG_Motort_tmegconne_[LM-TEMG-RH]_[CM-imcoh]_[FB-alpha].imcoh.Yeo11.pconnseries.nii');

cii = ft_read_cifti(filename);

% do some checks
assert(isequal(size(cii.brainordinate.pos), [8004 3]))
assert(numel(cii.brainordinate.brainstructure)==8004);
assert(numel(cii.brainordinate.parcellation)==8004);
assert(length(cii.brainordinate.brainstructurelabel)==max(cii.brainordinate.brainstructure)); % 2
assert(length(cii.brainordinate.parcellationlabel)==max(cii.brainordinate.parcellation)); % 107
