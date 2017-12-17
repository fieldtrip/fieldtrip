function test_bug2621

% MEM 2gb
% WALLTIME 00:30:00

% TEST ft_volumesegment


%read in the mri
mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/Subject01.mri'));

%segment the mri
cfg.output = 'brain';
seg  = ft_volumesegment(cfg, mri);


mriCM = ft_convert_units(mri,'cm');
segCM = ft_volumesegment(cfg, mriCM);

if (strcmp(seg.unit, 'mm') && strcmp(segCM.unit, 'mm'))
  fprintf('both segmentations had mm units\n')
else
  fprintf('one of the segmentations dooesn"t have mm units\n')
end

if (all(seg.brain(:) == segCM.brain(:)))
  fprintf('both segmentations have the same brain masks\n')
else
  fprintf('segmentations have differing brain masks\n')
end

if (all(seg.transform(:) == segCM.transform(:)))
  fprintf('both segmentations have the same transform\n')
else
  fprintf('segmentations have differing transforms\n')
end

% Comment: this is actually how it should behave, because the transform
% reflects the units. The difference is a factor of 10, which is as
% expected. Conclusion, this is not a bug.



