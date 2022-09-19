function test_bug3379

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_redefinetrial ft_fetch_data

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3379.mat'));


assert(isequal(tmpcfg.trl(:,1:2), testdata.sampleinfo(setdiff(1:93, [5 6]),:)));

% Since the tmpcfg.trl is the same as testdata.sampleinfo, the call to
% ft_redefinetrial is expected to lead to the same data, which it doesn't
% note: the data contains overlapping trials.
% Yet, the overlapping segments of the data are not identical, which should
% have led to an error in ft_fetch_data, which it didn't. This has now been
% fixed in ft_fetch_data, so ft_redefinetrial correctly returns an error
% now, since the discrepancy is detected in ft_fetch_data

try
  testdata2 = ft_redefinetrial(tmpcfg, testdata);
  % at this point, there should have been an informative error in FieldTrip
  failed = false;
catch
  failed = true;
end
assert(failed,'some of the requested samples occur twice in the data and have conflicting values');
