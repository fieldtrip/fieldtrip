function test_historical

% this function creates the datafiles for the historical fieldtrip versions

% it should not run unattended
% return

version = {
   '20111231'
   '20110630'
   '20101231'
   '20100630' % from this one on it works quite well
   '20091231' % this one and previous don't have inputfile/outputfile
   '20090630' % this one and previous don't have the ft-prefix
   '20081231'
   '20080630'
   '20071231'
   '20070630'
   '20061231'
   '20060630'
   '20051231'
   '20050630'
   '20040623'
   '20031128'
};

writeflag = true;
datainfo  = test_datasets;
datainfo = datainfo(10); % this is the ctf275 dataset

for i=1:length(version)
  restoredefaultpath
  clear functions

  % avoid using ft_defaults or fieldtripdefs
  addpath(genpath(sprintf('/home/common/matlab/fieldtrip-%s', version{i})));

  for j=1:length(datainfo)
    try
      test_ft_preprocessing(datainfo(j), writeflag, version{i});
    catch me
      disp(me)
    end % try
  end % for datainfo

  for j=1:length(datainfo)
    try
      test_ft_timelockanalysis(datainfo(j), writeflag, version{i});
    catch me
      disp(me)
    end % try
  end % for datainfo

  for j=1:length(datainfo)
    try
      test_ft_freqanalysis(datainfo(j), writeflag, version{i});
    catch me
      disp(me)
    end % try
  end % for datainfo

  for j=1:length(datainfo)
    try
      test_ft_sourceanalysis(datainfo(j), writeflag, version{i});
    catch me
      disp(me)
    end % try
  end % for datainfo

end % for version

