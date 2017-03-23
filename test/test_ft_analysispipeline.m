function test_ft_analysispipeline

% MEM 1500mb
% WALLTIME 02:53:11

% TEST ft_analysispipeline

% the style of this test script is also used in test_ft_datatype and test_bug2185

global ft_default
ft_default.trackusage = 'no'; % this calls ft_analysispipeline more than 6000 times, those should not all be tracked

dirlist = {
  dccnpath('/home/common/matlab/fieldtrip/data/test/latest')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20131231')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20130630')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20121231')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20120630')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20111231')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20110630')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20101231')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20100630')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20091231')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20090630')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20081231')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20080630')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20071231')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20070630')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20061231')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20060630')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20051231')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20050630')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20040623')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20031128')
  };

for j=1:length(dirlist)
  filelist = hcp_filelist(dirlist{j});
  
  [dummy, dummy, x] = cellfun(@fileparts, filelist, 'uniformoutput', false);
  sel = strcmp(x, '.mat');
  filelist = filelist(sel);
  clear p f x
  
  for i=1:length(filelist)
    
    % skip the large files
    d = dir(filelist{i});
    if d.bytes>50000000
      continue
    end
    
    try
      fprintf('processing data structure from %s\n', filelist{i});
      var = loadvar(filelist{i});
      disp(var)
    catch
      % some of the mat files are corrupt, this should not spoil the test
      disp(lasterr);
      continue
    end
    
    cfg = [];
    cfg.showcallinfo = 'no';
    ft_analysispipeline(cfg, var);
    set(gcf, 'Name', shortname(filelist{i}));
    drawnow
    close all
    
  end % for filelist
  
end % for dirlist
end % main function



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = shortname(str)
len = 50;
if length(str)>len
  begchar = length(str)-len+4;
  endchar = length(str);
  str = ['...' str(begchar:endchar)];
end
end % function shortname
