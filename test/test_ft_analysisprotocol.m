function test_ft_analysisprotocol

% TEST test_ft_analysisprotocol
% TEST ft_analysisprotocol

% the style of this test script is also used in test_ft_datatype

dirlist = {
  '/home/common/matlab/fieldtrip/data/test/latest'
  '/home/common/matlab/fieldtrip/data/test/20111231'
  '/home/common/matlab/fieldtrip/data/test/20110630'
  '/home/common/matlab/fieldtrip/data/test/20101231'
  '/home/common/matlab/fieldtrip/data/test/20100630'
  '/home/common/matlab/fieldtrip/data/test/20091231'
  '/home/common/matlab/fieldtrip/data/test/20090630'
  '/home/common/matlab/fieldtrip/data/test/20081231'
  '/home/common/matlab/fieldtrip/data/test/20080630'
  '/home/common/matlab/fieldtrip/data/test/20071231'
  '/home/common/matlab/fieldtrip/data/test/20070630'
  '/home/common/matlab/fieldtrip/data/test/20061231'
  '/home/common/matlab/fieldtrip/data/test/20060630'
  '/home/common/matlab/fieldtrip/data/test/20051231'
  '/home/common/matlab/fieldtrip/data/test/20050630'
  '/home/common/matlab/fieldtrip/data/test/20040623'
  '/home/common/matlab/fieldtrip/data/test/20031128'
  };

for j=1:length(dirlist)
  filelist = hcp_filelist(dirlist{j});
  
  [~, ~, x] = cellfun(@fileparts, filelist, 'uniformoutput', false);
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
    ft_analysisprotocol(cfg, var);
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
