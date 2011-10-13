function [type, fname] = test_ft_senstype

% TEST test_ft_senstype
% TEST ft_senstype

% ft_senstype can work on different input data structures. Here, use the
% preprocessed data on /home/common/matlab/fieldtrip/test/raw/

cnt = 0;

cd('/home/common/matlab/fieldtrip/data/test/raw/');
cd('eeg');
d = dir('*.mat');
for k = 1:numel(d)
  cnt = cnt+1;
  load(d(k).name);
  fname{cnt} = d(k).name;
  type{cnt}  = ft_senstype(data);
end

cd('/home/common/matlab/fieldtrip/data/test/raw/');
cd('meg');
d = dir('*.mat');
for k = 1:numel(d)
  cnt = cnt+1;
  load(d(k).name);
  fname{cnt} = d(k).name;
  type{cnt}  = ft_senstype(data);
end

cd('/home/common/matlab/fieldtrip/data/test/raw/');
cd('lfp');
d = dir('*.mat');
for k = 1:numel(d)
  cnt = cnt+1;
  load(d(k).name);
  fname{cnt} = d(k).name;
  type{cnt}  = ft_senstype(data);
end

%if any(strmatch('unknown', type))
%  error('ft_senstype returns ''unknown'' sensor type for some test files on /home/common');
%end
