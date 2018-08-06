function test_bug2513

% WALLTIME 00:40:00
% MEM 6gb

datapath = dccnpath('/home/common/matlab/fieldtrip/data/test');

% historical preprocessed data
datahist = dir([datapath filesep '2*']);
filelist = {};

%% loop recursively through file directories
for i=1:size(datahist,1);
  if ispc
    dirs = regexp(genpath([datapath filesep datahist(i,1).name filesep]),['[^;]*.'],'match');
  else
    dirs = regexp(genpath([datapath filesep datahist(i,1).name filesep]),['[^:]*.'],'match');
  end
  
  for j=1:size(dirs,2);
    subfiles = dir([dirs{1,j}(1:end-1) filesep '*.mat']);
    if ~isempty(subfiles);
      for h=1:size(subfiles,1);
        subfiles(h,1).name = [dirs{1,j}(1:end-1) filesep subfiles(h,1).name];
      end
      filelist=cat(1,filelist,subfiles.name);
    end
  end
end

%% collect all the fields
collection.timelock = {};
collection.freq = {};
collection.raw = {};
collection.comp = {};
collection.grad = {};
collection.elec = {};
collection.source = {};
collection.volume = {};
collection.unknown = {};

% collect the files that ft_datatype cannot parse
unknown = {};
failed = [];

%%
datapath = dccnpath('/home/common/matlab/fieldtrip/data/test');

for i=1:length(filelist)
  try
    tmp = load(filelist{i,1});
    fn = fieldnames(tmp);
    if length(fn)>1
      error('multiple structures in the mat file');
    else
      tmp = tmp.(fn{1});
    end
    
    dtype = ft_datatype(tmp);
    collection.(dtype) = cat(2, collection.(dtype), fieldnames(tmp)');
    collection.(dtype) = unique(collection.(dtype));
    if strcmp(dtype, 'unknown')
      unknown{end+1} = filelist{i,1};
    end
    
  catch
    fprintf('failed for %s\n', filelist{i});
    failed(end+1)=i;
  end % try
end % for


%% loop the failed files
 failed2 = [];
for i=1:length(failed)
  try
    tmp = load(fullfile(datapath,filelist{failed(i),1}));
    fn = fieldnames(tmp);
    if length(fn)>1
      error('multiple structures in the mat file');
    else
      tmp = tmp.(fn{1});
    end
    
    dtype = ft_datatype(tmp);
    collection.(dtype) = cat(2, collection.(dtype), fieldnames(tmp)');
    collection.(dtype) = unique(collection.(dtype));
    if strcmp(dtype, 'unknown')
      unknown{end+1} = filelist{failed(i),1};
    end
    
  catch
    fprintf('failed for %s\n', filelist{failed(i)});
    failed2(end+1)=failed(i);
  end % try
end % for
