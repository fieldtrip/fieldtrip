function test_issue1618

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_volumewrite ft_write_mri

% this function was related to github issue 1618

%%
mri      = [];
mri.dim  = [2 2 2];
mri.transform = eye(4);
mri.unit = 'mm';
mri.coordsys = 'ras';

mri.dat1 = randn(2,2,2);
mri.dat2 = reshape(1:8,[2 2 2]);
mri.dat3 = rand(2,2,2)>0.5;

datatype = {'logical', 'uint8', 'int8', 'int16', 'int32', 'single', 'double'};
filetype = {'analyze_old', 'nifti', 'nifti_img', 'nifti_spm', 'nifti_gz', 'analyze_spm', 'mgz', 'mgh', 'vmp', 'vmr'};
fname    = fullfile(tempdir, 'issue1618');

ok_table = true(numel(datatype), numel(filetype), 3);
for m = 1:numel(filetype)
  for k = 1:numel(datatype)
    try
      cfg = [];
      cfg.datatype  = datatype{k};
      cfg.filetype  = filetype{m};
      cfg.parameter = 'dat1';
      cfg.filename  = [fname, '_', datatype{k}, '_', filetype{m},'_', cfg.parameter];
      %if k==2, keyboard; end
      ft_volumewrite(cfg, mri);
    catch
      ok_table(k,m,1) = false;
    end
    fnames{k,m,1} = cfg.filename;
    
    try
      cfg = [];
      cfg.datatype  = datatype{k};
      cfg.filetype  = filetype{m};
      cfg.parameter = 'dat2';
      cfg.scaling   = 'no';
      cfg.filename = [fname, '_', datatype{k}, '_', filetype{m},'_', cfg.parameter];
      ft_volumewrite(cfg, mri);
    catch
      ok_table(k,m,2) = false;
    end
    fnames{k,m,2} = cfg.filename;
    
    try
      cfg = [];
      cfg.datatype  = datatype{k};
      cfg.filetype  = filetype{m};
      cfg.parameter = 'dat3';
      cfg.scaling   = 'no';
      cfg.filename = [fname, '_', datatype{k}, '_', filetype{m},'_', cfg.parameter];
      ft_volumewrite(cfg, mri);
    catch
      ok_table(k,m,3) = false;
    end
    fnames{k,m,3} = cfg.filename;
  end
end

data = cell(size(fnames));
for k = 1:size(fnames,1)
  for kk = 1:size(fnames,2)
    for kkk = 1:size(fnames,3)
      if ok_table(k,kk,kkk) && ~contains(fnames{k,kk,kkk},'vmr') && ~contains(fnames{k,kk,kkk},'vmp')
        % get the filename + extension but only if file has been written, and
        % if it's not a brainvoyager type, because no reading functionality
        % availablt
        d = dir([fnames{k,kk,kkk} '*']);
        if numel(d)>1
          % could be a hdr/img combo, take the first one
          d = d(1);
        end
        fname = fullfile(d.folder,d.name);
        data{k,kk,kkk} = ft_read_mri(fname);
      end
    end
  end
end

% FIXME this ideally should contain some check as to whether the values
% that have been written to disk make sense and/or are a scaled version of
% the original input data