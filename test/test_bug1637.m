function test_bug1637

% MEM 1500mb
% WALLTIME 00:10:00

% TEST megplanar_sincos channelconnectivity ft_prepare_neighbours ft_channelselection

% this function checks whether megplanar_sincos relies on a fixed channel
% order or whether this can be totally mixed up (it should be able to deal
% with that!)
%
% http://bugzilla.fcdonders.nl/show_bug.cgi?id=1637

% load neighbours
cfg = [];
cfg.method = 'template';
cfg.template = 'CTF275_neighb';
neighbours = ft_prepare_neighbours(cfg);

cd(dccnpath('/home/common/matlab/fieldtrip/data/test'))
load bug1637_hdr.mat
load bug1637_grad.mat

d = pwd;
ftPath = fileparts(mfilename('fullpath')); % get path, strip away 'ft_defaults'
ftPath = strrep(ftPath, '\', '\\');

cd(fullfile(ftPath, '..', 'private'));
errored = false;
try
  cfg = [];
  cfg.neighbours = neighbours;
  cfg.channel    = ft_channelselection('MEG', hdr.label);
  [neighbsel] = match_str({cfg.neighbours.label}, cfg.channel);
  cfg.neighbours = cfg.neighbours(neighbsel);
  cfg.neighbsel = channelconnectivity(cfg);
  
  montage_other = megplanar_sincos(cfg, hdr.grad);
  
  montage = [];
  montage{1} = megplanar_sincos(cfg, grad);
  
  [sel12, sel22] = match_str(montage{1}.labelold, montage_other.labelold);
  [sel11, sel21] = match_str(montage{1}.labelnew, montage_other.labelnew);
  if ~isequal(montage{1}.tra(sel11, sel12), montage_other.tra(sel21, sel22))
    warning('tra matrix differs - but that''s because the montage_other.tra consists of nans - checking whether nonzero elements are the same')
    idx1 = montage{1}.tra(sel11, sel12)~=0;
    idx2 = montage_other.tra(sel21, sel22)~=0;
    if ~isequal(idx1, idx2)
      errored = true;
      error('tra matrix qualitatively differ!');
    else
      disp('...passed!')
    end
  elseif ~(all(cellfun(@isequal, montage{1}.labelold(sel12), montage_other.labelold(sel22))))
    errored = true;
    error('labelold mismatch');
  elseif ~(all(cellfun(@isequal, montage{1}.labelnew(sel11), montage_other.labelnew(sel21))))
    errored = true;
    error('labelnew mismatch');
  end
  
  for i=1:10 % 10 is an arbitrary number here, just do it for some time
    labperm = randperm(numel(grad.label));
    % permute everything
    grad.chanori    = grad.chanori(labperm, :);
    grad.chanpos    = grad.chanpos(labperm, :);
    grad.chantype   = grad.chantype(labperm, :);
    grad.chanunit   = grad.chanunit(labperm, :);
    grad.coiloriri  = grad.coilori(labperm, :);
    grad.coilpos    = grad.coilpos(labperm, :);
    grad.label      = grad.label(labperm, :);
    grad.tra        = grad.tra(labperm, :);
    
    montage{2} = megplanar_sincos(cfg, grad);
    
    % since we now reordered everything according to cfg.channel, the
    % following should not be necessary - but hey!
    [sel12, sel22] = match_str(montage{1}.labelold, montage{2}.labelold);
    [sel11, sel21] = match_str(montage{1}.labelnew, montage{2}.labelnew);
    if ~(isequal(montage{1}.tra(sel11, sel12), montage{2}.tra(sel21, sel22)))
      errored = true;
      error('tra matrix differs')
    elseif ~(all(cellfun(@isequal, montage{1}.labelold(sel12), montage{2}.labelold(sel22))))
      errored = true;
      error('labelold mismatch');
    elseif ~(all(cellfun(@isequal, montage{1}.labelnew(sel11), montage{2}.labelnew(sel21))))
      errored = true;
      error('labelnew mismatch');
    end
  end
  
  
  for i=1:10 % 10 is an arbitrary number here, just do it for some time
    labperm = randperm(numel(grad.label));
    % permute everything
    grad.chanori    = hdr.grad.chanori(labperm, :);
    grad.chanpos    = hdr.grad.chanpos(labperm, :);
    grad.chantype   = hdr.grad.chantype(labperm, :);
    grad.chanunit   = hdr.grad.chanunit(labperm, :);
    grad.coiloriri  = hdr.grad.coilori(labperm, :);
    grad.coilpos    = hdr.grad.coilpos(labperm, :);
    grad.label      = hdr.grad.label(labperm, :);
    grad.tra        = hdr.grad.tra(labperm, :);
    
    montage{2} = megplanar_sincos(cfg, hdr.grad);
    
    % since we now reordered everything according to cfg.channel, the
    % following should not be necessary - but hey!
    [sel12, sel22] = match_str(montage{1}.labelold, montage{2}.labelold);
    [sel11, sel21] = match_str(montage{1}.labelnew, montage{2}.labelnew);
    idx1 = montage{1}.tra(sel11, sel12)~=0;
    idx2 = montage_other.tra(sel21, sel22)~=0;
    if ~isequal(idx1, idx2)
      errored = true;
      error('tra matrix qualitatively differ!');
    elseif ~(all(cellfun(@isequal, montage{1}.labelold(sel12), montage_other.labelold(sel22))))
      errored = true;
      error('labelold mismatch');
    elseif ~(all(cellfun(@isequal, montage{1}.labelnew(sel11), montage_other.labelnew(sel21))))
      errored = true;
      error('labelnew mismatch');
    end
  end
  
end
cd(d)

if errored
  error('channelorder matters whereas it should not!');
end

