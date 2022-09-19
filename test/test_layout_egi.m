function test_layout_egi

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY

[ftver, ftpath] = ft_version;
dir_elec = fullfile(ftpath, 'template', 'electrode');
dir_lay  = fullfile(ftpath, 'template', 'layout');

cd(dir_elec);
d = dir('*.sfp');
for k = 1:numel(d)
  elec{1,k} = ft_read_sens(d(k).name);
end

cd(dir_lay);

% create a basic layout, with outline of head + circle
cfg         = [];
cfg.elec    = elec{1};
tmplay1     = ft_prepare_layout(cfg);

dum     = load('CTF151_helmet.mat');
tmplay2 = dum.lay;


for k = 1:numel(elec)
  layout = [];
  layout.mask = tmplay1.mask;
  layout.outline = [tmplay1.outline(1) tmplay2.outline([1 3:end])];
  
  cfg      = [];
  cfg.elec = elec{k};
  tmplay   = ft_prepare_layout(cfg);
  
  layout.pos    = tmplay.pos./0.85;
  layout.label  = tmplay.label;
  layout.width  = tmplay.width;
  layout.height = tmplay.height;
  
  if contains(d(k).name, '256') || contains(d(k).name, '257')
    for m = 2:numel(layout.outline)
      layout.outline{m}(:,2) = layout.outline{m}(:,2)-0.04;
    end
  elseif contains(d(k).name, '64') || contains(d(k).name, '65')
    layout.pos(:,2) = layout.pos(:,2) + 0.05;
    for m = 2:numel(layout.outline)
      layout.outline{m}(:,2) = layout.outline{m}(:,2)-0.04;
    end
  elseif contains(d(k).name, '32')
    
  end
  
  figure
  ft_plot_layout(layout);
  title(d(k).name, 'interpreter', 'none');
  lay = layout;
  % save(strrep(d(k).name,'sfp','mat'), 'lay');
end
