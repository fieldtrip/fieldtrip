function inspect_pull1413

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_rejectvisual rejectvisual_summary

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/pull1413'));

%%

load datameg

%%

cfg = [];
% cfg.channel = 'MEG';
cfg.method = 'summary';
cfg.layout = 'ctf151_helmet.mat';
cfg.neighbours = 'ctf151_neighb.mat';
data_clean = ft_rejectvisual(cfg, data);

%%

load datanirs

%%

clear allrx alltx wavelength

for i=1:numel(data.label)
  tok = strsplit(data.label{i}, {'-', ' '});
  if numel(tok)==3
    allrx{i} = tok{1};
    alltx{i} = tok{2};
    wavelength{i} = tok{3};
  else
    allrx{i} = nan;
    alltx{i} = nan;
    wavelength{i} = nan;
  end
end


%%

clear neighbours

for i=1:numel(data.label)
  % see ft_prepare_neighbours
  
  % search for the same Rx and Tx optode
  thisrx = allrx{i};
  thistx = alltx{i};
  labrx = data.label(strcmp(allrx, thisrx));
  labtx = data.label(strcmp(alltx, thistx));
  labsame = data.label(strcmp(allrx, thisrx) & strcmp(alltx, thistx));
  
  
  % neighbours(i).neighblabel = [labrx; labtx];
  neighbours(i).neighblabel = labsame;
  % neighbours(i).neighblabel = setdiff([labrx; labtx], labsame);
  
  neighbours(i).label = data.label{i};
  % do not include the channel itself
  neighbours(i).neighblabel = setdiff(neighbours(i).neighblabel, neighbours(i).label);
  
end

%%

cfg = [];
cfg.method = 'summary';
cfg.neighbours = neighbours;
data_clean = ft_rejectvisual(cfg, data);
