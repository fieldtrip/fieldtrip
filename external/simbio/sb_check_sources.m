function [inside, cfg] = sb_check_sources(cfg, vol, gridin)

% SB_CHECK_SOURCES
%
% Input: cfg, FE mesh, grid positions
% Output: valid grid positions
%
% cfg.sourcelabel = cell-array containing the labels of the valid source
% compartments
%
% cfg.corr = 'delete' to delete invalid source positions, otherwise only a
% warning will be given
%
% Copyright (C) 2012, Felix Lucka, Johannes Vorwerk
%
% $Id$

cfg.sourcelabel = ft_getopt(cfg.grid,'sourcelabel');

% find connectivity information and write them to elem
if isfield(vol,'tet')
    if size(vol.tet,1) == 4
        mele = size(vol.tet,1);
        elem = vol.tet';
    elseif size(vol.tet,2) == 4
        mele = size(vol.tet,2);
        elem = vol.tet;
    else
        error('vol.tet has wrong dimensions!')
    end
elseif isfield(vol,'hex')
    if size(vol.hex,1) == 8
        mele = size(vol.hex,1);
        elem = vol.hex';
    elseif size(vol.hex,2) == 8
        mele = size(vol.hex,2);
        elem = vol.hex;
    else
        error('vol.hex has wrong dimensions!')
    end
else
    error('Could not find connectivity information!')
end

nnod = size(vol.pos,1);

nsource = size(gridin,1);

% compute node element assignment
node_ele_assignment   = revert_assignment(elem,0.5);

if isfield(cfg,'sourcelabel') && isfield(vol,'tissue') && isfield(vol,'tissuelabel')
         if length(vol.tissue) == size(elem,1)
             numsource = length(cfg.sourcelabel);         
             numlabels = length(vol.tissuelabel);
             sourcetissue = [];
             for i=1:numlabels
                 for j=1:numsource
                     if strcmpi(vol.tissuelabel{i},cfg.sourcelabel{j})
                         sourcetissue = [sourcetissue; i];
                     end
                 end
             end
             source_log_ind = ismember(vol.tissue,sourcetissue);
        else
            error('Dimensions of vol.tet/vol.hex and vol.tissue do not fit!');
         end
else
    error('Not all needed fields present!');
end

com_nodes = unique(elem(source_log_ind,:));
pos_nodes_log_ind = false(nnod,1);

for i=1:length(com_nodes)
    ind_ele =  nonzeros(node_ele_assignment(:,com_nodes(i)));
    if(all(ismember(vol.tissue(ind_ele),sourcetissue)))
        pos_nodes_log_ind(com_nodes(i)) = true;
    end
end

valid_nodes = find(pos_nodes_log_ind);

gridout = zeros(size(gridin));
%inside = zeros(size(gridin,1),1);

tic;

near_node_ind = zeros(nsource,1);

if (ft_hastoolbox('STATISTICS'))
    near_node_ind = knnsearch(vol.pos,gridin);
else
    warning('Statistics toolbox not available, using fallback routine. This might take a while.')
    parfor i=1:nsource
        [dummy,near_node_ind(i)] = min(sum(bsxfun(@minus,vol.pos,gridin(i,:)).^2,2));
    end
end

clear dummy;

toc

inside = 1:nsource;

inside(find(~pos_nodes_log_ind(near_node_ind))) = [];

% 
% griddel = griddel(griddel ~= 0);
% if strcmp(cfg.corr,'delete')
%     gridout(griddel,:)=[];
%     nmov = length(griddel);
% end
% 
% if(meanmov ~= 0)
%     if strcmp(cfg.corr,'delete')
%         warning('%d source positions have been deleted.', nmov);
%     elseif nmov > 0
%         warning('%d source positions are not inside the designated source compartments! This may lead to high numerical errors!', nmov);
%     end
% else
%     warning('All source positions ok.');
% end
end
