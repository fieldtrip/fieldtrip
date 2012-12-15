function [gridout, cfg] = sb_check_sources(cfg, vol, gridin)

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

cfg.sourcelabel = ft_getopt(cfg,'sourcelabel');
cfg.corr = ft_getopt(cfg, 'corr');

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

nnod = size(vol.pnt,1);

nsource = size(gridin,1);

% compute node element assignment
node_ele_assignment   = revert_assignment(elem,0.5);

if isfield(cfg,'sourcelabel') && isfield(vol,'tissue') && isfield(vol,'tissuelabel')
         if length(vol.tissue) == size(elem,1)
             numsource = length(cfg.sourcelabel);
             source2num = [];
             for i=1:numsource
                 source2num(i) = str2num(cfg.sourcelabel{i});
             end                 
             numlabels = length(vol.tissuelabel);
             sourcetissue = [];
             for i=1:numlabels
                 if ismember(str2num(vol.tissuelabel{i}),source2num)
                     sourcetissue = [sourcetissue; i];
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
valid_nodes_pnt = vol.pnt(valid_nodes,:);

gridout = zeros(size(gridin));
griddel = zeros(size(gridin,1),1);

maxmov = 0;
meanmov = 0;
nmov = 0;

time = clock;
timing = zeros(nsource,1);
for i=1:nsource
    [~,near_node_ind] = min(sum(bsxfun(@minus,vol.pnt,gridin(i,:)).^2,2));
    if(pos_nodes_log_ind(near_node_ind))
        gridout(i,:) = gridin(i,:);
    elseif strcmp(cfg.corr,'delete')
        griddel(i) = i;
    else
        nmov = nmov + 1;
    end
end

griddel = griddel(griddel ~= 0);
if strcmp(cfg.corr,'delete')
    gridout(griddel,:)=[];
    nmov = length(griddel);
end

if(meanmov ~= 0)
    if strcmp(cfg.corr,'delete')
        warning('%d source positions have been deleted.', nmov);
    elseif nmov > 0
        warning('%d source positions are not inside the designated source compartments! This may lead to high numerical errors!', nmov);
    end
else
    warning('All source positions ok.');
end
end
