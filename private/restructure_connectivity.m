function data = restructure_connectivity(data)

% this function restructures connectivity data in the chancmb format to the
% chan_chan format. This required chancmb to contain values for all pairs
% of channels. Connectivity can be either symetric (undirected) or 
% assymetric (directed). 

if strcmp(data.dimord(1:7), 'pos_pos') || strcmp(data.dimord(1:9), 'chan_chan');
    % data is in chan_chan format. No implimentation yet for converting to
    % chancmb format
     error('Not yet implimented!');
    format=0;
elseif isfield(data,'labelcmb');
    % data is in chancmb format
    format=1;
else
    error('Unknown data format!');
end

n_chan=length(data.label);
n_chancmb=length(data.labelcmb);

% expected number of channels for directed connectivity
n_chan_cmb_dir= n_chan*n_chan-n_chan;
% expected number of channels for undirected connectivity
n_chan_cmb_und=n_chan_cmb_dir./2;

% detrimne if connectivity is directed or undirected
if n_chancmb==n_chan_cmb_dir
    n_chan_cmb=n_chan_cmb_dir;
    isdirected=1;
elseif n_chancmb==n_chan_cmb_und
    n_chan_cmb=n_chan_cmb_und;
    isdirected=0;
else
    error('Number of channel combnations is not as expected.');
end


% find fields with the right dimensions.
flds=fieldnames(data);
confld_idx=[];
for i=1:length(flds)
    if isnumeric(data.(flds{i}));
        siz=size(data.(flds{i}));
        if siz(1)==n_chan_cmb
            confld_idx=[confld_idx i];
            
        end
    end
end

% get the fields to transform
con_flds=flds(confld_idx);


for f=1:length(con_flds)
    data_temp=data.(con_flds{f});
    siz=size(data_temp);
    
    % transform the data!
    data_new=zeros([n_chan n_chan siz(2:end)]);
    k=0;
    if isdirected
        for i=1:n_chan
            for j=1:n_chan
                if i==j
                    continue
                end
                k=k+1;
                data_new(j,i,:)=data_temp(k,:);
            end
        end
        
    elseif ~isdirected
        for i=1:n_chan
            for j=i+1:n_chan
                k=k+1;
                data_new(i,j,:)=data_temp(k,:);
                data_new(j,i,:)=data_temp(k,:);
            end
        end
    end
    
    % put the transformed data back into the data structure
    data.(con_flds{f}) = data_new;
end

% update dimord
if strcmp(data.dimord(1:4),'chan');
    data.dimord=['chan_chan' data.dimord(5:end)];
elseif strcmp(data.dimord(1:3),'pos');
    data.dimord=['pos_pos' data.dimord(4:end)];
end

% remove labelcmb
data=rmfield(data,'labelcmb');



        

    



