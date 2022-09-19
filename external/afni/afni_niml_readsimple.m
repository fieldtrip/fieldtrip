function S=afni_niml_readsimple(fn,full)
% Reads an AFNI NIML file and returns a 'simple' struct
%
% S=AFNI_NIML_READSIMPLE(FN) reads the AFNI NIML ASCII file FN, and returns
% a 'simple' struct S with the relevant data fields. (This is unlike
% the function AFNI_NIML_READ, which returns all data fields and resembles
% the original file content of FN better).
%
% If the file is present, S will contain the following fields:
%  .data           PxN data for P nodes and N columns (values per node).
%  .node_indices   Px1 indices of P nodes that data refers to (base 0)
%  .history        String with history information
%  .stats          1xN cell with strings of the statistics of the data
%                  (for example, Z or T score).
%  .labels         1xN cell with strings of the labels of the data columns
%  .dset_type      String with the data set type
%
% S=AFNI_NIML_READSIMPLE(FN,true) returns a struct S but with a
% full data matrix (rather than sparse) for all nodes; node indices not
% defined in the file are set to zero. S does not contain a field
% .node_indices.
%
% Please note that this function is *VERY EXPERIMENTAL*, and has only been
% tested with functional data NIML files.
%
% NNO Jan 2010 <n.oosterhof@bangor.ac.uk>

D=afni_niml_read(fn);

if iscell(D)
    if numel(D)==1
        D=D{1};
    else
        error('Cell with multiple elements is not supported');
    end
end

if ~isstruct(D) || ~isfield(D,'dset_type')
    error('Unrecognized input');
end

S=struct();
S.dset_type=D.dset_type;

nodecount=numel(D.nodes);
for k=1:nodecount
    Nk=D.nodes{k};
    if ~isfield(Nk,'name')
        error('Missing name in node element %d',k);
    end

    switch(Nk.name)
        case 'INDEX_LIST'
            S.node_indices=Nk.data;

        case 'SPARSE_DATA'
            S.data=Nk.data;

        case 'AFNI_atr'
            if ~isfield(Nk,'atr_name')
                error('Field %s in node %d has missing atr_name\n', ...
                            Nk.name, D);
            end

            switch Nk.atr_name
                case 'HISTORY_NOTE'
                    S.history=get_string_cell_data(Nk.data);

                case 'COLMS_STATSYM'
                    S.stats=get_string_cell_data(Nk.data);

                case 'COLMS_LABS'
                    S.labels=get_string_cell_data(Nk.data);

            end

        otherwise
            S.misc.(Nk.name)=Nk.data;
    end
end

if nargin>1 && full && isfield(S,'node_indices')
    data=zeros(max(S.node_indices)+1,size(S.data,2));
    data(S.node_indices+1,:)=S.data;
    S.data=data;
    S=rmfield(S,'node_indices');
end


function labels=get_string_cell_data(data)
    while iscell(data) && numel(data)==1 && ~ischar(data)
        data=data{1};
    end

    labels=split_string(data,';');
    if numel(labels{end})==0        % as above
        labels=labels(1:(end-1));
    end



