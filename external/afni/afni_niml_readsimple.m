function S=afni_niml_readsimple(fn)
% Reads an AFNI NIML file and returns a 'simple struct
%
% S=AFNI_NIML_READSIMPLE(FN) reads the AFNI NIML ASCII file FN, and returns
% a 'simple' struct S with the relevant data fields. (This is unlike
% the function AFNI_NIML_READ, which returns all data fields and resembles
% the original file content of FN better).
%
% If present, S will contain the following fields:
%  .data           PxN data for P nodes and N columns (values per node).
%  .node_indices   Px1 indices of P nodes that data refers to (base 0)
%  .history        String with history information
%  .stats          1xN cell with strings of the statistics of the data 
%                  (for example, Z or T score).
%  .labels         1xN cell with strings of the labels of the data columns
%  .dset_type      String with the data set type
%
% Please note that this function is *VERY EXPERIMENTAL*, and has only been
% tested with functional data NIML files.
% 
% NNO Jan 2010 <n.oosterhof@bangor.ac.uk>

D=afni_niml_read(fn);

if ~isstruct(D) || ~isfield(D,'dset_type')
    disp(D);
    error('Unrecognized data cell');
end

S=struct();
S.dset_type=D.dset_type;

nodecount=numel(D.nodes);
for k=1:nodecount
    Nk=D.nodes{k};
    if ~isfield(Nk,'name')
        disp(Dk);
        error('Missing name in node element %d',k);
    end
    
    switch(Nk.name)
        case 'INDEX_LIST'
            S.node_indices=Nk.data;
        case 'SPARSE_DATA'
            S.data=Nk.data;
        case 'AFNI_atr'
            if ~isfield(Nk,'atr_name')
                disp(Nk);
                error('Field %s in node element %d has missing atr_name\n', Nk.name, D);
            end
            
            switch Nk.atr_name
                case 'HISTORY_NOTE'
                    S.history=Nk.data;
                case 'COLMS_STATSYM'
                    S.stats=split_string(Nk.data,';');
                    if numel(S.stats{end})==0         %if final character is ';'
                        S.stats={S.stats{1:(end-1)}}; %then remove last element
                    end
                case 'COLMS_LABS'
                    S.labels=split_string(Nk.data,';');
                    if numel(S.labels{end})==0        % as above
                        S.labels={S.labels{1:(end-1)}};
                    end
            end
        otherwise
            error('Unsupported element %s in node element %d', Nk.name, D);
    end    
end
    




