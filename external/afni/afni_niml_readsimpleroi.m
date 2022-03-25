function S=afni_niml_readsimpleroi(fn)
% Reads a AFNI NIML ROI dataset (extension ".niml.roi")
%
% S=AFNI_NIML_READSIMPLEROI(FN) reads an AFNI NIML .niml.roi dataset and
% returns a struct S with the ROIs.
% If FN specifies N rois, then S is a Nx1 cell with structs with the
% following fields:
%   .Label        string of ROI label
%   .iLabel       integer value of label
%   .edge         1xP cell with node indices of ROI edges
%   .region       1xQ cell with node indices contained in ROI (typically Q==1)
%
% Node indices are returned in base 1 (first node has index 1)
%
% Please note that this function is *VERY EXPERIMENTAL*, and has only been
% tested with some ROI NIML files generated with AFNI SUMA.
%
% SEE ALSO: AFNI_NIML_READSIMPLE, AFNI_NIML_READ
%
% NNO Nov 2010

D=afni_niml_read(fn);
isastruct=isstruct(D);
if isastruct
    D={D};
end

% typenames of ROI - so far I only understand two of them
tpnames={'region','unknown2','unknown3','edge'};
ntpnames=numel(tpnames);

n=numel(D);
S=cell(n,1); % multiple ROIs

for k=1:n
    Dk=D{k};

    Sk=struct();
    Sk.Label=Dk.Label;
    Sk.iLabel=str2double(Dk.iLabel);

    idxs=sscanf(Dk.data,'%f');
    nidxs=numel(idxs);
    pos=0; % last position of previous list of indices
    while pos<nidxs
        tp=idxs(pos+2); % type: 4 is line segment, 1 is list of nodes in ROI
        nj=idxs(pos+3); % number of elements to follow
        nodeidxs=idxs(pos+3+(1:nj))+1; % convert to base 1
        if tp>ntpnames
            error('Unrecognized ROI type %d at pos %d', tp, pos+2);
        end
        tpname=tpnames{tp};
        if ~isfield(Sk,tpname)
            Sk.(tpname)=cell(0); % empty cell
        end
        Sk.(tpname){end+1}=nodeidxs;
        pos=pos+nj+3;
    end
    S{k}=Sk;
end


