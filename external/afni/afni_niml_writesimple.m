function D=afni_niml_writesimple(S,fn)
% writes surface data in a 'simple' struct in NIML format to a file.
%
% D=AFNI_NIML_WRITESIMPLE(fn,S) writes surface data S to the file FN.
% S should be a struct with a field S.data, and optionally S.stats,
% S.labels, S.history, S.types, and S.node_indices. 
% This function returns the more complicated NIML-struct D that is 
% written by afni_niml_write
%
% S.data            a PxN struct for P nodes and N datapoints per node.
% S.node_indices    (optional) a Nx1 vector with indices in base 0. If 
%                   omitted, then it is assumed that all nodes have data. 
% S.stats       }   Either a string (with ';'-separated elements for 
% S.labels      }   S.stats and .labels: , or a cell with strings. 
% S.history     }   In the latter case, if the number of elements in the 
%                   cell is less than P, then elements are taken in cycles.
%
% Example: 
%   S=struct()
%   S.data=randn(100002,4);
%   S.labels={'A-mean','A-Tscore','B-mean','B-Tscore'};
%   S.stats={'none','Ttest(15)'}; % each stat descriptor will be used twice
%   afni_niml_writesimple('test.niml.dset',S)         
% 
% NNO Feb 2010 <n.oosterhof@bangor.ac.uk>

if isnumeric(S)
    T=struct();
    T.data=S;
    S=T;
    clear T;
end

if ischar(S) && isstruct(fn)
    % swap struct and filename
    T=fn;
    fn=S;
    S=T;
    clear T;
end

if ~(isstruct(S) && isfield(S,'data'))
    error('Illegal input: expected struct S with field S.data\n');
end

% main header
D=struct();
if isfield(S,'dset_type')
    D.dset_type=S.dset_type;
else
    D.dset_type='Node_Bucket'; % this is the default
end

% make new self_idcode, to keep SUMA happy
rand('twister',sum(100*clock));
D.self_idcode=['XYZ_' char(rand(1,24)*26+65)];

[path,f,ext]=fileparts(fn);
D.filename=[f ext];
D.label=D.filename;
D.name='AFNI_dataset';
D.ni_form='ni_group';

nodes=cell(7,1);

%data
nd=struct();
nd.data_type='Node_Bucket_data';
nd.name='SPARSE_DATA';
nd.data=S.data;
[nverts,ncols]=size(S.data);
nodes{1}=nd;

% node indices
nd=struct();
nd.data_type='Node_Bucket_node_indices';
if isfield(S,'node_indices')
    idxs=S.node_indices;
else
    idxs=0:(nverts-1);
end

if issorted(idxs)
    nd.sorted_node_def='Yes';
else
    nd.sorted_node_def='No';
end

if size(idxs,1)==1
    idxs=idxs';
end

if numel(idxs) ~= nverts
    error('The number of node indices (%d) does not match the number of rows (%d) in the data.',nverts,numel(idxs));
end

nd.data=idxs;
nd.COLMS_RANGE=data2range(idxs,0);
nd.COLMS_LABS='Node Indices';
nd.COLMS_TYPE='Node_Index';
nd.name='INDEX_LIST';
nodes{2}=nd;

% colum range
nd=struct();
nd.atr_name='COLMS_RANGE';
nd.name='AFNI_atr';
nd.data=data2range(S.data);
nodes{3}=nd;

% default labels for columns
defaultlabels=cell(1,ncols);
for k=1:ncols
    defaultlabels{k}=sprintf('col_%d',k-1);
end

% default value for history; make call stack
st=dbstack();
stc=cell(numel(st));
for j=1:numel(st)
    stc{j}=sprintf('%s (%d)',st(j).name,st(j).line);
end
defaulthistory=sprintf('Written %s using: %s',datestr(clock),unsplit_string(stc,' <- '));
if isfield(S,'history')
    S.history=[S.history ';' defaulthistory];
end

nodes{4}=make_str_element('COLMS_LABS',S,'labels',defaultlabels,ncols);
nodes{5}=make_str_element('COLMS_TYPE',S,'types',{'Generic_Float'},ncols);
nodes{6}=make_str_element('COLMS_STATSYM',S,'stats',{'none'},ncols);
nodes{7}=make_str_element('HISTORY_NOTE',S,'history',defaulthistory,ncols);

D.nodes=nodes;

afni_niml_write(D,fn);

function elem=make_str_element(Nodefieldname,S,Sfieldname,default,ncols)
elem=struct();
elem.atr_name=Nodefieldname;
elem.name='AFNI_atr';
if isfield(S,Sfieldname)
    vals=S.(Sfieldname);
else
    vals=default;
end

if iscell(vals);
    valcount=numel(vals);
    v=cell(1,ncols);
    for k=1:ncols
        v{k}=vals{mod(k-1,valcount)+1}; %cyclically take elements
    end
    vals=unsplit_string(v,';');
elseif ~ischar(vals)
    error('Unrecognized data type for field %s (%s)', Nodefieldname, Sfieldname);
end

elem.data=vals;


function r=data2range(data,precision)
% returns a string that defines the range of the data, and nodes where this
% range is found.

[rows,cols]=size(data);
rs=cell(cols,1);

for k=1:cols
    datak=data(:,k);
    [minv,mini]=min(datak);
    [maxv,maxi]=max(datak);
    if nargin<2
        precision=get_required_precision(data(:,k));
    end
    pat=sprintf('%%.%df %%.%df %%d %%d', precision, precision);
    rs{k}=sprintf(pat, minv, maxv, mini-1, maxi-1);
end

r=unsplit_string(rs,';');