function o = nifti2struc
% Create a data structure describing NIFTI-2 headers
%__________________________________________________________________________
% Copyright (C) 2005-2012 Wellcome Trust Centre for Neuroimaging

%
% $Id: nifti2struc.m 6314 2015-01-23 17:00:51Z guillaume $


persistent org;
if ~isempty(org)
    o = org;
    return;
end
t = struct(...
    'conv',{@char  , @uint8 , @int16 , @int32 , @int64 , @single , @double },...
    'prec',{'uint8', 'uint8', 'int16', 'int32', 'int64', 'single', 'double'},...
    'size',{   1   ,    1   ,    2   ,    4   ,    8   ,    4    ,    8    });
c = t(1);
b = t(2);
s = t(3);
i = t(4);
l = t(5);
f = t(6);
d = t(7);

table = {...
    i,  1, 'sizeof_hdr', 540
    c,  8, 'magic', ['ni2' char(0) sprintf('\r\n\032\n')]
    s,  1, 'datatype', 2
    s,  1, 'bitpix', 8
    l,  8, 'dim', [3 0 0 0 1 1 1 1]
    d,  1, 'intent_p1', 0
    d,  1, 'intent_p2', 0
    d,  1, 'intent_p3', 0
    d,  8, 'pixdim', [0 1 1 1]
    l,  1, 'vox_offset', 0
    d,  1, 'scl_slope', 1
    d,  1, 'scl_inter', 0
    d,  1, 'cal_max', []
    d,  1, 'cal_min', []
    d,  1, 'slice_duration', []
    d,  1, 'toffset', []
    l,  1, 'slice_start', []
    l,  1, 'slice_end', []
    c, 80, 'descrip', 'NIFTI-2 Image'
    c, 24, 'aux_file', ''
    i,  1, 'qform_code', 0
    i,  1, 'sform_code', 0
    d,  1, 'quatern_b', 0
    d,  1, 'quatern_c', 0
    d,  1, 'quatern_d', 0
    d,  1, 'qoffset_x', 0
    d,  1, 'qoffset_y', 0
    d,  1, 'qoffset_z', 0
    d,  4, 'srow_x', [1 0 0 0]
    d,  4, 'srow_y', [0 1 0 0]
    d,  4, 'srow_z', [0 0 1 0]
    i,  1, 'slice_code', []
    i,  1, 'xyzt_units', 10
    i,  1, 'intent_code', 0
    c, 16, 'intent_name', ''
    b,  1, 'dim_info', []
    c, 15, 'unused_str', ''};


org = struct('label',table(:,3),'dtype',table(:,1),'len',table(:,2),...
    'offset',0,'def',table(:,4));
os  = 0;
for j=1:length(org)
    os  = org(j).dtype.size*ceil(os/org(j).dtype.size);
    fun = org(j).dtype.conv;
    if ischar(org(j).def), z = char(0); else z = 0; end
    def = [org(j).def repmat(z,1,org(j).len-length(org(j).def))];
    org(j).def    = feval(fun,def);
    org(j).offset = os;
    os  = os + org(j).len*org(j).dtype.size;
end
o = org;
