function [trans_head2mri] = fiff_read_coord_trans(transfile)
%
%  usage: [trans_head2mri] = fiff_read_coord_trans(transfile)
%
%  input:
%       transfile = name of transformation fif file (usually stored in
%       the subject's Freesurfer directory in /mri/T1-neuromag/sets).
%
%  output:
%       trans_head2mri = transformation structure from head to MRI coordinate systems.
%
%                    Note: the inverse transformation, from MRI to head coordinate systems
%                              can be obtained by just taking the inverse:
%                                   trans_mri2head.from=5; trans_mri2head.to=4;
%                                   trans_mri2head.trans=inv(trans_head2mri.trans);
%
%  author:  Rey Ramirez  email: rrramir@uw.edu

FIFF = fiff_define_constants;
[fid,~,dir] = fiff_open(transfile);
a = find([dir.kind] == FIFF.FIFF_COORD_TRANS);
pos = dir(a(1)).pos;
tag = fiff_read_tag(fid,pos);
trans_head2mri = tag.data;
