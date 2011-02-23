function fiff_write_ctf_comp(fid,comps)
%
% fiff_write_ctf_comp(fid,comps)
%
% Writes the CTF compensation data into a fif file
%
%     fid           An open fif file descriptor
%     comps         The compensation data to write
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.6  2006/09/08 19:27:13  msh
%   Added KIT coil type to mne_load_coil_def
%   Allow reading of measurement info by specifying just a file name.
%
%   Revision 1.5  2006/06/22 21:22:46  msh
%   Take into account the possibility of calibrated compensation matrices
%
%   Revision 1.4  2006/05/03 19:09:03  msh
%   Fixed even more compatibility issues.
%
%   Revision 1.3  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.2  2006/04/18 20:44:46  msh
%   Added reading of forward solution.
%   Use length instead of size when appropriate
%
%   Revision 1.1  2006/04/12 10:29:02  msh
%   Made evoked data writing compatible with the structures returned in reading.
%
%

me='MNE:fiff_write_ctf_comp';

if nargin ~= 2
    error(me,'Incorrect number of arguments');
end

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

if isempty(comps)
    return;
end
%
%  This is very simple in fact
%
fiff_start_block(fid,FIFF.FIFFB_MNE_CTF_COMP);
for k = 1:length(comps)
    comp = comps(k);
    fiff_start_block(fid,FIFF.FIFFB_MNE_CTF_COMP_DATA);
    %
    %    Write the compensation kind
    %
    fiff_write_int(fid,FIFF.FIFF_MNE_CTF_COMP_KIND, ...
        comp.ctfkind);
    fiff_write_int(fid,FIFF.FIFF_MNE_CTF_COMP_CALIBRATED,comp.save_calibrated);
    %
    %    Write an uncalibrated or calibrated matrix
    %
    comp.data.data = inv(diag(comp.rowcals))*comp.data.data*inv(diag(comp.colcals));
    fiff_write_named_matrix(fid,FIFF.FIFF_MNE_CTF_COMP_DATA,comp.data);
    fiff_end_block(fid,FIFF.FIFFB_MNE_CTF_COMP_DATA);
end
fiff_end_block(fid,FIFF.FIFFB_MNE_CTF_COMP);

return;

end
