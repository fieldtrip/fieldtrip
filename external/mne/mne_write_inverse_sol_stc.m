function mne_write_inverse_sol_stc(stem,inv,sol,tmin,tstep)
%
% function mne_write_inverse_sol_stc(stem,inv,sol,tmin,tstep)
%
% Save dynamic inverse solution data into stc files
%
% stem      - Stem for the stc files
% inv       - The inverse operator structure (can be the forward operator as well)
% sol       - A solution matrix (locations x time)
% tmin      - Time of the first data point in seconds
% tstep     - Time between data points in seconds
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.3  2009/03/04 03:47:15  msh
%   Fixed comments and added a note that the source spaces can come from a forward operator as well
%
%   Revision 1.2  2006/09/14 22:12:48  msh
%   Added output of the files written.
%
%   Revision 1.1  2006/05/05 03:50:40  msh
%   Added routines to compute L2-norm inverse solutions.
%   Added mne_write_inverse_sol_stc to write them in stc files
%   Several bug fixes in other files
%
%

me='MNE:mne_write_inverse_sol_stc';
global FIFF;
if isempty(FIFF)
   FIFF = fiff_define_constants();
end

if nargin ~= 5
    error(me,'Incorrect number of arguments');
end

if size(sol,1) ~= inv.nsource
    error(me,'The solution matrix cannot correspond to these source spaces');
end

off = 0;
for k = 1:length(inv.src)
    off = off + inv.src(k).nuse;
    if (inv.src(k).id ~= FIFF.FIFFV_MNE_SURF_LEFT_HEMI && ...
            inv.src(k).id ~= FIFF.FIFFV_MNE_SURF_RIGHT_HEMI)
        error(me,'Source space hemispheres not properly assigned.');
    end
end
if off ~= inv.nsource
    error(me,'The source spaces are inconsistent with other inverse/forward operator data');
end
%
%   Write a separate stc file for each source space
%   
off = 0;
for k = 1:length(inv.src)
    if (inv.src(k).id == FIFF.FIFFV_MNE_SURF_LEFT_HEMI)
        outname = sprintf('%s-lh.stc',stem);
    elseif (inv.src(k).id == FIFF.FIFFV_MNE_SURF_RIGHT_HEMI)
        outname = sprintf('%s-rh.stc',stem);
    end
    stc.tmin     = tmin;
    stc.tstep    = tstep;
    stc.vertices = inv.src(k).vertno - 1;
    stc.data     = sol(off+1:off+inv.src(k).nuse,:);
    mne_write_stc_file(outname,stc);
    off = off + inv.src(k).nuse;
    fprintf(1,'Wrote %s\n',outname);
end

return;



