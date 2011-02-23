function [res,proj,comp] = mne_ex_cancel_noise(data,dest_comp)
%
%   Do projection and compensation as needed
%   Return the appropriate operators
%   
%   [res,proj,comp] = mne_ex_cancel_noise(data,dest_comp)
%
%   res     - Data after noise cancellation
%   proj    - The projection operator applied
%   comp    - The compensator which brings uncompensated data to the
%             desired compensation grade (will be useful in forward
%             calculations)
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.4  2006/05/05 03:50:40  msh
%   Added routines to compute L2-norm inverse solutions.
%   Added mne_write_inverse_sol_stc to write them in stc files
%   Several bug fixes in other files
%
%   Revision 1.3  2006/04/27 20:57:41  msh
%   Fixed typo on line 76 of mne_ex_cancel_noise
%
%   Revision 1.2  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.1  2006/04/21 17:31:07  msh
%   Added the examples.
%   Modified the formatting of some informative output.
%
%

me='MNE:mne_ex_cancel_noise';

if nargin == 1
    dest_comp = 0;
elseif nargin ~= 2
    error(me,'Incorrect number of arguments');
end
%
%   Compensate the data and make a compensator for forward modelling
%
comp = [];
proj = [];
comp_now = mne_get_current_comp(data.info);
if comp_now == dest_comp
    res = data;
else
    try 
        res  = mne_compensate_to(data,dest_comp);
        fprintf(1,'The data are now compensated to grade %d.\n',dest_comp);
    catch
        error(me,'%s',mne_omit_first_line(lasterr));
    end
end
if dest_comp > 0
    comp = mne_make_compensator(res.info,0,dest_comp);
    fprintf(1,'Appropriate forward operator compensator created.\n');
else
    fprintf(1,'No forward operator compensator needed.\n');
end
%
%   Do the projection
%
if isempty(data.info.projs)
    fprintf(1,'No projector included with these data\n');
else
    %
    %   Activate the projection items
    %
    for k = 1:length(res.info.projs)
        res.info.projs(k).active = true;
    end
    %
    %   Create the projector
    %
    [proj,nproj] = mne_make_projector_info(res.info);
    if nproj == 0
        fprintf(1,'The projection vectors do not apply to these channels\n');
        proj = [];
    else
        fprintf(1,'Created an SSP operator (subspace dimension = %d)\n',nproj);
        res.evoked.epochs = proj*res.evoked.epochs;
        fprintf(1,'Projector applied to the data\n');
    end
end

return;

end
