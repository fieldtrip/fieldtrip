function [data,times] = mne_ex_read_raw(fname,from,to,in_samples,dest_comp)
%
%   Example of reading raw data
%
%   [ data, times ] = mne_ex_read_raw(fname,from,to,in_samples,dest_comp);
%
%   data        - The data read, compensated and projected, channel by
%                 channel
%   times       - The time points of the samples, in seconds
%
%
%   fname       - The name of the input file
%   from        - Starting time or sample
%   to          - Ending time or sample
%   in_samples  - Are from and to given in samples rather than in seconds
%                 (optional)
%   dest_comp   - Desired (CTF) compensation in the output data (optional)
%
%   NOTE: The purpose of this function is to demonstrate the raw data reading
%   routines. In real world, you probably make multiple calls to
%   fiff_read_raw_segment_times or fiff_read_raw_segment
%   between open and close
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.5  2008/11/18 02:38:51  msh
%   Modified mne_ex_read_evoked to apply projection and compensation
%   Modified mne_ex_read_raw to call mne_set_current_comp
%
%   Revision 1.4  2006/05/16 00:39:32  msh
%   Fixed error in mne_ex_read_raw: All projection items were not activated.
%   List the initial states of projection items when they are loaded.
%
%   Revision 1.3  2006/05/05 03:50:40  msh
%   Added routines to compute L2-norm inverse solutions.
%   Added mne_write_inverse_sol_stc to write them in stc files
%   Several bug fixes in other files
%
%   Revision 1.2  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.1  2006/04/21 17:31:07  msh
%   Added the examples.
%   Modified the formatting of some informative output.
%
%

%
%   Fiddle with the arguments
%
me='MNE:mne_ex_read_raw';

keep_comp = false;
if nargin == 3
    in_samples = false;
    keep_comp = true;
elseif nargin == 4
    keep_comp = true;
elseif nargin ~= 5
    error(me,'Incorrect number of arguments');
end
%
%   Setup for reading the raw data
%
try
    raw = fiff_setup_read_raw(fname);
catch
    error(me,'%s',mne_omit_first_line(lasterr));
end
%
%   Set up pick list: MEG + STI 014 - bad channels
%
include{1} = 'STI 014';
want_meg   = true;
want_eeg   = false;
want_stim  = false;
%
picks = fiff_pick_types(raw.info,want_meg,want_eeg,want_stim,include,raw.info.bads);
%
%   Set up projection
%
if isempty(raw.info.projs)
    fprintf(1,'No projector specified for these data\n');
    raw.proj = [];
else
    %
    %   Activate the projection items
    %
    for k = 1:length(raw.info.projs)
        raw.info.projs(k).active = true;
    end
    fprintf(1,'%d projection items activated\n',length(raw.info.projs));
    %
    %   Create the projector
    %
    [proj,nproj] = mne_make_projector_info(raw.info);
    if nproj == 0
        fprintf(1,'The projection vectors do not apply to these channels\n');
        raw.proj = [];
    else
        fprintf(1,'Created an SSP operator (subspace dimension = %d)\n',nproj);
        raw.proj = proj;
    end
end
%
%   Set up the CTF compensator
%
current_comp = mne_get_current_comp(raw.info);
if current_comp > 0
    fprintf(1,'Current compensation grade : %d\n',current_comp);
end
if keep_comp
    dest_comp = current_comp;
end
if current_comp ~= dest_comp
    try
        raw.comp = mne_make_compensator(raw.info,current_comp,dest_comp);
        raw.info.chs  = mne_set_current_comp(raw.info.chs,dest_comp);
        fprintf(1,'Appropriate compensator added to change to grade %d.\n',dest_comp);
    catch
        error(me,'%s',mne_omit_first_line(lasterr));
    end
end
%
%   Read a data segment
%   times output argument is optional
%
try
    if in_samples
        [ data, times ] = fiff_read_raw_segment(raw,from,to,picks);
    else
        [ data, times ] = fiff_read_raw_segment_times(raw,from,to,picks);
    end
catch
    fclose(raw.fid);
    error(me,'%s',mne_omit_first_line(lasterr));
end
fprintf(1,'Read %d samples.\n',size(data,2));

return;

end
