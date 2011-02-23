function [res] = mne_ex_read_evoked(fname,setno,apply_proj,dest_comp,use_ctf_head)
%
%   Load one evoked-response data set and do various kinds
%   of preprocessing
%
%   [res] = mne_ex_read_evoked(fname,setno,apply_proj,dest_comp,use_ctf_head)
%
%   fname           - Name of the data file
%   setno           - Data set number (default = 1)
%   apply_proj      - Apply SSP to the data (default = true)
%   dest_comp       - Desired (CTF/4D) compensation in the output data (default = keep the one in the file)
%   use_ctf_head    - Use the CTF/4D head coordinate system instead of the
%                     Neuromag one if available
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.3  2008/11/18 02:38:51  msh
%   Modified mne_ex_read_evoked to apply projection and compensation
%   Modified mne_ex_read_raw to call mne_set_current_comp
%
%   Revision 1.2  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.1  2006/04/21 17:31:07  msh
%   Added the examples.
%   Modified the formatting of some informative output.
%
%

me='MNE:mne_ex_read_evoked';

%
%  Different cases for the input arguments
%
if nargin == 1
    setno = 1;
    apply_proj   = true;
    use_ctf_head = false;
    keep_comp    = true;
elseif nargin == 2
    apply_proj   = true;
    use_ctf_head = false;
    keep_comp    = true;
elseif nargin == 3
    use_ctf_head = false;
    keep_comp    = true;
elseif nargin == 4
    use_ctf_head = false;
    keep_comp    = false;
elseif nargin == 5
    keep_comp     = false;
else
    error(me,'Incorrect number of arguments');
end
%
%   Load the desired data
%
fprintf(1,'\nLoading set %d from %s...\n\n',setno,fname);
%
try
    res = fiff_read_evoked(fname,setno);
catch
    error(me,'%s',mne_omit_first_line(lasterr));
end
%
%   Go ahead and do all kinds of interesting stuff
%
fprintf(1,'\nPreprocessing...\n');
%
%   Omit bad channels and pick MEG
%   (Change this according to your needs)
%
want_meg   = true;
want_eeg   = false;
want_stim  = false;
%
%include{1} = 'STI 014';
include = [];
%
res = fiff_pick_types_evoked(res,want_meg,want_eeg,want_stim,include,res.info.bads);
fprintf(1,'\t%d channels remain after picking\n',res.info.nchan);
%
%   Handle the software gradient compensation
%
comp=[];
current_comp = mne_get_current_comp(res.info);
if current_comp > 0
    fprintf(1,'\tCurrent compensation grade : %d\n',current_comp);
end
if keep_comp
    dest_comp = current_comp;
end
if current_comp ~= dest_comp
    try
        comp = mne_make_compensator(res.info,current_comp,dest_comp);
        res.info.chs  = mne_set_current_comp(res.info.chs,dest_comp);
        fprintf(1,'\tAppropriate compensator created to change to grade %d.\n',dest_comp);
    catch
        error(me,'%s',mne_omit_first_line(lasterr));
    end
end
if ~isempty(comp)
    for k = 1:length(res.evoked)
        res.evoked(k).epochs = comp*res.evoked(k).epochs;
    end
    fprintf(1,'\tThe data are now compensated to grade %d.\n',dest_comp);
end
%
%   Set up projection
%
if apply_proj
    if isempty(res.info.projs)
        fprintf(1,'\tNo projector specified in these data\n');
    else
        %
        %   Activate the projection items
        %
        for k = 1:length(res.info.projs)
            res.info.projs(k).active = true;
        end
        fprintf(1,'\t%d projection items activated\n',length(res.info.projs));
        %
        %   Create the projector
        %
        [res.info.proj,nproj] = mne_make_projector_info(res.info);
        if nproj == 0
            fprintf(1,'\tThe projection vectors do not apply to these channels\n');
        else
            fprintf(1,'\tCreated an SSP operator (subspace dimension = %d)\n',nproj);
            for k = 1:length(res.evoked)
                res.evoked(k).epochs = res.info.proj*res.evoked(k).epochs;
            end
            fprintf(1,'\tSSP operator applied to the evoked data\n');
        end
    end
end
%
%   Select the head coordinate system
%
if use_ctf_head
    if isempty(res.info.dev_ctf_t)
        error(me,'No CTF head transformation available');
    end
    meg_trans = res.info.dev_ctf_t;
    eeg_trans = fiff_invert_transform(res.info.ctf_head_t);
    fprintf(1,'\tEmploying the CTF/4D head coordinate system\n');
else
    meg_trans = res.info.dev_head_t;
    eeg_trans = [];
    fprintf(1,'\tEmploying the Neuromag head coordinate system\n');
end
%
%   Transform coil and electrode locations to the desired coordinate frame
%
res.info.chs = fiff_transform_meg_chs(res.info.chs,meg_trans);
res.info.chs = fiff_transform_eeg_chs(res.info.chs,eeg_trans);
%
%   Create the coil definitions
%
try
    accuracy = 1;       %   Use accuracy = 2 for accurate coil definitions
    res.info.chs = mne_add_coil_defs(res.info.chs,accuracy);
catch
    fprintf(1,'\tCoil definitions not added\n');
end
%
%   N.B. If a nonstandard (in MNE sense) coil def file is used, do
%
if false
    coil_def_file = 'whatever';
    templates = mne_load_coil_def(coil_def_file);
    res.info.chs = mne_add_coil_defs(res.info.chs,accuracy,templates);
end

fprintf(1,'\nReady.\n\n');

return;

end
