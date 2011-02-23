function [comp] = mne_make_compensator(info,from,to,exclude_comp_chs)
%
% [comp] = mne_make_compensator(info,from,to,exclude_comp_chs)
%
% info              - measurement info as returned by the fif reading routines
% from              - compensation in the input data
% to                - desired compensation in the output
% exclude_comp_chs  - exclude compensation channels from the output (optional)
%

%
% Create a compensation matrix to bring the data from one compensation
% state to another
%


%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%   Revision 1.8  2006/05/03 19:09:03  msh
%   Fixed even more compatibility issues.
%
%   Revision 1.7  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.6  2006/04/18 21:11:31  msh
%   Added option to omit compensation channels from the compensator
%
%   Revision 1.5  2006/04/18 20:44:46  msh
%   Added reading of forward solution.
%   Use length instead of size when appropriate
%
%   Revision 1.4  2006/04/17 11:52:15  msh
%   Added coil definition stuff
%
%   Revision 1.3  2006/04/15 12:21:00  msh
%   Several small improvements
%
%   Revision 1.2  2006/04/14 15:49:49  msh
%   Improved the channel selection code and added ch_names to measurement info.
%
%   Revision 1.1  2006/04/12 10:51:19  msh
%   Added projection writing and compensation routines
%
%

me='MNE:mne_make_compensator';

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

if nargin == 3
    exclude_comp_chs = false;
elseif nargin ~= 4
    error(me,'Incorrect number of arguments');
end

if from == to
    comp = zeros(info.nchan,info.nchan);
    return;
end

if from == 0
    C1 = zeros(info.nchan,info.nchan);
else
    try
        C1 = make_compensator(info,from);
    catch
        error(me,'Cannot create compensator C1 (%s)',mne_omit_first_line(lasterr));
    end
    if isempty(C1)
        error(me,'Desired compensation matrix (kind = %d) not found',from);
    end
end
if to == 0
    C2 = zeros(info.nchan,info.nchan);
else
    try
        C2 = make_compensator(info,to);
    catch
        error(me,'Cannot create compensator C2 (%s)',mne_omit_first_line(lasterr));
    end
    if isempty(C2)
        error(me,'Desired compensation matrix (kind = %d) not found',to);
    end
end
%
%   s_orig = s_from + C1*s_from = (I + C1)*s_from
%   s_to   = s_orig - C2*s_orig = (I - C2)*s_orig
%   s_to   = (I - C2)*(I + C1)*s_from = (I + C1 - C2 - C2*C1)*s_from
%
comp = eye(info.nchan,info.nchan) + C1 - C2 - C2*C1;

if exclude_comp_chs
    pick  = zeros(info.nchan);
    npick = 0;
    for k = 1:info.nchan
        if info.chs(k).kind ~= FIFF.FIFFV_REF_MEG_CH
            npick = npick + 1;
            pick(npick) = k;
        end
    end
    if npick == 0
        error(me,'Nothing remains after excluding the compensation channels');
    end
    comp = comp(pick(1:npick),:);
end

return;

    function this_comp = make_compensator(info,kind)

        for k = 1:length(info.comps)
            if info.comps(k).kind == kind
                this_data = info.comps(k).data;
                %
                %   Create the preselector
                %
                presel  = zeros(this_data.ncol,info.nchan);
                for col = 1:this_data.ncol
                    c = strmatch(this_data.col_names{col},info.ch_names,'exact');
                    if isempty(c)
                        error(me,'Channel %s is not available in data',this_data.col_names{col});
                    elseif length(c) > 1
                        error(me,'Ambiguous channel %s',mat.col_names{col});
                    end
                    presel(col,c) = 1.0;
                end
                %
                %   Create the postselector
                %
                postsel = zeros(info.nchan,this_data.nrow);
                for c = 1:info.nchan
                    row = strmatch(info.ch_names{c},this_data.row_names,'exact');
                    if length(row) > 1
                        error(me,'Ambiguous channel %s', info.ch_names{c});
                    elseif length(row) == 1
                        postsel(c,row) = 1.0;
                    end
                end
                this_comp = postsel*this_data.data*presel;
                return;
            end
        end
        this_comp = [];
        return;
    end

end
