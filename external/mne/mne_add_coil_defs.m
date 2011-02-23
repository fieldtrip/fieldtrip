function [res] = mne_add_coil_defs(chs,accuracy,templates)
%
% [res] = mne_add_coil_defs(chs,accuracy,coil_def_templates)
%
% Add transformed coil definitions to the channel info
%
% chs        - original channel definitions
% accuracy   - desired accuracy (0, 1, or 2, defaults to 1)
% templates  - coil definition templates
%              (defaults to $MNE_ROOT/setup/mne/coil_def.dat or $MNE_ROOT/share/mne/coil_def.dat)
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%   Revision 1.8  2006/05/03 18:53:05  msh
%   Approaching Matlab 6.5 backward compatibility
%
%   Revision 1.7  2006/04/26 00:14:24  msh
%   Return empty structure from fiff_pick_channels evoked and fiff_pick_types_evoked if there is no match.
%   Accuracy checking was incorrect in mne_add_coil_defs
%
%   Revision 1.6  2006/04/24 00:07:24  msh
%   Fixed error in the computation of the surface-based coordinate system.
%
%   Revision 1.5  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.4  2006/04/21 17:31:07  msh
%   Added the examples.
%   Modified the formatting of some informative output.
%
%   Revision 1.3  2006/04/18 20:44:46  msh
%   Added reading of forward solution.
%   Use length instead of size when appropriate
%
%   Revision 1.2  2006/04/17 15:01:34  msh
%   More small improvements.
%
%   Revision 1.1  2006/04/17 11:52:15  msh
%   Added coil definition stuff
%
%
me='MNE:mne_add_coil_defs';

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

if nargin == 2 || nargin == 1 || isempty(templates)
    %
    %   We need the templates
    %
    try
        templates = mne_load_coil_def;
    catch
        error(me,'%s',mne_omit_first_line(lasterr));
    end
    if nargin == 1
        accuracy = 1;
    end
elseif nargin ~= 3
    error(me,'Incorrect number of arguments');
end

fprintf(1,'\t%d coil definition templates available\n',length(templates));

if accuracy ~= 0 && accuracy ~= 1 && accuracy ~= 2
    error(me,'Accuracy should attain one of the values 0, 1, 2');
end

res    = chs;
nchan  = length(chs);
ntemp  = length(templates);
nmeg   = 0;
neeg   = 0;
eeg_cm = zeros(1,3);

first = true;
for k = 2:nchan
    if chs(k).kind == FIFF.FIFFV_EEG_CH || chs(k).kind == FIFF.FIFFV_MEG_CH || ...
            chs(k).kind == FIFF.FIFFV_REF_MEG_CH
        if first
            coord_frame = chs(k).coord_frame;
            first = false;
        elseif chs(k).coord_frame ~= coord_frame
            error(me,'All coils and electrodes are not in the same coordinate system (%s)',chs(k).ch_name);
        end
    end
end

lower_half = hex2dec('FFFF');
for k = 1:nchan
    res(k).coil_def = [];
    if chs(k).kind == FIFF.FIFFV_EEG_CH && ~isempty(chs(k).eeg_loc)
        %
        %   Set up a "coil definition" for an EEG electrode
        %
        temp             = templates(1);             %   Just pick one
        temp.class       = 1000;                     %   EEG electrode class
        temp.id          = chs(k).coil_type;         %   Copy from the channel info
        temp.accuracy    = 1;                        %   Normal accuracy no matter what
        temp.num_points  = size(chs(k).eeg_loc,2);   %   Number of points is normally two
        temp.size        = 10e-3;                    %   Rather arbitrary
        temp.baseline    = 0;                        %   Baseline does not apply
        temp.description = 'EEG electrode';          %   Useful description
        temp.coildefs    = zeros(temp.num_points,7); %   Location data
        temp.coildefs(1,1)   = 1;                    %   The electrode location
        temp.coildefs(1,2:4) = chs(k).eeg_loc(:,1)';
        nn  = chs(k).eeg_loc(:,1);                   %   The normal definition is probably immaterial
        len = nn'*nn;
        if len > 0
            temp.coildefs(1,5:7) = nn';
        end
        if temp.num_points > 1                       %   Add the reference electrode
            temp.coildefs(2,1)   = -1;
            temp.coildefs(2,2:4) = chs(k).eeg_loc(:,2);
            nn  = chs(k).eeg_loc(:,2);
            len = nn'*nn;
            if len > 0
                temp.coildefs(2,5:7) = nn';
            end
        end
        temp.FV          = [];                            %   Empty for now
        eeg_cm           = eeg_cm + temp.coildefs(1,2:4); %   Accumulate center of mass
        res(k).coil_def  = temp;                          %   There it is
        neeg = neeg + 1;
    elseif chs(k).kind == FIFF.FIFFV_MEG_CH || chs(k).kind == FIFF.FIFFV_REF_MEG_CH
        %
        %   Set up an MEG coil definition
        %
        coil_type = bitand(double(chs(k).coil_type),lower_half);
        temp = [];
        for t = 1:ntemp
            if templates(t).id == coil_type && templates(t).accuracy == accuracy
                temp = templates(t);
                break;
            end
        end
        %
        %  Did we find the template?
        %
        if isempty(temp)
            error(me,'Could not find an MEG coil template (coil type = %d accuracy = %d) for channel %s', ...
                coil_type,accuracy,chs(k).ch_name);
        end
        %
        %  Transform the template using the coil transformation
        %
        np = size(temp.coildefs,1);
        %
        %  Transform the integration points and the normals
        %
        trans = chs(k).coil_trans(1:3,:);
        temp.coildefs(:,2:4) = (trans*[ temp.coildefs(:,2:4)' ; ones(1,np) ])';
        temp.coildefs(:,5:7) = (trans*[ temp.coildefs(:,5:7)' ; zeros(1,np) ])';
        %
        %  Transform the vertices
        %
        if ~isempty(temp.FV)
            nvert = size(temp.FV.vertices,1);
            temp.FV.vertices = (trans*[temp.FV.vertices' ; ones(1,nvert)])';
        end
        %
        %   Attach the coil definition to the result
        %
        res(k).coil_def = temp;
        nmeg = nmeg + 1;
    end
end
%
%     Create the EEG coils drawings now that we have the center of mass
%
if neeg > 0
    eeg_cm = eeg_cm/neeg;
    for k = 1:nchan
        if res(k).kind == FIFF.FIFFV_EEG_CH && ~isempty(res(k).eeg_loc)
            res(k).coil_def.FV = draw_electrode(res(k).coil_def,eeg_cm);
        end
    end
end

fprintf(1,'\t%d MEG coil definitions and %d EEG electrodes set up\n',nmeg,neeg);

return;

    function FV = draw_electrode(def,cm)

        % create a patch for an EEG electrode

        if def.class == 1000
            % round disk
            Radius  = def.size/2;           % radius
            Len_cir = 15;                   % number of points for circle
            r0      = def.coildefs(1,2:4);
            %
            %   Make a template circle
            %
            rr(1,:) = cos(2*pi*[0:Len_cir]/Len_cir);
            rr(2,:) = sin(2*pi*[0:Len_cir]/Len_cir);
            rr(3,:) = 0.0;
            rr(4,:) = ones(1,Len_cir+1);
            %
            %    Scaling
            %
            T1 = [ [ Radius 0 0 0 ] ; [ 0 Radius 0 0 ] ; [ 0 0 1 0 ...
                ] ; [  0 0 0 1 ] ];
            %
            %   Rotation (this is only approximately correct...)
            %
            %   The normal is parallel to the vector from CM to the electrode location
            %
            nn     = (def.coildefs(1,2:4) - cm)';
            sizenn = norm(nn);
            if sizenn > 0
                nn = nn/sizenn;
                [ u, s, v ]  = svd(eye(3,3) - nn*nn');
                if nn'*u(:,3) < 0
                    u = - u;
                end
                T2 = [ u [ 0 0 0 ]' ; [ 0 0 0 1 ]];
            else
                T2 = eye(4,4);
            end
            %
            %   Translation
            %
            T3 = [ [ 1 0 0 r0(1) ] ; [ 0 1 0 r0(2) ] ; [ 0 0 1 r0(3) ...
                ] ; [  0 0 0 1 ] ];
            %
            %    Transfrom the data and create the structure
            %
            rr = T3*T2*T1*rr;
            FV.vertices = [[rr(1,:)' rr(2,:)' rr(3,:)']];
            FV.faces    = [1:length(FV.vertices)];
        end
        return;
    end

end
