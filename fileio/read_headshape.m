function [shape] = read_headshape(filename, varargin)

% READ_HEADSHAPE reads the fiducials and/or the measured headshape
% from a variety of files (like CTF and Polhemus). The headshape and
% fiducials can for example be used for coregistration.
%
% Use as
%   [shape] = read_headshape(filename)
%
% See also READ_VOL, READ_SENS

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: read_headshape.m,v $
% Revision 1.13  2009/09/26 10:46:25  vlalit
% Added 4d_pdf to the list of BTi formats
%
% Revision 1.12  2009/04/01 16:59:43  vlalit
% Slight fix for JMs fix to assign unique names to the other BTi fiducials. Also make
%  sure that fid.label is a column cell array.
%
% Revision 1.11  2009/03/31 12:18:48  jansch
% added additional fiducials to shape.fid for 4D hs_files
%
% Revision 1.10  2009/03/23 12:09:07  vlalit
% Minor changes to make the ctf_old option fully functional.
%
% Revision 1.9  2009/03/17 10:58:13  vlalit
% Switched to MNE reader as default for Neuromag in read_headshape
%
% Revision 1.8  2009/01/28 18:29:07  vlalit
% Added '4d' type to BTi case
%
% Revision 1.7  2009/01/23 10:32:55  vlalit
% New reader for Neuromag fif format using the MNE toolbox (http://www.nmr.mgh.harvard.edu/martinos/userInfo/data/sofMNE.php)  implemented by Laurence Hunt.
%
% Revision 1.6  2008/10/07 16:22:32  roboos
% added option to specify coordinates to be obtained from ctf hc file
%
% Revision 1.5  2008/05/22 14:33:18  vlalit
% Changes related to generalization of fiducials'  handling in SPM.
%
% Revision 1.4  2008/05/11 16:30:30  vlalit
% Improved the support of 4d and neuromag
%
% Revision 1.3  2008/04/16 08:04:03  roboos
% allow headshape to be extracted from BEM volume conduction model
%
% Revision 1.2  2008/04/14 20:52:11  roboos
% ensure consistent output for all  file formats (thanks to Vladimir)
% added convert_units
%
% Revision 1.1  2008/04/11 12:04:55  roboos
% new impoementation, required for clean interface towards SPM
%

% test whether the file exists
if ~exist(filename)
    error(sprintf('file ''%s'' does not exist', filename));
end

% get the options
fileformat = keyval('fileformat',  varargin);
coordinates = keyval('coordinates',  varargin); if isempty(coordinates), coordinates = 'head'; end

if isempty(fileformat)
    fileformat = filetype(filename);
end

% start with an empty structure
shape           = [];
shape.pnt       = [];
shape.fid.pnt   = [];
shape.fid.label = {};

switch fileformat
    case {'ctf_ds', 'ctf_hc', 'ctf_meg4', 'ctf_res4', 'ctf_old'}
        [p, f, x] = fileparts(filename);
        
        if strcmp(fileformat, 'ctf_old')
            fileformat = filetype(filename);
        end
        
        if strcmp(fileformat, 'ctf_ds')
            filename = fullfile(p, [f x], [f '.hc']);
        elseif strcmp(fileformat, 'ctf_meg4')
            filename = fullfile(p, [f '.hc']);
        elseif strcmp(fileformat, 'ctf_res4')
            filename = fullfile(p, [f '.hc']);
        end

        orig = read_ctf_hc(filename);
        switch coordinates
          case 'head'
            shape.fid.pnt = cell2mat(struct2cell(orig.head));
          case 'dewar'
            shape.fid.pnt = cell2mat(struct2cell(orig.dewar));
          otherwise
            error('incorrect coordinates specified');
        end
        shape.fid.label = fieldnames(orig.head);

    case 'ctf_shape'
        orig = read_ctf_shape(filename);
        shape.pnt = orig.pnt;
        shape.fid.label = {'NASION', 'LEFT_EAR', 'RIGHT_EAR'};
        for i = 1:numel(shape.fid.label)
            shape.fid.pnt = cat(1, shape.fid.pnt, ...
                getfield(orig.MRI_Info, shape.fid.label{i}));
        end

    case {'4d_xyz', '4d_m4d', '4d_hs', '4d', '4d_pdf'}
        [p, f, x] = fileparts(filename);
        if ~strcmp(fileformat, '4d_hs')
            filename = fullfile(p, 'hs_file');
        end
        [shape.pnt, fid] = read_bti_hs(filename);
        
        % I'm making some assumptions here
        % which I'm not sure will work on all 4D systems
        
        %fid = fid(1:3, :);
        
        [junk, NZ] = max(fid(1:3,1));
        [junk, L]  = max(fid(1:3,2));
        [junk, R]  = min(fid(1:3,2));
        rest       = setdiff(1:size(fid,1),[NZ L R]);

        shape.fid.pnt = fid([NZ L R rest], :);
        shape.fid.label = {'NZ', 'L', 'R'};
        if ~isempty(rest),
            for i = 4:size(fid,1)
                shape.fid.label{i} = ['fiducial' num2str(i)];
                %in a 5 coil configuration this corresponds with Cz and Inion
            end
        end
    case 'neuromag_mex'
        [co,ki,nu] = hpipoints(filename);
        fid = co(:,find(ki==1))';

        [junk, NZ] = max(fid(:,2));
        [junk, L] = min(fid(:,1));
        [junk, R] = max(fid(:,1));

        shape.fid.pnt = fid([NZ L R], :);
        shape.fid.label = {'NZ', 'L', 'R'};
        
    case {'neuromag_mne', 'neuromag_fif'}
        hdr = read_header(filename,'headerformat','neuromag_mne');
        nFid = size(hdr.orig.dig,2); %work out number of fiducials
        switch coordinates
            case 'head' % digitiser points should be stored in head coordinates by default
                
                fidN=1;
                pntN=1;
                for i=1:nFid %loop over fiducials
                    %check this point is in head coordinates:
                    if hdr.orig.dig(i).coord_frame~=4 % 4 is MNE constant for head coordinates
                        error(['Digitiser point (' num2str(i) ') not stored in head coordinates!']);
                    end
                    
                    
                    switch hdr.orig.dig(i).kind % constants defined in MNE - see p.215 of MNE manual
                        case 1 % Cardinal point (nasion, LPA or RPA)
                            %get location of fiducial:
                            shape.fid.pnt(fidN,1:3) = hdr.orig.dig(i).r*100; %multiply by 100 to convert to cm
                            switch hdr.orig.dig(i).ident
                                case 1 % LPA
                                    shape.fid.label{fidN} = 'LPA';
                                case 2 % nasion
                                    shape.fid.label{fidN} = 'Nasion';
                                case 3 % RPA
                                    shape.fid.label{fidN} = 'RPA';
                                otherwise
                                    error('Unidentified cardinal point in file!');                                        
                            end
                            fidN = fidN + 1;
                                
                        case 2 % HPI coil
                            shape.pnt(pntN,1:3) = hdr.orig.dig(i).r*100;
                            pntN = pntN + 1;
                        case 3 % EEG electrode location (or ECG)
                            shape.pnt(pntN,1:3) = hdr.orig.dig(i).r*100;
                            pntN = pntN + 1;
                        case 4 % Additional head point
                            shape.pnt(pntN,1:3) = hdr.orig.dig(i).r*100;
                            pntN = pntN + 1;
                        otherwise
                            warning('Unidentified digitiser point in file!');
                    end
                        
                end
                shape.fid.label=shape.fid.label';
                
            case 'dewar'
                error('Dewar coordinates not supported for headshape yet (MNE toolbox)');
            otherwise
                error('Incorrect coordinates specified');
        end
        
    case 'polhemus_fil'
        [shape.fid.pnt, shape.pnt, shape.fid.label] = read_polhemus_fil(filename, 0);
        
    case 'spmeeg_mat'
        tmp = load(filename);
        if isfield(tmp.D, 'fiducials') && ~isempty(tmp.D.fiducials)
            shape = tmp.D.fiducials;
        else
            error('no headshape found in SPM EEG file');
        end

    case 'matlab'
        tmp = load(filename);
        if isfield(tmp, 'shape')
            shape = tmp.shape;
        elseif isfield(tmp, 'elec')
            shape.fid.pnt   = tmp.elec.pnt;
            shape.fid.label = tmp.elec.label;
        else
            error('no headshape found in Matlab file');
        end

    otherwise

        success = 0;
        if ~success
            % try reading it as electrode positions
            % and treat those as fiducials
            try
                elec = read_sens(filename);
                if ~senstype(elec, 'eeg')
                    error('headshape information can not be read from MEG gradiometer file');
                else
                    shape.fid.pnt   = elec.pnt;
                    shape.fid.label = elec.label;
                    success = 1;
                end
            end
        end

        if ~success
            % try reading it as volume conductor
            % and treat the skin surface as headshape
            try
                vol = read_vol(filename);
                if ~voltype(vol, 'bem')
                    error('skin surface can only be extracted from boundary element model');
                else
                    if ~isfield(vol, 'skin')
                        vol.skin = find_outermost_boundary(vol.bnd);
                    end
                    shape.pnt = vol.bnd(vol.skin).pnt;
                    shape.tri = vol.bnd(vol.skin).tri; % also return the triangulation
                    success = 1;
                end
            end
        end

        if ~success
            error('unknown fileformat for head shape information');
        end
end

shape.fid.label = shape.fid.label(:);

% this will add the units to the head shape
shape = convert_units(shape);
