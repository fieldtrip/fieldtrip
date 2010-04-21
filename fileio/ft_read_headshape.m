function [shape] = ft_read_headshape(filename, varargin)

% FT_READ_HEADSHAPE reads the fiducials and/or the measured headshape
% from a variety of files (like CTF and Polhemus). The headshape and
% fiducials can for example be used for coregistration.
%
% Use as
%   [shape] = ft_read_headshape(filename)
%
% See also FT_READ_VOL, FT_READ_SENS

% Copyright (C) 2008-2010 Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% test whether the file exists
if ~exist(filename)
    error(sprintf('file ''%s'' does not exist', filename));
end

% get the options
fileformat = keyval('fileformat',  varargin);
coordinates = keyval('coordinates',  varargin); if isempty(coordinates), coordinates = 'head'; end

if isempty(fileformat)
    fileformat = ft_filetype(filename);
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
            fileformat = ft_filetype(filename);
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
        hdr = ft_read_header(filename,'headerformat','neuromag_mne');
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
        
    case 'yokogawa_coregis'
        
        in_str = textread(filename, '%s');
        nr_items = size(in_str,1);
        ind = 1;
        coil_ind = 1;
        while ind < nr_items
            if strcmp(in_str{ind},'MEG:x=')
                shape.fid.pnt = [shape.fid.pnt; str2num(strtok(in_str{ind+1},[',','['])) ...
                    str2num(strtok(in_str{ind+3},[',','['])) str2num(strtok(in_str{ind+5},[',','[']))];
                shape.fid.label = [shape.fid.label ; ['Marker',num2str(coil_ind)]];
                coil_ind = coil_ind + 1;
                ind = ind + 6;
            else
                ind = ind +1;
            end
        end
        if size(shape.fid.label,1) ~= 5
            error('Wrong number of coils');           
        end        
        
        sw_ind = [3 1 2];
        
        shape.fid.pnt(1:3,:)= shape.fid.pnt(sw_ind, :);
        shape.fid.label(1:3)= {'nas', 'lpa', 'rpa'};        
        
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
                elec = ft_read_sens(filename);
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
                vol = ft_read_vol(filename);
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
shape = ft_convert_units(shape);
