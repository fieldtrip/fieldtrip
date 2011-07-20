function source = nutmegBeam2fieldtripSource(cfg,beam)

% NUTMEGBEAM2FIELDTRIPSOURCE converts the source structure in NUTMEG called
% 'beam' to a valid FieldTrip source structure
%
% source = nutmegBeam2fieldtripSource(beam)
%
% beam:        output from Nutmeg (beamforming_gui, tfbf, or tfZ)
%              only scalar (not vector) results supported at the moment
%              May be either 1) the file name (s_beam*.mat) or 2) the beam structure
% cfg.
%    .out         's' (pos_freq_time) or 'trial' (pos_rpt)
%    .keepvol     =1 calls ft_read_mri into .vol field; =0 leave out
%
% source output will be array of structures (source{1}, etc)
%
% see ft_datatype_source of current FieldTrip version to see if any modifications

% Copyright (C) 2011, Johanna Zumer
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

cfg.ftsource='old';


if ~isstruct(beam) && exist(beam,'file')
    beam=load(beam);
elseif ~isstruct(beam)
    error('please enter either beam structure or beam filename')
end

tmp.pos=0.1*beam.voxels; % NM in CTF mm, FT in cm
tmp.inside=1:size(beam.voxels,1);
tmp.outside=[];

if cfg.keepvol
    tmp.vol=ft_read_mri(beam.coreg.mripath);
    tmp.vol.transform=0.1*beam.coreg.meg2mri_tfm*tmp.vol.transform; % now .transform goes from MRI voxel to MEG cm, and matches .pos
end

if strcmp(cfg.out,'s')
    if isfield(beam,'bands')
        tmp.freq=mean(beam.bands'); % center of each band.
    end
    if isfield(beam,'timepts')
        tmp.time=0.001*beam.timepts'; % NM in ms, FT in s
    end
    if strcmp(cfg.ftsource,'old')
        % must have separate outputs per freq band
        if size(beam.s{1},4)==1
            for ii=1:size(beam.s,2)
                for jj=1:size(beam.s{1},3)
                    source{ii,jj}=tmp;
                    source{ii,jj}.freq=source{ii,jj}.freq(jj);
                    for vv=1:size(beam.s{ii},1)
                        source{ii,jj}.avg(vv).pow=beam.s{ii}(vv,:,jj);
                    end
                end
            end
        elseif size(beam.s{1},4)>1
            for ii=1:size(beam.s,2)
                source{ii}=tmp;
                for kk=1:size(beam.s{ii},4)
                    for vv=1:size(beam.s{ii},1)
                        source{ii}.avg(vv).mom(jj,:)=beam.s{ii}(vv,:,1,jj);
                    end
                end
            end
        end
    elseif strcmp(cfg.ftsource,'new')
        if size(beam.s{1},2)>1 && size(beam.s{1},3)>1
            tmp.powdimord='pos_freq_time'; % FIXME is this valid dimord?
        elseif size(beam.s{1},3)>1
            tmp.powdimord='pos_freq';
        elseif size(beam.s{1},4)>1
            tmp.powdimord='pos_ori_time';
        elseif size(beam.s{1},2)>1
            tmp.powdimord='pos_time';
        else
            tmp.powdimord='pos';
        end
        if size(beam.s{1},2)>1
            if size(beam.s{1},4)==1
                for ii=1:size(beam.s,2)
                    source{ii}=tmp;
                    source{ii}.pow=squeeze(permute(beam.s{ii},[1 3 2]));
                end
            elseif size(beam.s{1},4)==1
                for ii=1:size(beam.s,2)
                    source{ii}=tmp;
                    source{ii}.pow=squeeze(permute(beam.s{ii},[1 4 3 2]));
                end
            end
        end
    end
elseif strcmp(cfg.out,'trial')
    if isfield(beam,'trial') % from tfZ keeptrials='yes'
        for ii=1:size(beam.trial,2)
            source{ii}=tmp;
            if strcmp(cfg.ftsource,'old')
                for jj=1:size(beam.trial{ii},2)
                    source{ii}.trial(jj).pow=beam.trial{ii}(:,jj);
                end
            elseif strcmp(cfg.ftsource,'new')
                source{ii}.powdimord='pos_rpt'; % FIXME is this a valid field?
                source{ii}.pow=beam.trial{ii};
            end
        end
    end
else
    error('not a valid cfg.out specified');
end
