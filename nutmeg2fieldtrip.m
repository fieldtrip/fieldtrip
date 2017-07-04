function [data, mri, grid] = nutmeg2fieldtrip(cfg,fileorstruct)

% NUTMEG2FIELDTRIP converts from NUTMEG either a sensor data structure
% ('nuts') to a valid FieldTrip 'raw' structure (plus 'grid' and 'mri' if
% available), OR a source structure 'beam' to a valid FieldTrip source
% structure
%
% Use as
%    [data, mri, grid] = nutmeg2fieldtrip(cfg, fileorstruct)
%
% Input:
%      cfg
%          .keepmri (required for either input): =1 calls ft_read_mri for 'mri' output; =0 not save out 'mri'
%          .out     (required for source input): 's' (pos_freq_time) or 'trial' (pos_rpt)
%      fileorstruct: may be one of following:
%         1) *.mat file containing nuts sensor structure
%         2) nuts sensor structure
%         3) s*.mat file containing beam source structure
%         4) beam source structure (output from Nutmeg (beamforming_gui, tfbf, or tfZ)
%               (only scalar not vector results supported at the moment)
%
% Output: depending on input, one of options
%         1) If nuts sensor structure input, then 'data' will be 'raw' and
%            optionally 'grid' if Lp present, or 'mri' if individual MRI present
%         2) If beam source structure input, then 'data' will be 'source'
%            (May be an array of source structures (source{1} etc))
%            'grid' and 'mri' may be output as well if present in beam structure
%
% See alo FT_DATATYPE_SOURCE, LORETA2FIELDTRIP, SPASS2FIELDTRIP, FIELDTRIP2SPSS

% Copyright (C) 2011, Johanna Zumer
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble provenance
ft_preamble trackconfig

if ~isstruct(fileorstruct) && exist(fileorstruct,'file')
  structin=load(fileorstruct);
elseif isstruct(fileorstruct)
  structin=fileorstruct;
elseif ~isstruct(fileorstruct)
  ft_error('please enter either a structure (beam or nuts) or filename containing either beam or nuts')
end
clear fileorstruct

if isfield(structin,'meg')
  nutsorbeam=1;
elseif isfield(structin,'nuts')
  structin=structin.nuts;
  nutsorbeam=1;
elseif isfield(structin,'s')
  nutsorbeam=2;
elseif isfield(structin,'sa')
  structin=nut_beam_legacy_compatibility(structin);
  nutsorbeam=2;
elseif isfield(structin,'beam')
  structin=structin.beam;
  if isfield(structin,'sa')
    structin=nut_beam_legacy_compatibility(structin);
  end
  nutsorbeam=2;
else
  ft_error('not sure of the input type. maybe you need to call nut_beam_legacy_compatibility?')
end

if cfg.keepmri
  mri=ft_read_mri(structin.coreg.mripath);
  if isfield(mri,'transform')
    mri.transform=0.1*structin.coreg.meg2mri_tfm*mri.transform; % now .transform goes from MRI voxel to MEG cm, and matches .pos
  else
    % not yet sure what to do with lack of mri.transform
  end
else
  mri=[];
end

if isfield(structin,'voxels')
  grid.pos=0.1*structin.voxels; % NM in CTF mm, FT in cm
  grid.unit='cm';
  grid.inside=1:size(structin.voxels,1);
  grid.outside=[];
end
if isfield(structin,'Lp')
  for kk=1:size(structin.Lp,3)
    grid.leadfield{kk}=structin.Lp(:,:,kk);
  end
end
if isfield(structin,'W')
  grid.filter=structin.W;
end

if nutsorbeam==1
  for kk=1:size(structin.meg.data,3)
    raw.trial{kk}=structin.meg.data(:,:,kk)';
    raw.time{kk}=0.001*structin.meg.latency';
  end
  raw.label=structin.meg.sensor_labels;
  raw.grad.label=structin.meg.sensor_labels;
  if isfield(structin.meg,'refSensorOrient')
    raw.grad.coilori=[structin.meg.sensorOrient; structin.meg.refSensorOrient];
    raw.grad.coilpos=[structin.meg.sensorCoord; structin.meg.refSensorCoord];
  else
    raw.grad.coilori=[structin.meg.sensorOrient];
    raw.grad.coilpos=[structin.meg.sensorCoord];
  end
  if size(raw.grad.coilori,1)<2*size(structin.meg.data,2)
    raw.grad.coilori=cat(1,raw.grad.coilori,raw.grad.coilori);
  end
  raw.grad.coilpos=reshape(permute(raw.grad.coilpos,[1 3 2]),size(raw.grad.coilpos,1)*size(raw.grad.coilpos,3),size(raw.grad.coilpos,2));
  if isfield(structin.meg,'chanmixMtx')
    raw.grad.tra=structin.meg.chanmixMtx{1};
  else
    raw.grad.tra=zeros(size(structin.meg.data,2),2*size(structin.meg.data,2)+size(structin.meg.refSensorCoord,1));
    raw.grad.tra(:,1:size(structin.meg.data,2))=eye(size(structin.meg.data,2));
    raw.grad.tra(:,size(structin.meg.data,2)+1:2*size(structin.meg.data,2))=eye(size(structin.meg.data,2));
    raw.grad.tra(:,2*size(structin.meg.data,2)+1:end)=structin.meg.Gcoef;
    ft_warning('FIXME: should Gcoef be added regardless of structin.meg.grad_order?')
  end

  if structin.meg.grad_order==0
    raw.grad.balance.current='none';
  elseif structin.meg.grad_order==1
    raw.grad.balance.current='G1BR';
  elseif structin.meg.grad_order==2
    raw.grad.balance.current='G2BR';
  elseif structin.meg.grad_order==3
    raw.grad.balance.current='G3BR';
  end

  raw.headmodel.o = structin.meg.lsc; % note this may include reference channels. may need to do str_match of lsc_sensor_labels with sensor_labels?
  data=raw;
  clear raw

elseif nutsorbeam==2
  cfg.ftsource='old'; % 'old' way of source structure in FT
  tmp.pos=0.1*structin.voxels; % NM in CTF mm, FT in cm
  tmp.inside=1:size(structin.voxels,1);
  tmp.outside=[];
  if strcmp(cfg.out,'s')
    if isfield(structin,'bands')
      tmp.freq=mean(structin.bands'); % center of each band.
    end
    if isfield(structin,'timepts')
      tmp.time=0.001*structin.timepts'; % NM in ms, FT in s
    end
    if strcmp(cfg.ftsource,'old')
      % must have separate outputs per freq band
      if size(structin.s{1},4)==1
        for ii=1:size(structin.s,2)
          for jj=1:size(structin.s{1},3)
            source{ii,jj}=tmp;
            source{ii,jj}.freq=source{ii,jj}.freq(jj);
            for vv=1:size(structin.s{ii},1)
              source{ii,jj}.avg(vv).pow=structin.s{ii}(vv,:,jj);
            end
          end
        end
      elseif size(structin.s{1},4)>1
        for ii=1:size(structin.s,2)
          source{ii}=tmp;
          for kk=1:size(structin.s{ii},4)
            for vv=1:size(structin.s{ii},1)
              source{ii}.avg(vv).mom(kk,:)=structin.s{ii}(vv,:,1,kk);
            end
          end
        end
      end
    elseif strcmp(cfg.ftsource,'new')
      if size(structin.s{1},2)>1 && size(structin.s{1},3)>1
        tmp.powdimord='pos_freq_time'; % FIXME is this valid dimord?
      elseif size(structin.s{1},3)>1
        tmp.powdimord='pos_freq';
      elseif size(structin.s{1},4)>1
        tmp.powdimord='pos_ori_time';
      elseif size(structin.s{1},2)>1
        tmp.powdimord='pos_time';
      else
        tmp.powdimord='pos';
      end
      if size(structin.s{1},2)>1
        if size(structin.s{1},4)==1
          for ii=1:size(structin.s,2)
            source{ii}=tmp;
            source{ii}.pow=squeeze(permute(structin.s{ii},[1 3 2]));
          end
        elseif size(structin.s{1},4)==1
          for ii=1:size(structin.s,2)
            source{ii}=tmp;
            source{ii}.pow=squeeze(permute(structin.s{ii},[1 4 3 2]));
          end
        end
      end
    end
  elseif strcmp(cfg.out,'trial')
    if isfield(structin,'trial') % from tfZ keeptrials='yes'
      for ii=1:size(structin.trial,2)
        source{ii}=tmp;
        if strcmp(cfg.ftsource,'old')
          for jj=1:size(structin.trial{ii},2)
            source{ii}.trial(jj).pow=structin.trial{ii}(:,jj);
          end
        elseif strcmp(cfg.ftsource,'new')
          source{ii}.powdimord='pos_rpt'; % FIXME is this a valid field?
          source{ii}.pow=structin.trial{ii};
        end
      end
    end
  else
    ft_error('not a valid cfg.out specified');
  end
  data=source;
  clear source
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
% save the output cfg in all three output data structures
ft_postamble history    data
ft_postamble history    mri
ft_postamble history    grid
