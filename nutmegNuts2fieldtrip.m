function [raw,grid,mri]=nutmegNuts2fieldtrip(nuts)

% NUTMEGBEAM2FIELDTRIP converts the sensor data structure in NUTMEG called
% 'nuts' to a valid FieldTrip 'raw' structure
%
% [raw,grid,mri]=nutmegNuts2fieldtrip(nuts)
%
% nuts:        output from Nutmeg session file 
%              May be either 1) the file name (e.g. nuts*.mat) or 2) the nuts structure
%
% output:
%         raw: FT style 'raw' data structrue 
%        grid: FT style grid with source .pos and .leadfield
%         mri: FT style 'mri' from ft_read_mri with appropriate 'transform'
%              to the grid.pos

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

if ~isstruct(nuts) && exist(nuts,'file')
    nuts=load(nuts);
elseif ~isstruct(nuts)
    error('please enter either nuts structure or nuts filename')
end

mri=ft_read_mri(nuts.coreg.mripath);
mri.transform=0.1*nuts.coreg.meg2mri_tfm*mri.transform; % now .transform goes from MRI voxel to MEG cm, and matches .pos

grid.pos=0.1*nuts.voxels; % NM in CTF mm, FT in cm
grid.unit='cm';
grid.inside=1:size(nuts.voxels,1);
grid.outside=[];
if isfield(nuts,'Lp')
for kk=1:size(nuts.Lp,3)
    grid.leadfield{kk}=nuts.Lp(:,:,kk);
end
end

for kk=1:size(nuts.meg.data,3)
    raw.trial{kk}=nuts.meg.data(:,:,kk)';
    raw.time{kk}=0.001*nuts.meg.latency';
end
raw.label=nuts.meg.sensor_labels;

raw.grad.label=nuts.meg.sensor_labels;
raw.grad.ori=[nuts.meg.sensorOrient; nuts.meg.refSensorOrient];
raw.grad.ori=cat(1,raw.grad.ori,raw.grad.ori);
raw.grad.pnt=[nuts.meg.sensorCoord; nuts.meg.refSensorCoord];
raw.grad.pnt=reshape(permute(raw.grad.pnt,[1 3 2]),size(raw.grad.pnt,1)*2,3);
if isfield(nuts.meg,'chanmixMtx')
    raw.grad.tra=nuts.meg.chanmixMtx{1};
else
    raw.grad.tra=zeros(size(nuts.meg.data,2),2*size(nuts.meg.data,2)+size(nuts.meg.refSensorCoord,1));
    raw.grad.tra(:,1:size(nuts.meg.data,2))=eye(size(nuts.meg.data,2));
    raw.grad.tra(:,size(nuts.meg.data,2)+1:2*size(nuts.meg.data,2))=eye(size(nuts.meg.data,2));
    raw.grad.tra(:,2*size(nuts.meg.data,2)+1:end)=nuts.meg.Gcoef;
    warning('FIXME: should Gcoef be added regardless of nuts.meg.grad_order?')
end

if nuts.meg.grad_order==0
    raw.grad.balance.current='none';
elseif nuts.meg.grad_order==1
    raw.grad.balance.current='G1BR';
elseif nuts.meg.grad_order==2
    raw.grad.balance.current='G2BR';
elseif nuts.meg.grad_order==3
    raw.grad.balance.current='G3BR';
end    

raw.vol.o=nuts.meg.lsc;
% raw.vol.r=9*size(raw.vol.o); % to appease ft_sourceanalysis


