function atlas = ft_read_atlas(filename, varargin)

% FT_READ_ATLAS reads an template/individual segmentation or parcellation from
% disk. The volumetric segmentation or the surface-based parcellation can
% either represent a fixed template atlas (eg. AAL or the Talairach Daemon),
% it can represent an individualized atlas (e.g. obtained from FreeSurfer) or
% it can represent an unlabeled parcellation obtained from the individual's
% DTi or resting state fMRI.
%
% Use as
%   atlas = ft_read_atlas(filename, ...)
% where the output atlas will be represented as structure according to
% FT_DATATYPE_SEGMENTATION or FT_DATATYPE_PARCELLATION.
%
% The "lines" and the "colorcube" colormaps are useful for plotting the
% different patches.
%
% See also FT_READ_MRI, FT_READ_HEADSHAPE

% Copyright (C) 2005-2012, Robert Oostenveld, Ingrid Nieuwenhuis
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


[p, f, x] = fileparts(filename);

if strcmp(f, 'TTatlas+tlrc')
  defaultformat = 'afni';
else
  defaultformat  = 'wfu';
end

% get the optional input arguments
atlasformat = ft_getopt(varargin, 'format', defaultformat);

switch atlasformat
  case 'afni'
    % check whether the required AFNI toolbox is available
    ft_hastoolbox('afni', 1);
    
    atlas = ft_read_mri(filename);
    
    % the AFNI atlas contains two volumes at 1mm resolution
    atlas.brick0 = atlas.anatomy(:,:,:,1);
    atlas.brick1 = atlas.anatomy(:,:,:,2);
    atlas        = rmfield(atlas, 'anatomy');
    atlas.coord  = 'tal';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % the following information is from https://afni.nimh.nih.gov/afni/doc/misc/ttatlas_tlrc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    atlas.descr.brick = [
      1
      1
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      1
      1
      1
      0
      0
      0
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      1
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      1
      0
      0
      ];
    
    atlas.descr.value = [
      68
      71
      20
      21
      22
      24
      25
      26
      27
      28
      29
      30
      31
      32
      33
      34
      35
      36
      37
      39
      40
      41
      42
      43
      44
      45
      46
      47
      48
      49
      50
      51
      52
      70
      72
      73
      74
      75
      76
      77
      124
      125
      126
      128
      129
      130
      131
      132
      133
      134
      135
      136
      137
      138
      144
      145
      151
      146
      147
      148
      149
      81
      82
      83
      84
      85
      86
      87
      88
      89
      90
      91
      93
      94
      95
      96
      97
      98
      99
      100
      101
      102
      103
      104
      105
      106
      107
      108
      109
      110
      111
      112
      113
      114
      115
      116
      117
      118
      119
      120
      121
      122
      123
      53
      54
      55
      56
      57
      58
      59
      60
      61
      62
      63
      66
      65
      127
      64
      67
      ];
    
    atlas.descr.name = {
      'Hippocampus'
      'Amygdala'
      'Posterior Cingulate'
      'Anterior Cingulate'
      'Subcallosal Gyrus'
      'Transverse Temporal Gyrus'
      'Uncus'
      'Rectal Gyrus'
      'Fusiform Gyrus'
      'Inferior Occipital Gyrus'
      'Inferior Temporal Gyrus'
      'Insula'
      'Parahippocampal Gyrus'
      'Lingual Gyrus'
      'Middle Occipital Gyrus'
      'Orbital Gyrus'
      'Middle Temporal Gyrus'
      'Superior Temporal Gyrus'
      'Superior Occipital Gyrus'
      'Inferior Frontal Gyrus'
      'Cuneus'
      'Angular Gyrus'
      'Supramarginal Gyrus'
      'Cingulate Gyrus'
      'Inferior Parietal Lobule'
      'Precuneus'
      'Superior Parietal Lobule'
      'Middle Frontal Gyrus'
      'Paracentral Lobule'
      'Postcentral Gyrus'
      'Precentral Gyrus'
      'Superior Frontal Gyrus'
      'Medial Frontal Gyrus'
      'Lentiform Nucleus'
      'Hypothalamus'
      'Red Nucleus'
      'Substantia Nigra'
      'Claustrum'
      'Thalamus'
      'Caudate'
      'Caudate Tail'
      'Caudate Body'
      'Caudate Head'
      'Ventral Anterior Nucleus'
      'Ventral Posterior Medial Nucleus'
      'Ventral Posterior Lateral Nucleus'
      'Medial Dorsal Nucleus'
      'Lateral Dorsal Nucleus'
      'Pulvinar'
      'Lateral Posterior Nucleus'
      'Ventral Lateral Nucleus'
      'Midline Nucleus'
      'Anterior Nucleus'
      'Mammillary Body'
      'Medial Globus Pallidus'
      'Lateral Globus Pallidus'
      'Putamen'
      'Nucleus Accumbens'
      'Medial Geniculum Body'
      'Lateral Geniculum Body'
      'Subthalamic Nucleus'
      'Brodmann area 1'
      'Brodmann area 2'
      'Brodmann area 3'
      'Brodmann area 4'
      'Brodmann area 5'
      'Brodmann area 6'
      'Brodmann area 7'
      'Brodmann area 8'
      'Brodmann area 9'
      'Brodmann area 10'
      'Brodmann area 11'
      'Brodmann area 13'
      'Brodmann area 17'
      'Brodmann area 18'
      'Brodmann area 19'
      'Brodmann area 20'
      'Brodmann area 21'
      'Brodmann area 22'
      'Brodmann area 23'
      'Brodmann area 24'
      'Brodmann area 25'
      'Brodmann area 27'
      'Brodmann area 28'
      'Brodmann area 29'
      'Brodmann area 30'
      'Brodmann area 31'
      'Brodmann area 32'
      'Brodmann area 33'
      'Brodmann area 34'
      'Brodmann area 35'
      'Brodmann area 36'
      'Brodmann area 37'
      'Brodmann area 38'
      'Brodmann area 39'
      'Brodmann area 40'
      'Brodmann area 41'
      'Brodmann area 42'
      'Brodmann area 43'
      'Brodmann area 44'
      'Brodmann area 45'
      'Brodmann area 46'
      'Brodmann area 47'
      'Uvula of Vermis'
      'Pyramis of Vermis'
      'Tuber of Vermis'
      'Declive of Vermis'
      'Culmen of Vermis'
      'Cerebellar Tonsil'
      'Inferior Semi-Lunar Lobule'
      'Fastigium'
      'Nodule'
      'Uvula'
      'Pyramis'
      'Culmen'
      'Declive'
      'Dentate'
      'Tuber'
      'Cerebellar Lingual'
      };
    
    % the following is basically the conversion from the 2005 format to the 2012 format
    sel0   = (atlas.descr.brick==0);
    label0 = atlas.descr.name(sel0);
    value0 = atlas.descr.value(sel0);
    % construct a new array with parcel or atlas values
    if numel(label0)<=intmax('uint8')
      new_brick0 = zeros(atlas.dim, 'uint8');
    elseif numel(label0)<=intmax('uint16')
      new_brick0 = zeros(atlas.dim, 'uint16');
    elseif numel(label0)<=intmax('uint32')
      new_brick0 = zeros(atlas.dim, 'uint32');
    else
      new_brick0 = zeros(atlas.dim);
    end
    for i=1:numel(label0)
      % replace the original values with numbers from 1 to N
      new_brick0(atlas.brick0==value0(i)) = i;
    end
    
    sel1   = (atlas.descr.brick==1);
    label1 = atlas.descr.name(sel1);
    value1 = atlas.descr.value(sel1);
    % construct a new array with parcel or atlas values
    if numel(label1)<=intmax('uint8')
      new_brick1 = zeros(atlas.dim, 'uint8');
    elseif numel(label1)<=intmax('uint16')
      new_brick1 = zeros(atlas.dim, 'uint16');
    elseif numel(label1)<=intmax('uint32')
      new_brick1 = zeros(atlas.dim, 'uint32');
    else
      new_brick1 = zeros(atlas.dim);
    end
    for i=1:numel(label1)
      % replace the original values with numbers from 1 to N
      new_brick1(atlas.brick1==value1(i)) = i;
    end
    
    atlas = rmfield(atlas, {'brick0', 'brick1', 'descr'});
    atlas.brick0 = new_brick0;
    atlas.brick0label = label0;
    atlas.brick1 = new_brick1;
    atlas.brick1label = label1;
    
  case 'wfu'
    
    atlas        = ft_read_mri(filename);
    brick0       = atlas.anatomy;
    atlas        = rmfield(atlas, 'anatomy');
    
    % FIXME the human WFU atlas contains a single volume at 2mm resolution in
    % MNI space, but the coordinates will not be guaranteed to be MNI
    % compatible, for example for the rhesus or mouse atlas.
    
    % atlas.coord  = 'mni';
    
    [p, f, x] = fileparts(filename);
    
    % this is a mat file that Ingrid apparently discovered somewhere
    % filename1 = fullfile(p, [f '_List.mat']);
    
    filename2 = fullfile(p, [f '.txt']);
    if ~exist(filename2, 'file')
      error('cannot locate %s', filename2);
    end
    
    % the download from http://fmri.wfubmc.edu comes with pairs of nii and txt files
    % the text file looks like this, with tabs between the colums
    % the section at the end with three times 191 does not always exist
    %
    % [ TD Labels]
    % 53	Angular Gyrus			191 191 191
    % 39	Anterior Cingulate		191 191 191
    % 36	Caudate				191 191 191
    % 27	Cerebellar Lingual		191 191 191
    % 2	Cerebellar Tonsil		191 191 191
    % 52	Cingulate Gyrus			191 191 191
    % ...
    
    
    fid = fopen(filename2);
    i = 1;
    value = [];
    label = {};
    while true
      tline = fgetl(fid);
      if ~ischar(tline), break, end
      % use TAB as deliniter
      [num, rem] = strtok(tline, 9);
      [str, rem] = strtok(rem, 9);
      num = str2double(num);
      if ~isnan(num)
        value(i) = num;
        label{i} = str;
        i = i+1;
      end % if num
    end % while
    fclose(fid);
    
    % this is loosely modeled after the AFNI implementation, hence the "brick0" naming
    if numel(label)<=intmax('uint8')
      new_brick0 = zeros(atlas.dim, 'uint8');
    elseif numel(label)<=intmax('uint16')
      new_brick0 = zeros(atlas.dim, 'uint16');
    elseif numel(label)<=intmax('uint32')
      new_brick0 = zeros(atlas.dim, 'uint32');
    else
      new_brick0 = zeros(atlas.dim);
    end
    for i=1:numel(label)
      % replace the original values with numbers from 1 to N
      new_brick0(brick0==value(i)) = i;
    end
    
    % replace the original brick with interspersed integers with one that contains contiguous integets
    atlas.brick0      = new_brick0;
    atlas.brick0label = label;
    
  otherwise
    error('unsupported atlas format %s', atlasformat);
end % case
