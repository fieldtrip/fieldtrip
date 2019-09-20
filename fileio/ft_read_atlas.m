function atlas = ft_read_atlas(filename, varargin)

% FT_READ_ATLAS reads an template/individual segmentation or parcellation from disk.
% The volumetric segmentation or the surface-based parcellation can either represent
% a template atlas (eg. AAL or the Talairach Daemon), it can represent an
% individualized atlas (e.g. obtained from FreeSurfer) or it can represent an
% unlabeled parcellation obtained from the individual's DTi or resting state fMRI.
%
% Use as
%   atlas = ft_read_atlas(filename, ...)
% or
%   atlas = ft_read_atlas({filenamelabels, filenamemesh}, ...)
%
% Additional options should be specified in key-value pairs and can include
%   'format'      = string, see below
%   'unit'        = string, e.g. 'mm' (default is to keep it in the native units of the file)
%   'map'         = string, 'maxprob' (default), or 'prob', for fsl-based
%                     atlases, providing either a probabilistic
%                     segmentation or a maximum a posterior probability map
%
% For individual surface-based atlases from FreeSurfer you should specify two
% filenames as a cell-array: the first points to the file that contains information
% with respect to the parcels' labels, the second points to the file that defines the
% mesh on which the parcellation is defined.
%
% The output atlas will be represented as structure according to
% FT_DATATYPE_SEGMENTATION or FT_DATATYPE_PARCELLATION.
%
% The "lines" and the "colorcube" colormaps are useful for plotting the different
% patches, for example using FT_PLOT_MESH.
%
% See also FT_READ_MRI, FT_READ_HEADSHAPE, FT_PREPARE_SOURCEMODEL, FT_SOURCEPARCELLATE, FT_PLOT_MESH

% Copyright (C) 2005-2019, Robert Oostenveld, Ingrid Nieuwenhuis,
% Jan-Mathijs Schoffelen, Arjen Stolk
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

% deal with multiple filenames
if isa(filename, 'cell')
  if numel(filename)==2
    filenamemesh = filename{2};
    filename     = filename{1};
  else
    ft_error('with multiple filenames, only 2 files are allowed');
  end % if precisely two input files
end % iscell

% optionally get the data from the URL and make a temporary local copy
filename = fetch_url(filename);

[p, f, x] = fileparts(filename);

if strcmp(f, 'TTatlas+tlrc')
  format = 'afni';
elseif strcmp(x, '.nii') && exist(fullfile(p, [f '.txt']), 'file')
  % This is a combination of nii+txt file, where the txt file may contain three columns like this
  %   FAG	Precentral_L	2001
  %   FAD	Precentral_R	2002
  %   F1G	Frontal_Sup_L	2101
  %   F1D	Frontal_Sup_R	2102
  %   F1OG	Frontal_Sup_Orb_L	2111
  %   ...
  % Alternatively, the txt file may contain a header line with the atlas name between square brackets
  % and then a variable number of column text info, where the first column
  % is the index, and the second column the label.
  labelfile = fullfile(p, [f '.txt']);
  fid = fopen_or_error(labelfile, 'rt');
  l1  = fgetl(fid);
  if strcmp(l1(1),'[') && strcmp(l1(end),']')
    format = 'aal_ext';
  elseif strcmp(l1,'Brainnetome Atlas')
    format= 'brainnetome';
  else
    format = 'aal';
  end
  fclose(fid);
elseif strcmp(x, '.mgz') && ~isempty(strfind(f, 'aparc')) || ~isempty(strfind(f, 'aseg'))
  % individual volume based segmentation from freesurfer
  format = 'freesurfer_volume';
elseif ft_filetype(filename, 'caret_label')
  % this is a gifti file that contains both the values for a set of
  % vertices as well as the labels.
  format = 'caret_label';
elseif ~isempty(strfind(filename, 'MPM'))
  % assume to be from the spm_anatomy toolbox
  format = 'spm_anatomy';
elseif strcmp(x, '.xml') && (isfolder(strtok(fullfile(p,f), '_')) || isfolder(strtok(fullfile(p,f), '-')))
  % fsl-format atlas, this is assumed to consist of an .xml file that
  % specifies the labels, as well as the filenames of the files with the actual data stored
  % in a directory with the of the strtok'ed (with '-' or '_') file name.
  format = 'fsl';
elseif strcmp(x, '.mat')
  format = 'mat';
elseif strcmp(x, '.nii') && ~isempty(strfind(f, 'Yeo2011_7Networks'))
  format = 'yeo7';
elseif strcmp(x, '.nii') && ~isempty(strfind(f, 'Yeo2011_17Networks'))
  format = 'yeo17';
else
  format = 'wfu';
end

% get the optional input arguments
fileformat = ft_getopt(varargin, 'format', format);
unit       = ft_getopt(varargin, 'unit');

switch fileformat
  case 'aal'
    labelfile = fullfile(p, [f '.txt']);
    fid = fopen_or_error(labelfile, 'rt');
    C = textscan(fid, '%s%s%d');
    lab = C{2};
    idx = C{3};
    fclose(fid);
    
    atlas = ft_read_mri(filename);
    atlas.tissue = atlas.anatomy;
    atlas = rmfield(atlas, 'anatomy');
    atlas.tissuelabel       = {};
    atlas.tissuelabel(idx)  = lab;
    atlas.coordsys = 'mni';
    % The original contains a rather sparse labeling, since not all indices
    % are being used (it starts at 2001) The question is whether it is more
    % important to keep the original numbers or to make the list with
    % labels compact. This could be made optional.
    compact = true;
    if compact
      [a, i, j] = unique(atlas.tissue);
      atlas.tissue = reshape(j-1, atlas.dim);
      atlas.tissuelabel = atlas.tissuelabel(a(a~=0));
    end
    
  case 'aal_ext'
    labelfile = fullfile(p, [f '.txt']);
    fid = fopen_or_error(labelfile, 'rt');
    C = textscan(fid, '%d%s%*[^\n]', 'HeaderLines', 1, 'Delimiter', '\t');
    lab = C{2};
    idx = C{1};
    fclose(fid);
    
    atlas = ft_read_mri(filename);
    atlas.tissue = atlas.anatomy;
    atlas = rmfield(atlas, 'anatomy');
    atlas.tissuelabel       = {};
    atlas.tissuelabel(idx)  = lab;
    atlas.coordsys = 'mni';
    % The original contains a rather sparse labeling, since not all indices
    % are being used (it starts at 2001) The question is whether it is more
    % important to keep the original numbers or to make the list with
    % labels compact. This could be made optional.
    compact = true;
    if compact
      [a, i, j] = unique(atlas.tissue);
      atlas.tissue = reshape(j-1, atlas.dim);
      atlas.tissuelabel = atlas.tissuelabel(a(a~=0));
    end
    
  case 'brainnetome'
    % Brainnetome Atlas
    % L. Fan, et al.The Human Brainnetome Atlas: A New Brain Atlas Based on
    % Connectional Architecture. Cereb Cortex 2016; 26 (8): 3508-3526.
    % doi: 10.1093/cercor/bhw157
    atlas = ft_read_mri(filename);
    atlas.tissue = atlas.anatomy;
    atlas = rmfield(atlas, 'anatomy');
    atlas.coordsys = 'mni';
    
    %Brainnetome atlas comes as radiological view convention.
    %change to neurological view: patient Left->image Left.
    atlas.transform(1,1)=-atlas.transform(1,1);
    atlas.transform(1,4)=-atlas.transform(1,4);
    
    %labels
    atlas.tissuelabel = cell(1,246);
    fid = fopen_or_error(labelfile, 'rt');
    lab  = fgetl(fid); %lab='Brainnetome Atlas'
    for label_i=1:246
      atlas.tissuelabel{1,label_i}=fgetl(fid);
    end
    fclose(fid);
    
  case 'afni'
    % check whether the required AFNI toolbox is available
    ft_hastoolbox('afni', 1);
    
    atlas = ft_read_mri(filename);
    
    % the AFNI atlas contains two volumes at 1mm resolution
    atlas.brick0   = atlas.anatomy(:,:,:,1);
    atlas.brick1   = atlas.anatomy(:,:,:,2);
    atlas          = rmfield(atlas, 'anatomy');
    atlas.dim      = atlas.dim([1 2 3]);
    atlas.coordsys = 'tal';
    
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
    atlas.brick0      = new_brick0;
    atlas.brick0label = label0;
    atlas.brick1      = new_brick1;
    atlas.brick1label = label1;
    
  case 'wfu'
    
    atlas        = ft_read_mri(filename);
    brick0       = atlas.anatomy;
    atlas        = rmfield(atlas, 'anatomy');
    
    % FIXME the human WFU atlas contains a single volume at 2mm resolution in
    % MNI space, but the coordinates will not be guaranteed to be MNI
    % compatible, for example for the rhesus or mouse atlas.
    
    % atlas.coordsys  = 'mni';
    atlas.coordsys = 'unknown';
    
    [p, f, x] = fileparts(filename);
    
    % if the original file was a .gz
    if isequal(x,'.gz')
      [p, f, x] = fileparts(filename(1:end-3));
    end
    
    % this is a mat file that Ingrid apparently discovered somewhere
    % filename1 = fullfile(p, [f '_List.mat']);
    
    filename2 = fullfile(p, [f '.txt']);
    if exist(filename2, 'file')
      
      % the download from http://fmri.wfubmc.edu comes with pairs of nii and txt files
      % the text file looks like this, with tabs between the colums
      % the section at the end with three times 191 does not always exist
      %
      % [ TD Labels]
      % 53	Angular Gyrus         191 191 191
      % 39	Anterior Cingulate		191 191 191
      % 36	Caudate               191 191 191
      % 27	Cerebellar Lingual		191 191 191
      % 2	  Cerebellar Tonsil		  191 191 191
      % 52	Cingulate Gyrus			  191 191 191
      % ...
      
      
      fid = fopen_or_error(filename2);
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
      
    else
      % the file does not exist
      ft_warning('cannot locate %s, making fake tissue labels', filename2);
      value = [];
      label = {};
      for i=1:max(brick0(:))
        value(i) = i;
        label{i} = sprintf('tissue %d', i);
      end
    end
    
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
    atlas.parcellation      = new_brick0;
    atlas.parcellationlabel = label(:);
    
  case {'freesurfer_volume'}
    % numeric values in the volume correspond to a label that can be found
    % in a lookup table. The lookup table is provided here and has been
    % copied from FreeSurferColorLUT.txt (v5.1)
    
    
    % Below is the color table for the cortical labels of the seg volume
    % created by mri_aparc2aseg (with --a2009s flag) in which the aseg
    % cortex label is replaced by the labels in the aparc.a2009s. The
    % cortical labels are the same as in Simple_surface_labels2009.txt,
    % except that left hemisphere has 11100 added to the index and the
    % right has 12100 added.  The label names are also prepended with
    % ctx_lh_, ctx_rh_, wm_lh_ and wm_rh_ (note usage of _ instead of -
    % to differentiate from a2005s labels).
    
    
    value =[
      1
      2
      3
      4
      5
      6
      7
      8
      9
      10
      11
      12
      13
      14
      15
      16
      17
      18
      19
      20
      21
      22
      23
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
      38
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
      64
      65
      66
      67
      68
      69
      70
      71
      72
      73
      74
      75
      76
      77
      78
      79
      80
      81
      82
      83
      84
      85
      86
      96
      97
      98
      250
      251
      252
      253
      254
      255
      1000
      1001
      1002
      1003
      1004
      1005
      1006
      1007
      1008
      1009
      1010
      1011
      1012
      1013
      1014
      1015
      1016
      1017
      1018
      1019
      1020
      1021
      1022
      1023
      1024
      1025
      1026
      1027
      1028
      1029
      1030
      1031
      1032
      1033
      1034
      1035
      1100
      1101
      1102
      1103
      1104
      1200
      1201
      1202
      1205
      1206
      1207
      1210
      1211
      1212
      1105
      1106
      1107
      1108
      1109
      1110
      1111
      1112
      1113
      1114
      1115
      1116
      1117
      1118
      1119
      1120
      1121
      1122
      1123
      1124
      1125
      1126
      1127
      1128
      1129
      1130
      1131
      1132
      1133
      1134
      1135
      1136
      1137
      1138
      1139
      1140
      1141
      1142
      1143
      1144
      1145
      1146
      1147
      1148
      1149
      1150
      1151
      1152
      1153
      1154
      1155
      1156
      1157
      1158
      1159
      1160
      1161
      1162
      1163
      1164
      1165
      1166
      1167
      1168
      1169
      1170
      1171
      1172
      1173
      1174
      1175
      1176
      1177
      1178
      1179
      1180
      1181
      2000
      2001
      2002
      2003
      2004
      2005
      2006
      2007
      2008
      2009
      2010
      2011
      2012
      2013
      2014
      2015
      2016
      2017
      2018
      2019
      2020
      2021
      2022
      2023
      2024
      2025
      2026
      2027
      2028
      2029
      2030
      2031
      2032
      2033
      2034
      2035
      2100
      2101
      2102
      2103
      2104
      2105
      2106
      2107
      2108
      2109
      2110
      2111
      2112
      2113
      2114
      2115
      2116
      2117
      2118
      2119
      2120
      2121
      2122
      2123
      2124
      2125
      2126
      2127
      2128
      2129
      2130
      2131
      2132
      2133
      2134
      2135
      2136
      2137
      2138
      2139
      2140
      2141
      2142
      2143
      2144
      2145
      2146
      2147
      2148
      2149
      2150
      2151
      2152
      2153
      2154
      2155
      2156
      2157
      2158
      2159
      2160
      2161
      2162
      2163
      2164
      2165
      2166
      2167
      2168
      2169
      2170
      2171
      2172
      2173
      2174
      2175
      2176
      2177
      2178
      2179
      2180
      2181
      2200
      2201
      2202
      2205
      2206
      2207
      2210
      2211
      2212
      11100
      11101
      11102
      11103
      11104
      11105
      11106
      11107
      11108
      11109
      11110
      11111
      11112
      11113
      11114
      11115
      11116
      11117
      11118
      11119
      11120
      11121
      11122
      11123
      11124
      11125
      11126
      11127
      11128
      11129
      11130
      11131
      11132
      11133
      11134
      11135
      11136
      11137
      11138
      11139
      11140
      11141
      11142
      11143
      11144
      11145
      11146
      11147
      11148
      11149
      11150
      11151
      11152
      11153
      11154
      11155
      11156
      11157
      11158
      11159
      11160
      11161
      11162
      11163
      11164
      11165
      11166
      11167
      11168
      11169
      11170
      11171
      11172
      11173
      11174
      11175
      12100
      12101
      12102
      12103
      12104
      12105
      12106
      12107
      12108
      12109
      12110
      12111
      12112
      12113
      12114
      12115
      12116
      12117
      12118
      12119
      12120
      12121
      12122
      12123
      12124
      12125
      12126
      12127
      12128
      12129
      12130
      12131
      12132
      12133
      12134
      12135
      12136
      12137
      12138
      12139
      12140
      12141
      12142
      12143
      12144
      12145
      12146
      12147
      12148
      12149
      12150
      12151
      12152
      12153
      12154
      12155
      12156
      12157
      12158
      12159
      12160
      12161
      12162
      12163
      12164
      12165
      12166
      12167
      12168
      12169
      12170
      12171
      12172
      12173
      12174
      12175];
    
    label = {
      'Left-Cerebral-Exterior'
      'Left-Cerebral-White-Matter'
      'Left-Cerebral-Cortex'
      'Left-Lateral-Ventricle'
      'Left-Inf-Lat-Vent'
      'Left-Cerebellum-Exterior'
      'Left-Cerebellum-White-Matter'
      'Left-Cerebellum-Cortex'
      'Left-Thalamus'
      'Left-Thalamus-Proper'
      'Left-Caudate'
      'Left-Putamen'
      'Left-Pallidum'
      '3rd-Ventricle'
      '4th-Ventricle'
      'Brain-Stem'
      'Left-Hippocampus'
      'Left-Amygdala'
      'Left-Insula'
      'Left-Operculum'
      'Line-1'
      'Line-2'
      'Line-3'
      'CSF'
      'Left-Lesion'
      'Left-Accumbens-area'
      'Left-Substancia-Nigra'
      'Left-VentralDC'
      'Left-undetermined'
      'Left-vessel'
      'Left-choroid-plexus'
      'Left-F3orb'
      'Left-lOg'
      'Left-aOg'
      'Left-mOg'
      'Left-pOg'
      'Left-Stellate'
      'Left-Porg'
      'Left-Aorg'
      'Right-Cerebral-Exterior'
      'Right-Cerebral-White-Matter'
      'Right-Cerebral-Cortex'
      'Right-Lateral-Ventricle'
      'Right-Inf-Lat-Vent'
      'Right-Cerebellum-Exterior'
      'Right-Cerebellum-White-Matter'
      'Right-Cerebellum-Cortex'
      'Right-Thalamus'
      'Right-Thalamus-Proper'
      'Right-Caudate'
      'Right-Putamen'
      'Right-Pallidum'
      'Right-Hippocampus'
      'Right-Amygdala'
      'Right-Insula'
      'Right-Operculum'
      'Right-Lesion'
      'Right-Accumbens-area'
      'Right-Substancia-Nigra'
      'Right-VentralDC'
      'Right-undetermined'
      'Right-vessel'
      'Right-choroid-plexus'
      'Right-F3orb'
      'Right-lOg'
      'Right-aOg'
      'Right-mOg'
      'Right-pOg'
      'Right-Stellate'
      'Right-Porg'
      'Right-Aorg'
      '5th-Ventricle'
      'Left-Interior'
      'Right-Interior'
      'Left-Lateral-Ventricles'
      'Right-Lateral-Ventricles'
      'WM-hypointensities'
      'Left-WM-hypointensities'
      'Right-WM-hypointensities'
      'non-WM-hypointensities'
      'Left-non-WM-hypointensities'
      'Right-non-WM-hypointensities'
      'Left-F1'
      'Right-F1'
      'Optic-Chiasm'
      'Corpus_Callosum'
      'Left-Amygdala-Anterior'
      'Right-Amygdala-Anterior'
      'Dura'
      'Fornix'
      'CC_Posterior'
      'CC_Mid_Posterior'
      'CC_Central'
      'CC_Mid_Anterior'
      'CC_Anterior'
      'ctx-lh-unknown'
      'ctx-lh-bankssts'
      'ctx-lh-caudalanteriorcingulate'
      'ctx-lh-caudalmiddlefrontal'
      'ctx-lh-corpuscallosum'
      'ctx-lh-cuneus'
      'ctx-lh-entorhinal'
      'ctx-lh-fusiform'
      'ctx-lh-inferiorparietal'
      'ctx-lh-inferiortemporal'
      'ctx-lh-isthmuscingulate'
      'ctx-lh-lateraloccipital'
      'ctx-lh-lateralorbitofrontal'
      'ctx-lh-lingual'
      'ctx-lh-medialorbitofrontal'
      'ctx-lh-middletemporal'
      'ctx-lh-parahippocampal'
      'ctx-lh-paracentral'
      'ctx-lh-parsopercularis'
      'ctx-lh-parsorbitalis'
      'ctx-lh-parstriangularis'
      'ctx-lh-pericalcarine'
      'ctx-lh-postcentral'
      'ctx-lh-posteriorcingulate'
      'ctx-lh-precentral'
      'ctx-lh-precuneus'
      'ctx-lh-rostralanteriorcingulate'
      'ctx-lh-rostralmiddlefrontal'
      'ctx-lh-superiorfrontal'
      'ctx-lh-superiorparietal'
      'ctx-lh-superiortemporal'
      'ctx-lh-supramarginal'
      'ctx-lh-frontalpole'
      'ctx-lh-temporalpole'
      'ctx-lh-transversetemporal'
      'ctx-lh-insula'
      'ctx-lh-Unknown'
      'ctx-lh-Corpus_callosum'
      'ctx-lh-G_and_S_Insula_ONLY_AVERAGE'
      'ctx-lh-G_cingulate-Isthmus'
      'ctx-lh-G_cingulate-Main_part'
      'ctx-lh-G_cingulate-caudal_ACC'
      'ctx-lh-G_cingulate-rostral_ACC'
      'ctx-lh-G_cingulate-posterior'
      'ctx-lh-S_cingulate-caudal_ACC'
      'ctx-lh-S_cingulate-rostral_ACC'
      'ctx-lh-S_cingulate-posterior'
      'ctx-lh-S_pericallosal-caudal'
      'ctx-lh-S_pericallosal-rostral'
      'ctx-lh-S_pericallosal-posterior'
      'ctx-lh-G_cuneus'
      'ctx-lh-G_frontal_inf-Opercular_part'
      'ctx-lh-G_frontal_inf-Orbital_part'
      'ctx-lh-G_frontal_inf-Triangular_part'
      'ctx-lh-G_frontal_middle'
      'ctx-lh-G_frontal_superior'
      'ctx-lh-G_frontomarginal'
      'ctx-lh-G_insular_long'
      'ctx-lh-G_insular_short'
      'ctx-lh-G_and_S_occipital_inferior'
      'ctx-lh-G_occipital_middle'
      'ctx-lh-G_occipital_superior'
      'ctx-lh-G_occipit-temp_lat-Or_fusiform'
      'ctx-lh-G_occipit-temp_med-Lingual_part'
      'ctx-lh-G_occipit-temp_med-Parahippocampal_part'
      'ctx-lh-G_orbital'
      'ctx-lh-G_paracentral'
      'ctx-lh-G_parietal_inferior-Angular_part'
      'ctx-lh-G_parietal_inferior-Supramarginal_part'
      'ctx-lh-G_parietal_superior'
      'ctx-lh-G_postcentral'
      'ctx-lh-G_precentral'
      'ctx-lh-G_precuneus'
      'ctx-lh-G_rectus'
      'ctx-lh-G_subcallosal'
      'ctx-lh-G_subcentral'
      'ctx-lh-G_temporal_inferior'
      'ctx-lh-G_temporal_middle'
      'ctx-lh-G_temp_sup-G_temp_transv_and_interm_S'
      'ctx-lh-G_temp_sup-Lateral_aspect'
      'ctx-lh-G_temp_sup-Planum_polare'
      'ctx-lh-G_temp_sup-Planum_tempolare'
      'ctx-lh-G_and_S_transverse_frontopolar'
      'ctx-lh-Lat_Fissure-ant_sgt-ramus_horizontal'
      'ctx-lh-Lat_Fissure-ant_sgt-ramus_vertical'
      'ctx-lh-Lat_Fissure-post_sgt'
      'ctx-lh-Medial_wall'
      'ctx-lh-Pole_occipital'
      'ctx-lh-Pole_temporal'
      'ctx-lh-S_calcarine'
      'ctx-lh-S_central'
      'ctx-lh-S_central_insula'
      'ctx-lh-S_cingulate-Main_part_and_Intracingulate'
      'ctx-lh-S_cingulate-Marginalis_part'
      'ctx-lh-S_circular_insula_anterior'
      'ctx-lh-S_circular_insula_inferior'
      'ctx-lh-S_circular_insula_superior'
      'ctx-lh-S_collateral_transverse_ant'
      'ctx-lh-S_collateral_transverse_post'
      'ctx-lh-S_frontal_inferior'
      'ctx-lh-S_frontal_middle'
      'ctx-lh-S_frontal_superior'
      'ctx-lh-S_frontomarginal'
      'ctx-lh-S_intermedius_primus-Jensen'
      'ctx-lh-S_intraparietal-and_Parietal_transverse'
      'ctx-lh-S_occipital_anterior'
      'ctx-lh-S_occipital_middle_and_Lunatus'
      'ctx-lh-S_occipital_superior_and_transversalis'
      'ctx-lh-S_occipito-temporal_lateral'
      'ctx-lh-S_occipito-temporal_medial_and_S_Lingual'
      'ctx-lh-S_orbital-H_shapped'
      'ctx-lh-S_orbital_lateral'
      'ctx-lh-S_orbital_medial-Or_olfactory'
      'ctx-lh-S_paracentral'
      'ctx-lh-S_parieto_occipital'
      'ctx-lh-S_pericallosal'
      'ctx-lh-S_postcentral'
      'ctx-lh-S_precentral-Inferior-part'
      'ctx-lh-S_precentral-Superior-part'
      'ctx-lh-S_subcentral_ant'
      'ctx-lh-S_subcentral_post'
      'ctx-lh-S_suborbital'
      'ctx-lh-S_subparietal'
      'ctx-lh-S_supracingulate'
      'ctx-lh-S_temporal_inferior'
      'ctx-lh-S_temporal_superior'
      'ctx-lh-S_temporal_transverse'
      'ctx-rh-unknown'
      'ctx-rh-bankssts'
      'ctx-rh-caudalanteriorcingulate'
      'ctx-rh-caudalmiddlefrontal'
      'ctx-rh-corpuscallosum'
      'ctx-rh-cuneus'
      'ctx-rh-entorhinal'
      'ctx-rh-fusiform'
      'ctx-rh-inferiorparietal'
      'ctx-rh-inferiortemporal'
      'ctx-rh-isthmuscingulate'
      'ctx-rh-lateraloccipital'
      'ctx-rh-lateralorbitofrontal'
      'ctx-rh-lingual'
      'ctx-rh-medialorbitofrontal'
      'ctx-rh-middletemporal'
      'ctx-rh-parahippocampal'
      'ctx-rh-paracentral'
      'ctx-rh-parsopercularis'
      'ctx-rh-parsorbitalis'
      'ctx-rh-parstriangularis'
      'ctx-rh-pericalcarine'
      'ctx-rh-postcentral'
      'ctx-rh-posteriorcingulate'
      'ctx-rh-precentral'
      'ctx-rh-precuneus'
      'ctx-rh-rostralanteriorcingulate'
      'ctx-rh-rostralmiddlefrontal'
      'ctx-rh-superiorfrontal'
      'ctx-rh-superiorparietal'
      'ctx-rh-superiortemporal'
      'ctx-rh-supramarginal'
      'ctx-rh-frontalpole'
      'ctx-rh-temporalpole'
      'ctx-rh-transversetemporal'
      'ctx-rh-insula'
      'ctx-rh-Unknown'
      'ctx-rh-Corpus_callosum'
      'ctx-rh-G_and_S_Insula_ONLY_AVERAGE'
      'ctx-rh-G_cingulate-Isthmus'
      'ctx-rh-G_cingulate-Main_part'
      'ctx-rh-G_cuneus'
      'ctx-rh-G_frontal_inf-Opercular_part'
      'ctx-rh-G_frontal_inf-Orbital_part'
      'ctx-rh-G_frontal_inf-Triangular_part'
      'ctx-rh-G_frontal_middle'
      'ctx-rh-G_frontal_superior'
      'ctx-rh-G_frontomarginal'
      'ctx-rh-G_insular_long'
      'ctx-rh-G_insular_short'
      'ctx-rh-G_and_S_occipital_inferior'
      'ctx-rh-G_occipital_middle'
      'ctx-rh-G_occipital_superior'
      'ctx-rh-G_occipit-temp_lat-Or_fusiform'
      'ctx-rh-G_occipit-temp_med-Lingual_part'
      'ctx-rh-G_occipit-temp_med-Parahippocampal_part'
      'ctx-rh-G_orbital'
      'ctx-rh-G_paracentral'
      'ctx-rh-G_parietal_inferior-Angular_part'
      'ctx-rh-G_parietal_inferior-Supramarginal_part'
      'ctx-rh-G_parietal_superior'
      'ctx-rh-G_postcentral'
      'ctx-rh-G_precentral'
      'ctx-rh-G_precuneus'
      'ctx-rh-G_rectus'
      'ctx-rh-G_subcallosal'
      'ctx-rh-G_subcentral'
      'ctx-rh-G_temporal_inferior'
      'ctx-rh-G_temporal_middle'
      'ctx-rh-G_temp_sup-G_temp_transv_and_interm_S'
      'ctx-rh-G_temp_sup-Lateral_aspect'
      'ctx-rh-G_temp_sup-Planum_polare'
      'ctx-rh-G_temp_sup-Planum_tempolare'
      'ctx-rh-G_and_S_transverse_frontopolar'
      'ctx-rh-Lat_Fissure-ant_sgt-ramus_horizontal'
      'ctx-rh-Lat_Fissure-ant_sgt-ramus_vertical'
      'ctx-rh-Lat_Fissure-post_sgt'
      'ctx-rh-Medial_wall'
      'ctx-rh-Pole_occipital'
      'ctx-rh-Pole_temporal'
      'ctx-rh-S_calcarine'
      'ctx-rh-S_central'
      'ctx-rh-S_central_insula'
      'ctx-rh-S_cingulate-Main_part_and_Intracingulate'
      'ctx-rh-S_cingulate-Marginalis_part'
      'ctx-rh-S_circular_insula_anterior'
      'ctx-rh-S_circular_insula_inferior'
      'ctx-rh-S_circular_insula_superior'
      'ctx-rh-S_collateral_transverse_ant'
      'ctx-rh-S_collateral_transverse_post'
      'ctx-rh-S_frontal_inferior'
      'ctx-rh-S_frontal_middle'
      'ctx-rh-S_frontal_superior'
      'ctx-rh-S_frontomarginal'
      'ctx-rh-S_intermedius_primus-Jensen'
      'ctx-rh-S_intraparietal-and_Parietal_transverse'
      'ctx-rh-S_occipital_anterior'
      'ctx-rh-S_occipital_middle_and_Lunatus'
      'ctx-rh-S_occipital_superior_and_transversalis'
      'ctx-rh-S_occipito-temporal_lateral'
      'ctx-rh-S_occipito-temporal_medial_and_S_Lingual'
      'ctx-rh-S_orbital-H_shapped'
      'ctx-rh-S_orbital_lateral'
      'ctx-rh-S_orbital_medial-Or_olfactory'
      'ctx-rh-S_paracentral'
      'ctx-rh-S_parieto_occipital'
      'ctx-rh-S_pericallosal'
      'ctx-rh-S_postcentral'
      'ctx-rh-S_precentral-Inferior-part'
      'ctx-rh-S_precentral-Superior-part'
      'ctx-rh-S_subcentral_ant'
      'ctx-rh-S_subcentral_post'
      'ctx-rh-S_suborbital'
      'ctx-rh-S_subparietal'
      'ctx-rh-S_supracingulate'
      'ctx-rh-S_temporal_inferior'
      'ctx-rh-S_temporal_superior'
      'ctx-rh-S_temporal_transverse'
      'ctx-rh-G_cingulate-caudal_ACC'
      'ctx-rh-G_cingulate-rostral_ACC'
      'ctx-rh-G_cingulate-posterior'
      'ctx-rh-S_cingulate-caudal_ACC'
      'ctx-rh-S_cingulate-rostral_ACC'
      'ctx-rh-S_cingulate-posterior'
      'ctx-rh-S_pericallosal-caudal'
      'ctx-rh-S_pericallosal-rostral'
      'ctx-rh-S_pericallosal-posterior'
      'ctx_lh_Unknown'
      'ctx_lh_G_and_S_frontomargin'
      'ctx_lh_G_and_S_occipital_inf'
      'ctx_lh_G_and_S_paracentral'
      'ctx_lh_G_and_S_subcentral'
      'ctx_lh_G_and_S_transv_frontopol'
      'ctx_lh_G_and_S_cingul-Ant'
      'ctx_lh_G_and_S_cingul-Mid-Ant'
      'ctx_lh_G_and_S_cingul-Mid-Post'
      'ctx_lh_G_cingul-Post-dorsal'
      'ctx_lh_G_cingul-Post-ventral'
      'ctx_lh_G_cuneus'
      'ctx_lh_G_front_inf-Opercular'
      'ctx_lh_G_front_inf-Orbital'
      'ctx_lh_G_front_inf-Triangul'
      'ctx_lh_G_front_middle'
      'ctx_lh_G_front_sup'
      'ctx_lh_G_Ins_lg_and_S_cent_ins'
      'ctx_lh_G_insular_short'
      'ctx_lh_G_occipital_middle'
      'ctx_lh_G_occipital_sup'
      'ctx_lh_G_oc-temp_lat-fusifor'
      'ctx_lh_G_oc-temp_med-Lingual'
      'ctx_lh_G_oc-temp_med-Parahip'
      'ctx_lh_G_orbital'
      'ctx_lh_G_pariet_inf-Angular'
      'ctx_lh_G_pariet_inf-Supramar'
      'ctx_lh_G_parietal_sup'
      'ctx_lh_G_postcentral'
      'ctx_lh_G_precentral'
      'ctx_lh_G_precuneus'
      'ctx_lh_G_rectus'
      'ctx_lh_G_subcallosal'
      'ctx_lh_G_temp_sup-G_T_transv'
      'ctx_lh_G_temp_sup-Lateral'
      'ctx_lh_G_temp_sup-Plan_polar'
      'ctx_lh_G_temp_sup-Plan_tempo'
      'ctx_lh_G_temporal_inf'
      'ctx_lh_G_temporal_middle'
      'ctx_lh_Lat_Fis-ant-Horizont'
      'ctx_lh_Lat_Fis-ant-Vertical'
      'ctx_lh_Lat_Fis-post'
      'ctx_lh_Medial_wall'
      'ctx_lh_Pole_occipital'
      'ctx_lh_Pole_temporal'
      'ctx_lh_S_calcarine'
      'ctx_lh_S_central'
      'ctx_lh_S_cingul-Marginalis'
      'ctx_lh_S_circular_insula_ant'
      'ctx_lh_S_circular_insula_inf'
      'ctx_lh_S_circular_insula_sup'
      'ctx_lh_S_collat_transv_ant'
      'ctx_lh_S_collat_transv_post'
      'ctx_lh_S_front_inf'
      'ctx_lh_S_front_middle'
      'ctx_lh_S_front_sup'
      'ctx_lh_S_interm_prim-Jensen'
      'ctx_lh_S_intrapariet_and_P_trans'
      'ctx_lh_S_oc_middle_and_Lunatus'
      'ctx_lh_S_oc_sup_and_transversal'
      'ctx_lh_S_occipital_ant'
      'ctx_lh_S_oc-temp_lat'
      'ctx_lh_S_oc-temp_med_and_Lingual'
      'ctx_lh_S_orbital_lateral'
      'ctx_lh_S_orbital_med-olfact'
      'ctx_lh_S_orbital-H_Shaped'
      'ctx_lh_S_parieto_occipital'
      'ctx_lh_S_pericallosal'
      'ctx_lh_S_postcentral'
      'ctx_lh_S_precentral-inf-part'
      'ctx_lh_S_precentral-sup-part'
      'ctx_lh_S_suborbital'
      'ctx_lh_S_subparietal'
      'ctx_lh_S_temporal_inf'
      'ctx_lh_S_temporal_sup'
      'ctx_lh_S_temporal_transverse'
      'ctx_rh_Unknown'
      'ctx_rh_G_and_S_frontomargin'
      'ctx_rh_G_and_S_occipital_inf'
      'ctx_rh_G_and_S_paracentral'
      'ctx_rh_G_and_S_subcentral'
      'ctx_rh_G_and_S_transv_frontopol'
      'ctx_rh_G_and_S_cingul-Ant'
      'ctx_rh_G_and_S_cingul-Mid-Ant'
      'ctx_rh_G_and_S_cingul-Mid-Post'
      'ctx_rh_G_cingul-Post-dorsal'
      'ctx_rh_G_cingul-Post-ventral'
      'ctx_rh_G_cuneus'
      'ctx_rh_G_front_inf-Opercular'
      'ctx_rh_G_front_inf-Orbital'
      'ctx_rh_G_front_inf-Triangul'
      'ctx_rh_G_front_middle'
      'ctx_rh_G_front_sup'
      'ctx_rh_G_Ins_lg_and_S_cent_ins'
      'ctx_rh_G_insular_short'
      'ctx_rh_G_occipital_middle'
      'ctx_rh_G_occipital_sup'
      'ctx_rh_G_oc-temp_lat-fusifor'
      'ctx_rh_G_oc-temp_med-Lingual'
      'ctx_rh_G_oc-temp_med-Parahip'
      'ctx_rh_G_orbital'
      'ctx_rh_G_pariet_inf-Angular'
      'ctx_rh_G_pariet_inf-Supramar'
      'ctx_rh_G_parietal_sup'
      'ctx_rh_G_postcentral'
      'ctx_rh_G_precentral'
      'ctx_rh_G_precuneus'
      'ctx_rh_G_rectus'
      'ctx_rh_G_subcallosal'
      'ctx_rh_G_temp_sup-G_T_transv'
      'ctx_rh_G_temp_sup-Lateral'
      'ctx_rh_G_temp_sup-Plan_polar'
      'ctx_rh_G_temp_sup-Plan_tempo'
      'ctx_rh_G_temporal_inf'
      'ctx_rh_G_temporal_middle'
      'ctx_rh_Lat_Fis-ant-Horizont'
      'ctx_rh_Lat_Fis-ant-Vertical'
      'ctx_rh_Lat_Fis-post'
      'ctx_rh_Medial_wall'
      'ctx_rh_Pole_occipital'
      'ctx_rh_Pole_temporal'
      'ctx_rh_S_calcarine'
      'ctx_rh_S_central'
      'ctx_rh_S_cingul-Marginalis'
      'ctx_rh_S_circular_insula_ant'
      'ctx_rh_S_circular_insula_inf'
      'ctx_rh_S_circular_insula_sup'
      'ctx_rh_S_collat_transv_ant'
      'ctx_rh_S_collat_transv_post'
      'ctx_rh_S_front_inf'
      'ctx_rh_S_front_middle'
      'ctx_rh_S_front_sup'
      'ctx_rh_S_interm_prim-Jensen'
      'ctx_rh_S_intrapariet_and_P_trans'
      'ctx_rh_S_oc_middle_and_Lunatus'
      'ctx_rh_S_oc_sup_and_transversal'
      'ctx_rh_S_occipital_ant'
      'ctx_rh_S_oc-temp_lat'
      'ctx_rh_S_oc-temp_med_and_Lingual'
      'ctx_rh_S_orbital_lateral'
      'ctx_rh_S_orbital_med-olfact'
      'ctx_rh_S_orbital-H_Shaped'
      'ctx_rh_S_parieto_occipital'
      'ctx_rh_S_pericallosal'
      'ctx_rh_S_postcentral'
      'ctx_rh_S_precentral-inf-part'
      'ctx_rh_S_precentral-sup-part'
      'ctx_rh_S_suborbital'
      'ctx_rh_S_subparietal'
      'ctx_rh_S_temporal_inf'
      'ctx_rh_S_temporal_sup'
      'ctx_rh_S_temporal_transverse'};
    
    % read in the volume
    atlas = ft_read_mri(filename);
    dat   = atlas.anatomy;
    atlas = rmfield(atlas, 'anatomy');
    
    % get the unique values to save time later on
    uval = unique(dat(:));
    sel  = find(ismember(value, uval));
    fprintf('subselecting %d labels from the total list of %d\n', numel(sel), numel(label));
    value = value(sel);
    label = label(sel);
    
    % remap the values in the data
    aparc = zeros(size(dat));
    cnt   = 0;
    for k = 1:numel(value)
      sel = dat==value(k);
      if sum(sel(:))
        cnt = cnt+1;
        fprintf('re-indexing label %s to a value of %d (was %d)\n', label{k}, cnt, value(k));
        aparc(sel)      = cnt;
        aparclabel{cnt,1} = label{k};
      end
    end
    atlas.aparc      = aparc;
    atlas.aparclabel = aparclabel;
    
  case {'freesurfer_a2009s' 'freesurfer_aparc' 'freesurfer_ba'}
    % ensure freesurfer on the path and get the info how to get from value to label
    ft_hastoolbox('freesurfer', 1);
    
    if strcmp(fileformat, 'freesurfer_a2009s')
      parcelfield = 'a2009s';
    elseif strcmp(fileformat, 'freesurfer_aparc')
      parcelfield = 'aparc';
    elseif strcmp(fileformat, 'freesurfer_ba')
      parcelfield = 'BA';
    else
      ft_error('unknown freesurfer parcellation method requested');
      %[index, label, rgb] = read_fscolorlut(lookuptable);
      %label = cellstr(label);
      %rgb = rand(length(label),3);
    end
    %rgb = rgb(:,1) + rgb(:,2)*256 + rgb(:,3)*256*256;
    
    % read the labels
    switch ft_filetype(filename)
      %case 'caret_label'
      %  p = gifti(filename);
      %  p = p.cdata;
      case 'freesurfer_annot'
        [v, p, c] = read_annotation(filename);
        
        label = c.struct_names;
        rgba  = c.table(:,1:4);
        rgb   = c.table(:,5); % compound value that is used for the indexing in vector p
        index = ((1:c.numEntries)-1)';
      otherwise
        ft_error('unsupported fileformat for parcel file');
    end
    
    switch ft_filetype(filenamemesh)
      %case {'caret_surf' 'gifti'}
      %  tmp = gifti(filenamemesh);
      %  bnd.pos = ft_warp_apply(tmp.mat, tmp.vertices);
      %  bnd.tri = tmp.faces;
      %  reindex = false;
      case 'freesurfer_triangle_binary'
        [pos, tri] = read_surf(filenamemesh);
        
        % ensure the triangles to be 1-indexed
        if min(tri(:))==0 && max(tri(:))==size(pos,1)-1
          tri = tri+1;
        end
        
        bnd.pos    = pos;
        bnd.tri    = tri;
        reindex    = true;
      otherwise
        ft_error('unsupported fileformat for surface mesh');
    end
    
    % check the number of vertices
    if size(bnd.pos,1) ~= numel(p)
      ft_error('the number of vertices in the mesh does not match the number of elements in the parcellation');
    end
    
    % reindex the parcels, if needed: I am not fully sure about this, but the caret
    % label files seem to have the stuff numbered with normal numbers, with unknown
    % being -1. assuming them to be in order;
    if reindex
      % this is then freesurfer convention, coding in rgb
      newp = zeros(size(p));
      for k = 1:numel(label)
        newp(p==rgb(k)) = index(k)+1;
      end
    else
      uniquep = unique(p);
      if uniquep(1)<0,
        p(p<0) = 0;
      end
      newp   = p;
    end
    atlas       = [];
    atlas.pos   = bnd.pos;
    atlas.tri   = bnd.tri;
    atlas.(parcelfield)            = newp;
    atlas.([parcelfield, 'label']) = label;
    atlas.rgba  = rgba;
    atlas       = ft_determine_units(atlas);
    
  case 'caret_label'
    ft_hastoolbox('gifti', 1);
    g = gifti(filename);
    
    rgba = [];
    if isfield(g, 'labels'),
      label = g.labels.name(:);
      key   = g.labels.key(:);
      if isfield(g.labels, 'rgba'),
        rgba = g.labels.rgba; % I'm not sure whether this always exists
      end
    else
      label = g.private.label.name(:);
      key   = g.private.label.key(:);
      if isfield(g.private.label, 'rgba')
        rgba = g.private.label.rgba; % I'm not sure whether this always exists
      end
    end
    
    %label = g.private.label.name; % provides the name of the parcel
    %key   = g.private.label.key;  % maps value to name
    
    % Store each column in cdata as an independent parcellation, because
    % each vertex can have multiple values in principle
    
    atlas = [];
    for k = 1:size(g.cdata,2)
      tmporig  = g.cdata(:,k);
      tmpnew   = nan(size(tmporig));
      tmplabel = cell(0,1);
      tmprgba  = zeros(0,4);
      cnt = 0;
      for m = 1:numel(label)
        sel = find(tmporig==key(m));
        if ~isempty(sel)
          cnt = cnt+1;
          if any(strcmp(tmplabel,deblank(label{m})))
            % one feature of a gifti can be that the same labels can exist (and
            % are treated as a different parcel) while they should be the same
            % parcel, i.e. when there's an empty space at the end of the label
            val = find(strcmp(tmplabel, deblank(label{m})));
          else
            % add as a new label
            tmplabel{end+1,1} = label{m};
            if ~isempty(rgba), tmprgba(end+1,:)  = rgba(m,:); end
            val = cnt;
          end
          tmpnew(tmporig==key(m)) = val;
        end
      end
      
      % there is some additional meta data that may be useful, but for now
      % stick to the rather uninformative parcellation1/2/3 etc.
      %
      % if strcmp(g.private.data{k}.metadata(1).name, 'Name')
      %   parcelfield = fixname(g.private.data{k}.metadata(1).value);
      % else
      %   ft_error('could not determine parcellation name');
      % end
      
      if size(g.cdata,2)>1
        parcelfield = ['parcellation' num2str(k)];
      else
        parcelfield = 'parcellation';
      end
      
      atlas.(parcelfield)           = tmpnew;
      atlas.([parcelfield 'label']) = tmplabel;
      if ~isempty(tmprgba), atlas.rgba = tmprgba; end
    end
    
    if exist('filenamemesh', 'var')
      tmp       = ft_read_headshape(filenamemesh);
      atlas.pos = tmp.pos;
      atlas.tri = tmp.tri;
      atlas     = ft_determine_units(atlas);
    elseif ~isfield(atlas, 'coordsys')
      atlas.coordsys = 'unknown';
    end
    
  case 'spm_anatomy'
    ft_hastoolbox('spm8up', 1);
    
    % load the map, this is assumed to be the struct-array MAP
    load(filename);
    [p,f,e]      = fileparts(filename);
    mrifilename  = fullfile(p,[strrep(f, '_MPM',''),'.img']);
    atlas        = ft_read_mri(mrifilename, 'dataformat', 'analyze_img');
    tissue       = round(atlas.anatomy); % I don't know why the values are non-integer
    atlas        = rmfield(atlas, 'anatomy');
    label        = {MAP.name}';
    idx          = [MAP.GV]';
    
    % check whether all labels are present
    if numel(intersect(idx,unique(tissue(:))))<numel(idx)
      fprintf('there are fewer labels in the volume than in the list\n');
    end
    
    % remap the values of the labels to run from 1-numel(idx)
    newtissue = zeros(size(tissue));
    for k = 1:numel(idx)
      newtissue(tissue==idx(k)) = k;
    end
    atlas.tissue      = newtissue;
    atlas.tissuelabel = label;
    atlas.coordsys    = 'spm'; % I think this is safe to assume
    
    clear tissue newtissue;
    
  case 'fsl'
    map = ft_getopt(varargin, 'map', 'maxprob');
    switch map
      case 'prob'
        imagefile = 'imagefile';
      case 'maxprob'
        imagefile = 'summaryimagefile';
      otherwise
        error('unknown map requested');
    end
    
    ft_hastoolbox('gifti', 1);
    hdr = xmltree(filename);
    hdr = convert(hdr);
    
    % get the full path
    [p, f , x]  = fileparts(filename);
    
    mrifilename = sprintf('%s.nii.gz', hdr.header.images{1}.(imagefile));
    if isequal(mrifilename(1), '/')
      mrifilename = mrifilename(2:end); % remove first backslash to avoid error if no full path was given
    end
    
    % this uses the thresholded image
    switch map
      case 'maxprob'
        atlas        = ft_read_mri(fullfile(p, mrifilename));
        atlas.tissue = atlas.anatomy;
        atlas        = rmfield(atlas, 'anatomy');
        atlas.tissuelabel = hdr.data.label(:);
        atlas.coordsys    = 'mni';
      case 'prob'
        tmp   = ft_read_mri(fullfile(p, mrifilename));
        atlas = removefields(tmp, 'anatomy');
        for m = 1:numel(hdr.data.label)
          fn = strrep(hdr.data.label{m}, ' ' , '_');
          fn = strrep(fn, '-', '_');
          fn = strrep(fn, '(', '');
          fn = strrep(fn, ')', '');
          
          atlas.(fn) = tmp.anatomy(:,:,:,m);
        end
        atlas.coordsys = 'mni';
    end
    
    
  case 'mat'
    tmp = load(filename);
    if isfield(tmp, 'Vertices') && isfield(tmp, 'Atlas')
      % this applies to BrainStorm *.mat files containing cortical meshes
      atlas.pos = tmp.Vertices;
      % copy some optional fields over with a new name
      atlas = copyfields(tmp, atlas, {'Faces', 'Curvature', 'SulciMap'});
      atlas = renamefields(atlas, {'Faces', 'Curvature', 'SulciMap'}, {'tri', 'curv', 'sulc'});
      for i=1:numel(tmp.Atlas)
        name   = fixname(tmp.Atlas(i).Name);
        nlabel = numel(tmp.Atlas(i).Scouts);
        label  = cell(1, nlabel);
        index  = zeros(size(atlas.pos,1), 1);
        for j=1:nlabel
          index(tmp.Atlas(i).Scouts(j).Vertices) = j;
          label{j} = tmp.Atlas(i).Scouts(j).Label;
        end
        atlas.( name         ) = index;
        atlas.([name 'label']) = label;
      end
    elseif isfield(tmp, 'atlas')
      % this applies to most FieldTrip *.mat files
      atlas = tmp.atlas;
    elseif numel(fieldnames(tmp))==1
      % just take whatever variable is contained in the file
      fn = fieldnames(tmp);
      atlas = tmp.(fn{1});
      if isstruct(atlas)
        ft_warning('assuming that the variable "%s" in "%s" represents the atlas', fn{1}, filename);
      else
        ft_error('cannot read atlas structure from "%s"', filename);
      end
    else
      ft_error('the mat-file %s does not contain a variable called ''atlas''',filename);
    end
    
  case 'yeo7'
    % this uses Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask_colin27.nii, which is
    % the 7 network parcellation from https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation_Yeo2011 
    % aligned to the colin27 template (skull-stripped version of single_subj_T1_1mm.nii) 
    % using AFNI's 3dQwarp and 3dNwarpApply
    atlas = ft_read_mri(filename);
    atlas.tissue = atlas.anatomy;
    atlas = rmfield(atlas, 'anatomy');
    atlas.tissuelabel = {
      '7Networks_1'
      '7Networks_2'
      '7Networks_3'
      '7Networks_4'
      '7Networks_5'
      '7Networks_6'
      '7Networks_7'
      };
    atlas.coordsys = 'mni';
    colors = [
      120 18 134;
      70 130 180;
      0 118 14;
      196 58 250;
      220 248 164;
      230 148 34;
      205 62 78
      ]; % not used
    
  case 'yeo17'
    % this uses Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask_colin27.nii, which is
    % the 17 network parcelation from https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation_Yeo2011 
    % aligned to the colin27 template (skull-stripped version of single_subj_T1_1mm.nii) 
    % using AFNI's 3dQwarp and 3dNwarpApply
    atlas = ft_read_mri(filename);
    atlas.tissue = atlas.anatomy;
    atlas = rmfield(atlas, 'anatomy');
    atlas.tissuelabel = {
      '17Networks_1'
      '17Networks_2'
      '17Networks_3'
      '17Networks_4'
      '17Networks_5'
      '17Networks_6'
      '17Networks_7'
      '17Networks_8'
      '17Networks_9'
      '17Networks_10'
      '17Networks_11'
      '17Networks_12'
      '17Networks_13'
      '17Networks_14'
      '17Networks_15'
      '17Networks_16'
      '17Networks_17'
      };
    atlas.coordsys = 'mni';
    colors = [
      120 18 134;
      255 0 0;
      70 130 180;
      42 204 164;
      74 155 60;
      0 118 14;
      196 58 250;
      255 152 213;
      220 248 164;
      122 135 50;
      119 140 176;
      230 148 34;
      135 50 74;
      12 48 255;
      0 0 130;
      255 255 0;
      205 62 78
      ]; % not used
    
  otherwise
    ft_error('unsupported format "%s"', fileformat);
    
end % switch fileformat


if ~isempty(unit)
  % ensure the atlas is in the desired units
  atlas = ft_convert_units(atlas, unit);
else
  % ensure the units of the atlas are specified
  try
    atlas = ft_determine_units(atlas);
  catch
    % ft_determine_units will fail for triangle-only gifties.
  end
end
