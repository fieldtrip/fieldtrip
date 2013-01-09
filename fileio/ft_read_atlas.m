function atlas = ft_read_atlas(filename, varargin)

% FT_READ_ATLAS reads an template/individual segmentation or parcellation
% from disk. The volumetric segmentation or the surface-based parcellation
% can either represent a template atlas (eg. AAL or the Talairach Daemon),
% it can represent an individualized atlas (e.g. obtained from FreeSurfer)
% or it can represent an unlabeled parcellation obtained from the
% individual's DTi or resting state fMRI.
%
% Use as
%   atlas = ft_read_atlas(filename, ...)
%   atlas = ft_read_atlas({filenamelabels, filenamemesh}, ...)
%
% For individual surface based atlases two filenames are needed:
%   Filenamelabels points to the file that contains information with respect to
%     the parcels' labels.
%   Filenamemesh points to the file that defines the mesh on which the
%     parcellation is defined.
%
% The output atlas will be represented as structure according to
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

% deal with multiple filenames
if isa(filename, 'cell')
  if numel(filename)==2
    filenamemesh = filename{2};
    filename     = filename{1};
  else
    error('with multiple filenames, only 2 files are allowed');
  end % if precisely two input files
end % iscell

% optionally get the data from the URL and make a temporary local copy
filename = fetch_url(filename);

[p, f, x] = fileparts(filename);

if strcmp(f, 'TTatlas+tlrc')
  defaultformat = 'afni';
elseif strcmp(x, '.nii') && exist(fullfile(p, [f '.txt']))
  % this is a combination of nii+txt file, where the txt file contains three columns like this
  %   FAG	Precentral_L	2001
  %   FAD	Precentral_R	2002
  %   F1G	Frontal_Sup_L	2101
  %   F1D	Frontal_Sup_R	2102
  %   F1OG	Frontal_Sup_Orb_L	2111
  %   ...
  defaultformat  = 'aal';
else
  defaultformat  = 'wfu';
end

% get the optional input arguments
atlasformat = ft_getopt(varargin, 'format', defaultformat);

switch atlasformat
  case 'aal'
    labelfile = fullfile(p, [f '.txt']);
    fid = fopen(labelfile, 'rt');
    C = textscan(fid, '%s%s%d');
    lab = C{2};
    idx = C{3};
    fclose(fid);
    
    atlas = ft_read_mri(filename);
    atlas.tissue = atlas.anatomy;
    atlas = rmfield(atlas, 'anatomy');
    atlas.tissuelabel       = {};
    atlas.tissuelabel(idx)  = lab;
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
    
  case {'freesurfer_a2009s' 'freesurfer_aparc' 'freesurfer_ba'}
    % ensure freesurfer on the path and get the info how to get from value to label
    ft_hastoolbox('freesurfer', 1);
    
    if strcmp(atlasformat, 'freesurfer_a2009s')
      lookuptable = 'Simple_surface_labels2009.txt';
      parcelfield = 'a2009s';
      
      index = (0:75)';
      
      label = {'Unknown'
        'G_and_S_frontomargin'
        'G_and_S_occipital_inf'
        'G_and_S_paracentral'
        'G_and_S_subcentral'
        'G_and_S_transv_frontopol'
        'G_and_S_cingul-Ant'
        'G_and_S_cingul-Mid-Ant'
        'G_and_S_cingul-Mid-Post'
        'G_cingul-Post-dorsal'
        'G_cingul-Post-ventral'
        'G_cuneus'
        'G_front_inf-Opercular'
        'G_front_inf-Orbital'
        'G_front_inf-Triangul'
        'G_front_middle'
        'G_front_sup'
        'G_Ins_lg_and_S_cent_ins'
        'G_insular_short'
        'G_occipital_middle'
        'G_occipital_sup'
        'G_oc-temp_lat-fusifor'
        'G_oc-temp_med-Lingual'
        'G_oc-temp_med-Parahip'
        'G_orbital'
        'G_pariet_inf-Angular'
        'G_pariet_inf-Supramar'
        'G_parietal_sup'
        'G_postcentral'
        'G_precentral'
        'G_precuneus'
        'G_rectus'
        'G_subcallosal'
        'G_temp_sup-G_T_transv'
        'G_temp_sup-Lateral'
        'G_temp_sup-Plan_polar'
        'G_temp_sup-Plan_tempo'
        'G_temporal_inf'
        'G_temporal_middle'
        'Lat_Fis-ant-Horizont'
        'Lat_Fis-ant-Vertical'
        'Lat_Fis-post'
        'Medial_wall'
        'Pole_occipital'
        'Pole_temporal'
        'S_calcarine'
        'S_central'
        'S_cingul-Marginalis'
        'S_circular_insula_ant'
        'S_circular_insula_inf'
        'S_circular_insula_sup'
        'S_collat_transv_ant'
        'S_collat_transv_post'
        'S_front_inf'
        'S_front_middle'
        'S_front_sup'
        'S_interm_prim-Jensen'
        'S_intrapariet_and_P_trans'
        'S_oc_middle_and_Lunatus'
        'S_oc_sup_and_transversal'
        'S_occipital_ant'
        'S_oc-temp_lat'
        'S_oc-temp_med_and_Lingual'
        'S_orbital_lateral'
        'S_orbital_med-olfact'
        'S_orbital-H_Shaped'
        'S_parieto_occipital'
        'S_pericallosal'
        'S_postcentral'
        'S_precentral-inf-part'
        'S_precentral-sup-part'
        'S_suborbital'
        'S_subparietal'
        'S_temporal_inf'
        'S_temporal_sup'
        'S_temporal_transverse'};
      
      rgb   = [  0   0   0
        23 220  60
        23  60 180
        63 100  60
        63  20 220
        13   0 250
        26  60   0
        26  60  75
        26  60 150
        25  60 250
        60  25  25
        180  20  20
        220  20 100
        140  60  60
        180 220 140
        140 100 180
        180  20 140
        23  10  10
        225 140 140
        180  60 180
        20 220  60
        60  20 140
        220 180 140
        65 100  20
        220  60  20
        20  60 220
        100 100  60
        220 180 220
        20 180 140
        60 140 180
        25  20 140
        20  60 100
        60 220  20
        60  60 220
        220  60 220
        65 220  60
        25 140  20
        220 220 100
        180  60  60
        61  20 220
        61  20  60
        61  60 100
        25  25  25
        140  20  60
        220 180  20
        63 180 180
        221  20  10
        221  20 100
        221  60 140
        221  20 220
        61 220 220
        100 200 200
        10 200 200
        221 220  20
        141  20 100
        61 220 100
        141  60  20
        143  20 220
        101  60 220
        21  20 140
        61  20 180
        221 140  20
        141 100 220
        221 100  20
        181 200  20
        101  20  20
        101 100 180
        181 220  20
        21 140 200
        21  20 240
        21  20 200
        21  20  60
        101  60  60
        21 180 180
        223 220  60
        221  60  60];
      
    elseif strcmp(atlasformat, 'freesurfer_aparc')
      lookuptable = 'colortable_desikan_killiany.txt';
      parcelfield = 'aparc';
      
      index = (0:35)';
      
      label = {'unknown'
        'bankssts'
        'caudalanteriorcingulate'
        'caudalmiddlefrontal'
        'corpuscallosum'
        'cuneus'
        'entorhinal'
        'fusiform'
        'inferiorparietal'
        'inferiortemporal'
        'isthmuscingulate'
        'lateraloccipital'
        'lateralorbitofrontal'
        'lingual'
        'medialorbitofrontal'
        'middletemporal'
        'parahippocampal'
        'paracentral'
        'parsopercularis'
        'parsorbitalis'
        'parstriangularis'
        'pericalcarine'
        'postcentral'
        'posteriorcingulate'
        'precentral'
        'precuneus'
        'rostralanteriorcingulate'
        'rostralmiddlefrontal'
        'superiorfrontal'
        'superiorparietal'
        'superiortemporal'
        'supramarginal'
        'frontalpole'
        'temporalpole'
        'transversetemporal'
        'insula'};
      
      rgb   = [ 25   5  25
        25 100  40
        125 100 160
        100  25   0
        120  70  50
        220  20 100
        220  20  10
        180 220 140
        220  60 220
        180  40 120
        140  20 140
        20  30 140
        35  75  50
        225 140 140
        200  35  75
        160 100  50
        20 220  60
        60 220  60
        220 180 140
        20 100  50
        220  60  20
        120 100  60
        220  20  20
        220 180 220
        60  20 220
        160 140 180
        80  20 140
        75  50 125
        20 220 160
        20 180 140
        140 220 220
        80 160  20
        100   0 100
        70  20 170
        150 150 200
        255 192  32];
      
    elseif strcmp(atlasformat, 'freesurfer_ba')
      lookuptable = 'colortable_BA.txt';
      parcelfield = 'BA';
      
      index = (0:12)';
      
      label = {'unknown'
        'BA1'
        'BA2'
        'BA3a'
        'BA3b'
        'BA4a'
        'BA4p'
        'BA6'
        'BA44'
        'BA45'
        'V1'
        'V2'
        'MT'};
      
      rgb   = [25  5   25
        0   92  23
        131 148 255
        0   0   255
        255 102 51
        196 255 20
        255 51  204
        1   38  153
        153 0   38
        115 153 0
        153 15  0
        0   214 129
        155 0   153];
      
    else
      error('unknown freesurfer parcellation method requested');
    end
    
    %[index, label, rgb] = read_fscolorlut(lookuptable);
    %label = cellstr(label);
    rgb   = rgb(:,1) + rgb(:,2)*256 + rgb(:,3)*256*256;
    
    % read in the file
    switch ft_filetype(filename)
      case 'caret_label'
        p = gifti(filename);
        p = p.cdata;
      case 'freesurfer_annot'
        [v, p, c] = read_annotation(filename);
      otherwise
        error('unsupported fileformat for parcel file');
    end
    
    % read in the mesh
    switch ft_filetype(filenamemesh)
      case {'caret_surf' 'gifti'}
        tmp = gifti(filenamemesh);
        bnd.pnt = warp_apply(tmp.mat, tmp.vertices);
        bnd.tri = tmp.faces;
      case 'freesurfer_triangle_binary'
        [pnt, tri] = read_surf(filenamemesh);
        bnd.pnt    = pnt;
        bnd.tri    = tri;
      otherwise
        error('unsupported fileformat for surface mesh');
    end
    
    % check the number of vertices
    if size(bnd.pnt,1) ~= numel(p)
      error('the number of vertices in the mesh does not match the number of elements in the parcellation');
    end
    
    % reindex the parcels
    newp = zeros(size(p));
    for k = 1:numel(label)
      newp(p==rgb(k)) = index(k);
    end
    
    atlas       = [];
    atlas.pos   = bnd.pnt;
    atlas.tri   = bnd.tri;
    atlas.(parcelfield)            = newp;
    atlas.([parcelfield, 'label']) = label(2:end);
    atlas       = ft_convert_units(atlas);
    
  otherwise
    error('unsupported atlas format %s', atlasformat);
end % case
