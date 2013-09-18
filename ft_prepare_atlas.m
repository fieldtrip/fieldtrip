function [atlas, cfg] = ft_prepare_atlas(cfg)

% FT_PREPARE_ATLAS reads in a specified atlas with coordinates and
% anatomical labels. It either uses the AFNI brik file that is available
% from http://afni.nimh.nih.gov/afni/doc/misc/ttatlas_tlrc, or it
% uses one of the WFU atlasses available from http://fmri.wfubmc.edu.
%
% This function is called by other FieldTrip functions that make
% use of an atlas, for example for plotting or for selection of an
% anatomical region of interest.
%
% Use as
%   [atlas] = ft_prepare_atlas(cfg)
%
% where the configuration should contain
%   cfg.atlas      string, filename of the atlas to use
%
% See also FT_VOLUMELOOKUP, FT_SOURCEPLOT, FT_SOURCESTATISTICS

% Copyright (C) 2005-2011, Robert Oostenveld, Ingrid Nieuwenhuis
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

revision = '$Id$';

if ischar(cfg)
  % prior to 7 December 2011, this function was called with the filename as input
  % for consistency with other fieldtrip functions, it is now in a cfg field
  cfg.atlas = cfg;
end

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble provenance

useafni = 0;
usewfu  = 0;

[p, f, x] = fileparts(cfg.atlas);

if strcmp(f, 'TTatlas+tlrc')
  useafni = 1;
else
  usewfu = 1;
end

if useafni
  % check whether the required AFNI toolbox is available
  ft_hastoolbox('afni', 1);
  
  atlas = ft_read_mri(cfg.atlas);
  
  % the AFNI atlas contains two volumes at 1mm resolution
  atlas.brick0 = atlas.anatomy(:,:,:,1);
  atlas.brick1 = atlas.anatomy(:,:,:,2);
  atlas = rmfield(atlas, 'anatomy');
  atlas.coord = 'tal';
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % the following information is from https://afni.nimh.nih.gov/afni/doc/misc/ttatlas_tlrc
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % column 1 contains sub-brick
  % column 2 contains value
  % column 3 contains stucture name
  % 1   68  Hippocampus
  % 1   71  Amygdala
  % 0   20  Posterior Cingulate
  % 0   21  Anterior Cingulate
  % 0   22  Subcallosal Gyrus
  % 0   24  Transverse Temporal Gyrus
  % 0   25  Uncus
  % 0   26  Rectal Gyrus
  % 0   27  Fusiform Gyrus
  % 0   28  Inferior Occipital Gyrus
  % 0   29  Inferior Temporal Gyrus
  % 0   30  Insula
  % 0   31  Parahippocampal Gyrus
  % 0   32  Lingual Gyrus
  % 0   33  Middle Occipital Gyrus
  % 0   34  Orbital Gyrus
  % 0   35  Middle Temporal Gyrus
  % 0   36  Superior Temporal Gyrus
  % 0   37  Superior Occipital Gyrus
  % 0   39  Inferior Frontal Gyrus
  % 0   40  Cuneus
  % 0   41  Angular Gyrus
  % 0   42  Supramarginal Gyrus
  % 0   43  Cingulate Gyrus
  % 0   44  Inferior Parietal Lobule
  % 0   45  Precuneus
  % 0   46  Superior Parietal Lobule
  % 0   47  Middle Frontal Gyrus
  % 0   48  Paracentral Lobule
  % 0   49  Postcentral Gyrus
  % 0   50  Precentral Gyrus
  % 0   51  Superior Frontal Gyrus
  % 0   52  Medial Frontal Gyrus
  % 0   70  Lentiform Nucleus
  % 1   72  Hypothalamus
  % 1   73  Red Nucleus
  % 1   74  Substantia Nigra
  % 0   75  Claustrum
  % 0   76  Thalamus
  % 0   77  Caudate
  % 1   124 Caudate Tail
  % 1   125 Caudate Body
  % 1   126 Caudate Head
  % 1   128 Ventral Anterior Nucleus
  % 1   129 Ventral Posterior Medial Nucleus
  % 1   130 Ventral Posterior Lateral Nucleus
  % 1   131 Medial Dorsal Nucleus
  % 1   132 Lateral Dorsal Nucleus
  % 1   133 Pulvinar
  % 1   134 Lateral Posterior Nucleus
  % 1   135 Ventral Lateral Nucleus
  % 1   136 Midline Nucleus
  % 1   137 Anterior Nucleus
  % 1   138 Mammillary Body
  % 1   144 Medial Globus Pallidus
  % 1   145 Lateral Globus Pallidus
  % 1   151 Putamen
  % 1   146 Nucleus Accumbens
  % 1   147 Medial Geniculum Body
  % 1   148 Lateral Geniculum Body
  % 1   149 Subthalamic Nucleus
  % 1   81  Brodmann area 1
  % 1   82  Brodmann area 2
  % 1   83  Brodmann area 3
  % 1   84  Brodmann area 4
  % 1   85  Brodmann area 5
  % 1   86  Brodmann area 6
  % 1   87  Brodmann area 7
  % 1   88  Brodmann area 8
  % 1   89  Brodmann area 9
  % 1   90  Brodmann area 10
  % 1   91  Brodmann area 11
  % 1   93  Brodmann area 13
  % 1   94  Brodmann area 17
  % 1   95  Brodmann area 18
  % 1   96  Brodmann area 19
  % 1   97  Brodmann area 20
  % 1   98  Brodmann area 21
  % 1   99  Brodmann area 22
  % 1   100 Brodmann area 23
  % 1   101 Brodmann area 24
  % 1   102 Brodmann area 25
  % 1   103 Brodmann area 27
  % 1   104 Brodmann area 28
  % 1   105 Brodmann area 29
  % 1   106 Brodmann area 30
  % 1   107 Brodmann area 31
  % 1   108 Brodmann area 32
  % 1   109 Brodmann area 33
  % 1   110 Brodmann area 34
  % 1   111 Brodmann area 35
  % 1   112 Brodmann area 36
  % 1   113 Brodmann area 37
  % 1   114 Brodmann area 38
  % 1   115 Brodmann area 39
  % 1   116 Brodmann area 40
  % 1   117 Brodmann area 41
  % 1   118 Brodmann area 42
  % 1   119 Brodmann area 43
  % 1   120 Brodmann area 44
  % 1   121 Brodmann area 45
  % 1   122 Brodmann area 46
  % 1   123 Brodmann area 47
  % 0   53  Uvula of Vermis
  % 0   54  Pyramis of Vermis
  % 0   55  Tuber of Vermis
  % 0   56  Declive of Vermis
  % 0   57  Culmen of Vermis
  % 0   58  Cerebellar Tonsil
  % 0   59  Inferior Semi-Lunar Lobule
  % 0   60  Fastigium
  % 0   61  Nodule
  % 0   62  Uvula
  % 0   63  Pyramis
  % 0   66  Culmen
  % 0   65  Declive
  % 1   127 Dentate
  % 0   64  Tuber
  % 0   67  Cerebellar Lingual
  
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
  
elseif usewfu
  atlas = ft_read_mri(cfg.atlas);
  atlas.brick0 = atlas.anatomy(:,:,:);
  atlas = rmfield(atlas, 'anatomy');
  atlas.coord = 'mni';
  
  % the WFU atlas contains a single atlas volume at 2mm resolution
  % to keep it compatible with the existing code, add a dummy atlas volume
  atlas.brick1 = zeros(size(atlas.brick0));
  
  [p, f, x] = fileparts(cfg.atlas);
  filename1 = fullfile(p, [f '_List.mat']);
  filename2 = fullfile(p, [f '.txt']);
  
  if exist(filename1, 'file')
    % this is a mat file that Ingrid apparently discovered somewhere
    load(filename1);
    atlas.descr = [];
    atlas.descr.brick = zeros(length(ROI),1);
    atlas.descr.value = [ROI.ID]';
    atlas.descr.name  = {ROI.Nom_C}'; % what is difference between Nom_C and Nom_L??
    
  elseif exist(filename2, 'file')
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
    
    atlas.descr = [];
    atlas.descr.brick = [];
    atlas.descr.value = [];
    atlas.descr.name  = {};
    
    fid = fopen(filename2);
    i = 1;
    iswfu = false;
    while 1
      tline = fgetl(fid);
      if ~ischar(tline), break, end
      
      if i==1 && strcmp(tline(1), '[')
          iswfu = true;
      end
      % split into separate strings
      C = textscan(tline,'%s');
      if iswfu
        % the wfu version of the aal atlas the following format
        % num \t label
        num = C{1}{1}; 
        str = C{1}{2};
      else
        % the original aal atlas has the following format 
        % id \t label \t num
        num = C{1}{3}; 
        str = C{1}{2};
      end % if iswfu
      
      num = str2double(num);
      if ~isnan(num)
        atlas.descr.brick(i) = 0;
        atlas.descr.value(i) = num;
        atlas.descr.name{i}  = str;
        i = i+1;
      end % if num
    end
    fclose(fid);
    
  else
    error('cannot locate the file that maps the labels to the numbers');
  end
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble provenance
ft_postamble history atlas
