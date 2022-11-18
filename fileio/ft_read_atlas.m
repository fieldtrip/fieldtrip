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
%   'map'         = string, 'maxprob' (default), or 'prob', for FSL-based atlases, providing 
%                   either a probabilistic segmentation or a maximum a posterior probability map
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

% Copyright (C) 2005-2019, Robert Oostenveld, Ingrid Nieuwenhuis, Jan-Mathijs Schoffelen, Arjen Stolk
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

if contains(filename, 'BRIK') || contains(filename, 'HEAD')
  % this is robust for both compressed or uncompressed afni atlases. 
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
elseif strcmp(x, '.mgz') && contains(f, 'aparc') || contains(f, 'aseg')
  % individual volume based segmentation from freesurfer
  format = 'freesurfer_volume';
elseif ft_filetype(filename, 'caret_label')
  % this is a gifti file that contains both the values for a set of
  % vertices as well as the labels.
  format = 'caret_label';
elseif contains(filename, 'MPM')
  % assume to be from the spm_anatomy toolbox
  format = 'spm_anatomy';
elseif strcmp(x, '.xml') && (isfolder(strtok(fullfile(p,f), '_')) || isfolder(strtok(fullfile(p,f), '-')))
  % fsl-format atlas, this is assumed to consist of an .xml file that
  % specifies the labels, as well as the filenames of the files with the actual data stored
  % in a directory with the of the strtok'ed (with '-' or '_') file name.
  format = 'fsl';
elseif strcmp(x, '.mat')
  format = 'mat';
elseif strcmp(x, '.nii') && contains(f, 'Yeo2011_7Networks')
  format = 'yeo7';
elseif strcmp(x, '.nii') && contains(f, 'Yeo2011_17Networks')
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
    if compact && ~isempty(atlas.tissuelabel)
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
    if compact && ~isempty(atlas.tissuelabel)
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
    
    tmp     = ft_read_mri(filename);
    if isfield(tmp, 'coordsys') && ~strcmp(tmp.coordsys, 'unknown')
      coordsys = tmp.coordsys;
    elseif isfield(tmp.hdr, 'TEMPLATE_SPACE') && ~isempty(tmp.hdr.TEMPLATE_SPACE)
      coordsys = lower(tmp.hdr.TEMPLATE_SPACE); % FIXME this is based on AFNI conventions, not easily decodable by FT
    else
      coordsys = 'tal'; % FIXME could be different in other atlases
    end

    if isfield(tmp.hdr, 'ATLAS_LABEL_TABLE') && ~isempty(tmp.hdr.ATLAS_LABEL_TABLE)
      if isfield(tmp.hdr.ATLAS_LABEL_TABLE(1), 'sb_label') && ~all(tmp.anatomy(:)==round(tmp.anatomy(:)))
        % probabilistic atlas
        isprobabilistic = true;
      else
        % indexed atlas
        isprobabilistic = false;
      end
      labels  = {tmp.hdr.ATLAS_LABEL_TABLE.struct}';
      values  = [tmp.hdr.ATLAS_LABEL_TABLE.val]';
      
    elseif contains(filename, 'TTatlas+tlrc')
      isprobabilistic = false;

      % the following information is from https://sscc.nimh.nih.gov/afni/doc/misc/afni_ttatlas/index_html
      values = [
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
      
      labels = {
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

    else
      ft_error('no information about the atlas labels is available');
    end
    
    atlas     = [];
    atlas.dim = tmp.dim(1:3);
    atlas.transform = tmp.transform;
    atlas.hdr = tmp.hdr;
    atlas.coordsys = coordsys;
    
    nbrick  = size(tmp.anatomy,4);
    for k = 1:nbrick
      
      if ~isprobabilistic
        brickname = sprintf('brick%d',k-1);
        brick     = tmp.anatomy(:,:,:,k);
        ulabel    = setdiff(unique(brick(:)), 0);
        label     = cell(size(ulabel));
        nlabel    = numel(label);
        
        % renumber the brick from 1:N and keep track of the label
        newbrick  = zeros(size(brick));
        for i = 1:nlabel
          sel = find(values==ulabel(i));
          if ~isempty(sel)
            label(i) = labels(sel);
            newbrick(brick==ulabel(i)) = i;
          else
            ft_warning('the value %d does not have a label according to the ATLAS_LABEL_TABLE and will be discarded', ulabel(i));
          end
        end
        atlas.(brickname) = newbrick;
        atlas.([brickname 'label']) = label;
      else
        atlas.(labels{k}) = tmp.anatomy(:,:,:,values(k)+1); % indexing is 0-based in this case
      end
        
    end

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
      ft_warning('cannot locate %s, making default tissue labels', filename2);
      
      value = [];
      label = {};
      for i=1:max(brick0(:))
        % this is consistent with FIXSEGMENTATION
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
    fscolorlut = ft_getopt(varargin, 'fscolorlut', []);
    if isempty(fscolorlut)
      [ftver, ftpath] = ft_version;
      fscolorlut  = fullfile(ftpath, 'external/freesurfer', 'FreeSurferColorLUT.txt');
    end
    [value, label, rgbv] = read_fscolorlut(fscolorlut);
    label = cellstr(label);
    
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
    rgbv  = rgbv(sel,:);

    % remap the values in the data
    aparc = zeros(size(dat));
    rgba  = zeros(0,4);
    cnt   = 0;
    for k = 1:numel(value)
      sel = dat==value(k);
      if sum(sel(:))
        cnt = cnt+1;
        fprintf('re-indexing label %s to a value of %d (was %d)\n', label{k}, cnt, value(k));
        aparc(sel)        = cnt;
        aparclabel{cnt,1} = label{k};
        rgba(cnt,:)       = rgbv(k,:);
      end
    end
    atlas.aparc      = aparc;
    atlas.aparclabel = aparclabel;
    atlas.rgba       = rgba;
    
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
    if isfield(g, 'labels')
      label = g.labels.name(:);
      key   = g.labels.key(:);
      if isfield(g.labels, 'rgba')
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
          if ~isa(atlas.(fn), 'double') && ~isa(atlas.(fn), 'single')
            % ensure that the probabilistic values are either double or
            % single precision, do single precision to save memory
            atlas.(fn) = single(atlas.(fn));
          end
          if any(atlas.(fn)(:)>1) & all(atlas.(fn)(:)<=100)
            % convert to probability values, assuming 100 to be max
            atlas.(fn) = atlas.(fn)./100;
          end
        end
        atlas.coordsys = 'mni';
        atlas.dim      = atlas.dim(1:3);
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
