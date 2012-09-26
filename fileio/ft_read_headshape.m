function [shape] = ft_read_headshape(filename, varargin)

% FT_READ_HEADSHAPE reads the fiducials and/or the measured headshape
% from a variety of files (like CTF and Polhemus). The headshape and
% fiducials can for example be used for coregistration.
%
% Use as
%   [shape] = ft_read_headshape(filename)
%
% Optional input arguments should be specified as key-value pairs and
% should include
%   format		= string, see below
%   coordinates = string, e.g. 'head' or 'dewar' (CTF)
%   unit        = string, e.g. 'cm'
%
% Supported input formats are
%   'ctf_*'
%   '4d_*'
%   'itab_asc'
%   'neuromag_*'
%   'mne_source'
%   'yokogawa_*'
%   'polhemus_*'
%   'spmeeg_mat'
%   'matlab'
%   'freesurfer_*'
%   'stl'          STereoLithography file format (often supported by
%                  CAD/generic 3D mesh editing programs)
%   'off'
%   'mne_*'        MNE surface description in ascii format ('mne_tri')
%                  or MNE source grid in ascii format, described as 3D
%                  points ('mne_pos')
%   'netmeg'
%   'vista'
%   'tet'
%   'tetgen'
%
% See also FT_READ_VOL, FT_READ_SENS, FT_WRITE_HEADSHAPE

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

% optionally get the data from the URL and make a temporary local copy
filename = fetch_url(filename);

% Check the input, if filename is a cell-array, call ft_read_headshape recursively and combine the outputs.
% This is used to read the left and right hemisphere of a Freesurfer cortical segmentation.
if iscell(filename)
  for i=1:numel(filename)
    tmpshape = ft_read_headshape(filename{i}, varargin{:});
    [path,name,ext] = fileparts(filename{i});
    if strcmp(ext, '.inflated') && strcmp(name, 'lh')
      % assume freesurfer inflated mesh in mm, mni space
      % move the mesh a bit to the left, to avoid overlap with the right
      % hemisphere
      tmpshape.pnt(:,1) = tmpshape.pnt(:,1) - max(tmpshape.pnt(:,1)) - 10;
      
    elseif strcmp(ext, '.inflated') && strcmp(name, 'rh')
      % id.
      % move the mesh a bit to the right, to avoid overlap with the left
      % hemisphere
      tmpshape.pnt(:,1) = tmpshape.pnt(:,1) - min(tmpshape.pnt(:,1)) + 10;
    end
    
    % try to read info regarding the sulcal pattern
    ft_hastoolbox('freesurfer', 1);
    try,
      tmpsulc = read_curv(fullfile(path, [name,'.sulc']));
    catch
      tmpsulc = [];
    end
    
    if i==1,
      shape = tmpshape;
      if ~isempty(tmpsulc)
        shape.sulc = tmpsulc;
      end
    else
      tmpshape  = ft_convert_units(tmpshape, shape.unit);
      npnt      = size(shape.pnt,1);
      shape.pnt = cat(1, shape.pnt, tmpshape.pnt);
      if isfield(shape, 'tri') && isfield(tmpshape, 'tri')
        shape.tri = cat(1, shape.tri, tmpshape.tri + npnt);
      elseif ~isfield(shape, 'tri') && ~isfield(tmpshape, 'tri')
        % this is ok
      else
        error('not all input files seem to contain a triangulation');
      end
      if ~isempty(tmpsulc)
        shape.sulc = cat(1, shape.sulc, tmpsulc);
      end
    end
  end
  return
end

% get the options
fileformat  = ft_getopt(varargin, 'format');
coordinates = ft_getopt(varargin, 'coordinates', 'head'); % the alternative for CTF coil positions is dewar
unit        = ft_getopt(varargin, 'unit'); % the default for yokogawa is cm, see below
annotationfile = ft_getopt(varargin, 'annotationfile');

if isempty(fileformat)
  % only do the autodetection if the format was not specified
  fileformat = ft_filetype(filename);
end

if ~isempty(annotationfile) && ~strcmp(fileformat, 'mne_source')
  error('at present extracting annotation information only works in conjunction with mne_source files');
end

% test whether the file exists
% FIXME why is this check not done for tetgen?
if ~strcmp(fileformat, 'tetgen') && ~exist(filename)
  error('file ''%s'' does not exist', filename);
end

% start with an empty structure
shape           = [];
shape.pnt       = [];
% FIXME is it required that it has an empty fiducial substructure? -> perhaps for SPM
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
    
    % fid = fid(1:3, :);
    
    [junk, NZ] = max(fid(1:3,1));
    [junk, L]  = max(fid(1:3,2));
    [junk, R]  = min(fid(1:3,2));
    rest       = setdiff(1:size(fid,1),[NZ L R]);
    
    shape.fid.pnt = fid([NZ L R rest], :);
    shape.fid.label = {'NZ', 'L', 'R'};
    if ~isempty(rest),
      for i = 4:size(fid,1)
        shape.fid.label{i} = ['fiducial' num2str(i)];
        % in a 5 coil configuration this corresponds with Cz and Inion
      end
    end
    
  case 'itab_asc'
    shape = read_itab_asc(filename);
    
  case 'gifti'
    ft_hastoolbox('gifti', 1);
    g = gifti(filename);
    if ~isfield(g, 'vertices')
      error('%s does not contain a tesselated surface', filename);
    end
    shape.pnt = warp_apply(g.mat, g.vertices);
    shape.tri = g.faces;
    if isfield(g, 'cdata')
      shape.mom = g.cdata;
    end
    
  case 'neuromag_mex'
    [co,ki,nu] = hpipoints(filename);
    fid = co(:,find(ki==1))';
    
    [junk, NZ] = max(fid(:,2));
    [junk, L]  = min(fid(:,1));
    [junk, R]  = max(fid(:,1));
    
    shape.fid.pnt = fid([NZ L R], :);
    shape.fid.label = {'NZ', 'L', 'R'};
    
  case 'mne_source'
    % read the source space from an MNE file
    ft_hastoolbox('mne', 1);
    
    src = mne_read_source_spaces(filename, 1);
    
    if ~isempty(annotationfile)
      ft_hastoolbox('freesurfer', 1);
      if numel(annotationfile)~=2
        error('two annotationfiles expected, one for each hemisphere');
      end
      for k = 1:numel(annotationfile)
        [v{k}, label{k}, c(k)] = read_annotation(annotationfile{k}, 1);
      end
      
      % match the annotations with the src structures
      if src(1).np == numel(label{1}) && src(2).np == numel(label{2})
        src(1).labelindx = label{1};
        src(2).labelindx = label{2};
      elseif src(1).np == numel(label{2}) && src(1).np == numel(label{1})
        src(1).labelindx = label{2};
        src(2).labelindx = label{1};
      else
        warning('incompatible annotation with triangulations, not using annotation information');
      end
      if ~isequal(c(1),c(2))
        error('the annotation tables differ, expecting equal tables for the hemispheres');
      end
      c = c(1);
    end

    shape = [];
    % only keep the points that are in use
    inuse1 = src(1).inuse==1;
    inuse2 = src(2).inuse==1;
    shape.pnt=[src(1).rr(inuse1,:); src(2).rr(inuse2,:)];
    
    % only keep the triangles that are in use; these have to be renumbered
    newtri1 = src(1).use_tris;
    newtri2 = src(2).use_tris;
    for i=1:numel(src(1).vertno)
      newtri1(newtri1==src(1).vertno(i)) = i;
    end
    for i=1:numel(src(2).vertno)
      newtri2(newtri2==src(2).vertno(i)) = i;
    end
    shape.tri  = [newtri1; newtri2 + numel(src(1).vertno)];
    shape.orig.pnt = [src(1).rr; src(2).rr];
    shape.orig.tri = [src(1).tris; src(2).tris + src(1).np];
    shape.orig.inuse = [src(1).inuse src(2).inuse]';
    if isfield(src(1), 'labelindx')
      shape.orig.labelindx = [src(1).labelindx;src(2).labelindx];
      shape.labelindx      = [src(1).labelindx(inuse1); src(2).labelindx(inuse2)];
%      ulabelindx = unique(c.table(:,5));
%       for k = 1:c.numEntries
%         % the values are really high (apart from the 0), so I guess it's safe to start
%         % numbering from 1
%         shape.orig.labelindx(shape.orig.labelindx==ulabelindx(k)) = k;
%         shape.labelindx(shape.labelindx==ulabelindx(k)) = k;
%       end
% FIXME the above screws up the interpretation of the labels, because the
% color table is not sorted
      shape.label = c.struct_names;
      shape.annotation = c.orig_tab; % to be able to recover which one
      shape.ctable = c.table;
    end
    
  case {'neuromag_mne', 'neuromag_fif'}
    % read the headshape and fiducials from an MNE file
    hdr = ft_read_header(filename,'headerformat','neuromag_mne');
    nFid = size(hdr.orig.dig,2); %work out number of fiducials
    switch coordinates
      case 'head' % digitiser points should be stored in head coordinates by default
        
        fidN=1;
        pntN=1;
        for i=1:nFid % loop over fiducials
          % check that this point is in head coordinates
          % 0 is unknown
          % 4 is fiducial system, i.e. head coordinates
          if hdr.orig.dig(i).coord_frame~=4
            warning(['Digitiser point (' num2str(i) ') not stored in head coordinates!']);
          end
          
          switch hdr.orig.dig(i).kind % constants defined in MNE - see p.215 of MNE manual
            case 1 % Cardinal point (Nasion, LPA or RPA)
              % get location of fiducial:
              shape.fid.pnt(fidN,1:3) = hdr.orig.dig(i).r*100; % multiply by 100 to convert to cm
              switch hdr.orig.dig(i).ident
                case 1 % LPA
                  shape.fid.label{fidN} = 'LPA';
                case 2 % nasion
                  shape.fid.label{fidN} = 'Nasion';
                case 3 % RPA
                  shape.fid.label{fidN} = 'RPA';
                otherwise
                  error('Unidentified cardinal point in file');
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
    
  case {'yokogawa_mrk', 'yokogawa_ave', 'yokogawa_con', 'yokogawa_raw' }
    if ft_hastoolbox('yokogawa_meg_reader')
      hdr = read_yokogawa_header_new(filename);
      marker = hdr.orig.coregist.hpi;
    else
      hdr = read_yokogawa_header(filename);
      marker = hdr.orig.matching_info.marker;
    end
    
    % markers 1-3 identical to zero: try *.mrk file
    if ~any([marker(:).meg_pos])
      [p, f, x] = fileparts(filename);
      filename = fullfile(p, [f '.mrk']);
      if exist(filename, 'file')
        if ft_hastoolbox('yokogawa_meg_reader')
          hdr = read_yokogawa_header_new(filename);
          marker = hdr.orig.coregist.hpi;
        else
          hdr = read_yokogawa_header(filename);
          marker = hdr.orig.matching_info.marker;
        end
      end
    end
    
    % non zero markers 1-3
    if any([marker(:).meg_pos])
      shape.fid.pnt = cat(1, marker(1:5).meg_pos);
      sw_ind = [3 1 2];
      shape.fid.pnt(1:3,:)= shape.fid.pnt(sw_ind, :);
      shape.fid.label = {'nas'; 'lpa'; 'rpa'; 'Marker4'; 'Marker5'};
    else
      error('no coil information found in Yokogawa file');
    end
    
    % convert to the units of the grad, the desired default for yokogawa is centimeter.
    shape = ft_convert_units(shape, 'cm');
    
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
    
  case 'yokogawa_hsp'
    fid = fopen(filename, 'rt');
    
    fidstart = false;
    hspstart = false;
    
    % try to locate the fiducial positions
    while ~fidstart && ~feof(fid)
      line = fgetl(fid);
      if ~isempty(strmatch('//Position of fiducials', line))
        fidstart = true;
      end
    end
    if fidstart
      line_xpos = fgetl(fid);
      line_ypos = fgetl(fid);
      line_yneg = fgetl(fid);
      xpos = sscanf(line_xpos(3:end), '%f');
      ypos = sscanf(line_ypos(3:end), '%f');
      yneg = sscanf(line_yneg(3:end), '%f');
      shape.fid.pnt = [
        xpos(:)'
        ypos(:)'
        yneg(:)'
        ];
      shape.fid.label = {
        'X+'
        'Y+'
        'Y-'
        };
    end
    
    % try to locate the fiducial positions
    while ~hspstart && ~feof(fid)
      line = fgetl(fid);
      if ~isempty(strmatch('//No of rows', line))
        hspstart = true;
      end
    end
    if hspstart
      line = fgetl(fid);
      siz = sscanf(line, '%f');
      shape.pnt = zeros(siz(:)');
      for i=1:siz(1)
        line = fgetl(fid);
        shape.pnt(i,:) = sscanf(line, '%f');
      end
    end
    
    fclose(fid);
    
  case 'polhemus_fil'
    [shape.fid.pnt, shape.pnt, shape.fid.label] = read_polhemus_fil(filename, 0);
    
  case 'polhemus_pos'
    [shape.fid.pnt, shape.pnt, shape.fid.label] = read_ctf_pos(filename);
    
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
    elseif isfield(tmp, 'bnd')
      % the variable in the file is most likely a precomputed triangulation of some
      % sort
      shape = tmp.bnd;
    elseif isfield(tmp, 'elec')
      tmp.elec        = ft_datatype_sens(tmp.elec);
      shape.fid.pnt   = tmp.elec.chanpos;
      shape.fid.label = tmp.elec.label;
    else
      error('no headshape found in Matlab file');
    end
    
  case {'freesurfer_triangle_binary', 'freesurfer_quadrangle'}
    % the freesurfer toolbox is required for this
    ft_hastoolbox('freesurfer', 1);
    [pnt, tri] = read_surf(filename);
    if min(tri(:)) == 0
      % start counting from 1
      tri = tri + 1;
    end
    shape.pnt = pnt;
    shape.tri = tri;
    shape = rmfield(shape, 'fid');
    
  case 'stl'
    [pnt, tri, nrm] = read_stl(filename);
    shape.pnt = pnt;
    shape.tri = tri;
    
  case 'off'
    [pnt, plc] = read_off(filename);
    shape.pnt  = pnt;
    shape.tri  = plc;
    
  case 'mne_tri'
    % FIXME this should be implemented, consistent with ft_write_headshape
    keyboard
    
  case 'mne_pos'
    % FIXME this should be implemented, consistent with ft_write_headshape
    keyboard
    
  case 'netmeg'
    hdr = ft_read_header(filename);
    if isfield(hdr.orig, 'headshapedata')
      shape.pnt = hdr.orig.Var.headshapedata;
    else
      error('the NetMEG file "%s" does not contain headshape data', filename);
    end
    
  case 'vista'
    ft_hastoolbox('simbio', 1);
    [nodes,elements,labels] = read_vista_mesh(filename);
    shape.pnt     = nodes;
    if size(elements,2)==8
      shape.hex     = elements;
    elseif size(elements,2)==4
      shape.tet = elements;
    else
      error('unknown elements format')
    end
    shape.index = labels;
    
  case 'tet'
    % the toolbox from Gabriel Peyre has a function for this
    ft_hastoolbox('toolbox_graph', 1);
    [vertex, face] = read_tet(filename);
    %     'vertex' is a '3 x nb.vert' array specifying the position of the vertices.
    %     'face' is a '4 x nb.face' array specifying the connectivity of the tet mesh.
    shape.pnt = vertex';
    shape.tet = face';
    
  case 'tetgen'
    % reads in the tetgen format and rearranges according to FT conventions
    % tetgen files also return a 'faces' field (not used here)
    shape = rmfield(shape,'fid');
    IMPORT = importdata([filename '.1.ele'],' ',1);
    shape.tet = IMPORT.data(:,2:5);
    IMPORT = importdata([filename '.1.node'],' ',1);
    shape.pnt = IMPORT.data(:,2:4);
    
  otherwise
    % try reading it from an electrode of volume conduction model file
    success = false;
    
    if ~success
      % try reading it as electrode positions
      % and treat those as fiducials
      try
        elec = ft_read_sens(filename);
        if ~ft_senstype(elec, 'eeg')
          error('headshape information can not be read from MEG gradiometer file');
        else
          shape.fid.pnt   = elec.chanpos;
          shape.fid.label = elec.label;
          success = 1;
        end
      catch
        success = false;
      end % try
    end
    
    if ~success
      % try reading it as volume conductor
      % and treat the skin surface as headshape
      try
        vol = ft_read_vol(filename);
        if ~ft_voltype(vol, 'bem')
          error('skin surface can only be extracted from boundary element model');
        else
          if ~isfield(vol, 'skin')
            vol.skin = find_outermost_boundary(vol.bnd);
          end
          shape.pnt = vol.bnd(vol.skin).pnt;
          shape.tri = vol.bnd(vol.skin).tri; % also return the triangulation
          success = 1;
        end
      catch
        success = false;
      end % try
    end
    
    if ~success
      error('unknown fileformat "%s" for head shape information', fileformat);
    end
end

if isfield(shape, 'fid') && isfield(shape.fid, 'label')
  % ensure that it is a column
  shape.fid.label = shape.fid.label(:);
end

% this will add the units to the head shape and optionally convert
if ~isempty(unit)
  shape = ft_convert_units(shape, unit);
else
  shape = ft_convert_units(shape);
end
