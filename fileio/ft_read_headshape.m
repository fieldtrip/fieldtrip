function [shape] = ft_read_headshape(filename, varargin)

% FT_READ_HEADSHAPE reads the fiducials and/or the measured headshape from a variety
% of files (like CTF and Polhemus). The headshape and fiducials can for example be
% used for coregistration.
%
% Use as
%   [shape] = ft_read_headshape(filename, ...)
% or
%   [shape] = ft_read_headshape({filename1, filename2}, ...)
%
% If you specify the filename as a cell-array, the following situations are supported:
%  - a two-element cell-array with the file names for the left and
%    right hemisphere, e.g. FreeSurfer's {'lh.orig' 'rh.orig'}, or
%    Caret's {'X.L.Y.Z.surf.gii' 'X.R.Y.Z.surf.gii'}
%  - a two-element cell-array points to files that represent
%    the coordinates and topology in separate files, e.g.
%    Caret's {'X.L.Y.Z.coord.gii' 'A.L.B.C.topo.gii'};
% By default all information from the two files will be concatenated (i.e. assumed to
% be the shape of left and right hemispeheres). The option 'concatenate' can be set
% to 'no' to prevent them from being concatenated in a single structure.
%
% Additional options should be specified in key-value pairs and can include
%   'format'      = string, see below
%   'coordsys'    = string, e.g. 'head' or 'dewar' (only supported for CTF)
%   'unit'        = string, e.g. 'mm' (default is the native units of the file)
%   'concatenate' = 'no' or 'yes' (default = 'yes')
%   'image'       = path to .jpeg file
%
% Supported input file formats include
%   'matlab'       containing FieldTrip or BrainStorm headshapes or cortical meshes
%   'stl'          STereoLithography file format, for use with CAD and/or generic 3D mesh editing programs
%   'vtk'          Visualization ToolKit file format, for use with Paraview
%   'mne_*'        MNE surface description in ASCII format ('mne_tri') or MNE source grid in ascii format, described as 3D points ('mne_pos')
%   'ctf_*'
%   '4d_*'
%   'itab_asc'
%   'neuromag_*'
%   'mne_source'
%   'yokogawa_*'
%   'polhemus_*'
%   'spmeeg_mat'
%   'off'
%   'freesurfer_*'
%   'netmeg'
%   'vista'
%   'tet'
%   'tetgen_ele'
%   'gifti'
%   'caret_surf'
%   'caret_coord'
%   'caret_topo'
%   'caret_spec'
%   'brainvisa_mesh'
%   'brainsuite_dfs'
%   'obj'           Wavefront .obj file obtained with the structure.io
%
% See also FT_READ_VOL, FT_READ_SENS, FT_READ_ATLAS, FT_WRITE_HEADSHAPE

% Copyright (C) 2008-2016 Robert Oostenveld
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

% get the options
annotationfile = ft_getopt(varargin, 'annotationfile');
concatenate    = ft_getopt(varargin, 'concatenate', 'yes');
coordsys       = ft_getopt(varargin, 'coordsys', 'head');    % for ctf or neuromag_mne coil positions, the alternative is dewar
fileformat     = ft_getopt(varargin, 'format');
unit           = ft_getopt(varargin, 'unit');
image          = ft_getopt(varargin, 'image',[100, 100 ,100]);               % path to .jpeg file

% Check the input, if filename is a cell-array, call ft_read_headshape recursively and combine the outputs.
% This is used to read the left and right hemisphere of a Freesurfer cortical segmentation.
if iscell(filename)
  
  for i = 1:numel(filename)
    tmp       = ft_read_headshape(filename{i}, varargin{:});
    haspos(i) = isfield(tmp, 'pos') && ~isempty(tmp.pos);
    hastri(i) = isfield(tmp, 'tri') && ~isempty(tmp.tri);
    if ~haspos(i), tmp.pos = []; end
    if ~hastri(i), tmp.tri = []; end
    if ~isfield(tmp, 'unit'), tmp.unit = 'unknown'; end
    bnd(i) = tmp;
  end
  
  % Concatenate the bnds (only if 'concatenate' = 'yes' ) and if all
  % structures have non-empty pnts and tris. If not, the input filenames
  % may have been caret-style coord and topo, which needs combination of
  % the pos and tri.
  
  if  numel(filename)>1 && all(haspos==1) && strcmp(concatenate, 'yes')
    if length(bnd)>2
      error('Cannot concatenate more than two files') % no more than two files are taken for cancatenation
    else
      fprintf('Concatenating the meshes in %s and %s\n', filename{1}, filename{2});
      
      shape     = [];
      shape.pos = cat(1, bnd.pos);
      npos      = size(bnd(1).pos,1);
      
      if isfield(bnd(1), 'tri')  && isfield(bnd(2), 'tri')
        shape.tri = cat(1, bnd(1).tri, bnd(2).tri + npos);
      elseif ~isfield(bnd(1), 'tri') && ~isfield(bnd(2), 'tri')
        % this is ok
      else
        error('not all input files seem to contain a triangulation');
      end
      
      % concatenate any other fields
      fnames = {'sulc' 'curv' 'area' 'thickness' 'atlasroi'};
      for k = 1:numel(fnames)
        if isfield(bnd(1), fnames{k}) && isfield(bnd(2), fnames{k})
          shape.(fnames{k}) = cat(1, bnd.(fnames{k}));
        elseif ~isfield(bnd(1), fnames{k}) && ~isfield(bnd(2), fnames{k})
          % this is ok
        else
          error('not all input files seem to contain a "%s"', fnames{k});
        end
      end
      
      shape.brainstructure = []; % keeps track of the order of files in concatenation
      for h = 1:length(bnd)
        shape.brainstructure  = [shape.brainstructure; h*ones(length(bnd(h).pos), 1)];
        [p,f,e]               = fileparts(filename{h});
        
        % do an educated guess, otherwise default to the filename
        iscortexright = ~isempty(strfind(f,'rh'));
        iscortexright = iscortexright || ~isempty(strfind(f,'.R.'));
        iscortexright = iscortexright || ~isempty(strfind(f,'Right'));
        iscortexright = iscortexright || ~isempty(strfind(f,'RIGHT'));
        
        iscortexleft = ~isempty(strfind(f,'lh'));
        iscortexleft = iscortexleft || ~isempty(strfind(f,'.L.'));
        iscortexleft = iscortexleft || ~isempty(strfind(f,'Left'));
        iscortexleft = iscortexleft || ~isempty(strfind(f,'LEFT'));
        
        if iscortexright && iscortexleft
          % something strange is going on, default to the filename and let the user take care of this
          shape.brainstructurelabel{h,1} = f;
        elseif iscortexleft
          shape.brainstructurelabel{h,1} = 'CORTEX_LEFT';
        elseif iscortexright
          shape.brainstructurelabel{h,1} = 'CORTEX_RIGHT';
        else
          % nothing to be guessed
          shape.brainstructurelabel{h,1} = f;
        end
      end
      
    end
  elseif numel(filename)>1 && ~all(haspos==1)
    if numel(bnd)>2
      error('Cannot combine more than two files') % no more than two files are taken for cancatenation
    else
      shape = [];
      if sum(haspos==1)==1
        fprintf('Using the vertex positions from %s\n', filename{find(haspos==1)});
        shape.pos  = bnd(haspos==1).pos;
        shape.unit = bnd(haspos==1).unit;
      else
        error('Don''t know what to do');
      end
      if sum(hastri==1)==1
        fprintf('Using the faces definition from %s\n', filename{find(hastri==1)});
        shape.tri = bnd(hastri==1).tri;
      end
      if max(shape.tri(:))~=size(shape.pos,1)
        error('mismatch in number of points in pos and tri');
      end
    end
    
  else
    % in case numel(filename)==1, or strcmp(concatenate, 'no')
    shape = bnd;
  end
  
  return
end % if iscell

% checks if there exists a .jpg file of 'filename'
  [pathstr,name]  = fileparts(filename);  
  if exist(fullfile(pathstr,[name,'.jpg']))
    image    = fullfile(pathstr,[name,'.jpg']);
    hasimage = 1;
  end  
    

% optionally get the data from the URL and make a temporary local copy
filename = fetch_url(filename);

if isempty(fileformat)
  % only do the autodetection if the format was not specified
  fileformat = ft_filetype(filename);
end

if ~isempty(annotationfile) && ~strcmp(fileformat, 'mne_source')
  error('at present extracting annotation information only works in conjunction with mne_source files');
end

% start with an empty structure
shape           = [];
shape.pos       = [];

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
    switch coordsys
      case 'head'
        shape.fid.pos = cell2mat(struct2cell(orig.head));
        shape.coordsys = 'ctf';
      case 'dewar'
        shape.fid.pos = cell2mat(struct2cell(orig.dewar));
        shape.coordsys = 'dewar';
      otherwise
        error('incorrect coordsys specified');
    end
    shape.fid.label = fieldnames(orig.head);
    
  case 'ctf_shape'
    orig = read_ctf_shape(filename);
    shape.pos = orig.pos;
    shape.fid.label = {'NASION', 'LEFT_EAR', 'RIGHT_EAR'};
    shape.fid.pos = zeros(0,3); % start with an empty array
    for i = 1:numel(shape.fid.label)
      shape.fid.pos = cat(1, shape.fid.pos, getfield(orig.MRI_Info, shape.fid.label{i}));
    end
    
  case {'4d_xyz', '4d_m4d', '4d_hs', '4d', '4d_pdf'}
    [p, f, x] = fileparts(filename);
    if ~strcmp(fileformat, '4d_hs')
      filename = fullfile(p, 'hs_file');
    end
    [shape.pos, fid] = read_bti_hs(filename);
    
    % I'm making some assumptions here
    % which I'm not sure will work on all 4D systems
    
    % fid = fid(1:3, :);
    
    [junk, NZ] = max(fid(1:3,1));
    [junk, L]  = max(fid(1:3,2));
    [junk, R]  = min(fid(1:3,2));
    rest       = setdiff(1:size(fid,1),[NZ L R]);
    
    shape.fid.pos = fid([NZ L R rest], :);
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
    shape.pos = ft_warp_apply(g.mat, g.vertices);
    shape.tri = g.faces;
    shape.unit = 'mm';  % defined in the GIFTI standard to be milimeter
    if isfield(g, 'cdata')
      shape.mom = g.cdata;
    end
    
  case {'caret_surf' 'caret_topo' 'caret_coord'}
    ft_hastoolbox('gifti', 1);
    g = gifti(filename);
    if ~isfield(g, 'vertices') && strcmp(fileformat, 'caret_topo')
      try
        % do a clever guess by replacing topo with coord
        g2 = gifti(strrep(filename, '.topo.', '.coord.'));
        vertices  = ft_warp_apply(g2.mat, g2.vertices);
      catch
        vertices  = [];
      end
    else
      vertices  = ft_warp_apply(g.mat, g.vertices);
    end
    if ~isfield(g, 'faces') && strcmp(fileformat, 'caret_coord')
      try
        % do a clever guess by replacing topo with coord
        g2 = gifti(strrep(filename, '.coord.', '.topo.'));
        faces = g2.faces;
      catch
        faces = [];
      end
    else
      faces = g.faces;
    end
    
    shape.pos = vertices;
    shape.tri = faces;
    if isfield(g, 'cdata')
      shape.mom = g.cdata;
    end
    
    % check whether there is curvature info etc
    filename    = strrep(filename, '.surf.', '.shape.');
    [p,f,e]     = fileparts(filename);
    tok         = tokenize(f, '.');
    if length(tok)>2
      tmpfilename = strrep(filename, tok{3}, 'sulc');
      if exist(tmpfilename, 'file'), g = gifti(tmpfilename); shape.sulc = g.cdata; end
      if exist(strrep(tmpfilename, 'sulc', 'curvature'), 'file'),  g = gifti(strrep(tmpfilename, 'sulc', 'curvature')); shape.curv = g.cdata; end
      if exist(strrep(tmpfilename, 'sulc', 'thickness'), 'file'),  g = gifti(strrep(tmpfilename, 'sulc', 'thickness')); shape.thickness = g.cdata; end
      if exist(strrep(tmpfilename, 'sulc', 'atlasroi'),  'file'),  g = gifti(strrep(tmpfilename, 'sulc', 'atlasroi'));  shape.atlasroi  = g.cdata; end
		end
    
  case 'caret_spec'
    [spec, headerinfo] = read_caret_spec(filename);
    fn = fieldnames(spec)
    
    % concatenate the filenames that contain coordinates
    % concatenate the filenames that contain topologies
    coordfiles = {};
    topofiles  = {};
    for k = 1:numel(fn)
      if ~isempty(strfind(fn{k}, 'topo'))
        topofiles = cat(1,topofiles, spec.(fn{k}));
      end
      if ~isempty(strfind(fn{k}, 'coord'))
        coordfiles = cat(1,coordfiles, spec.(fn{k}));
      end
    end
    [selcoord, ok] = listdlg('ListString',coordfiles,'SelectionMode','single','PromptString','Select a file describing the coordinates');
    [seltopo, ok]  = listdlg('ListString',topofiles,'SelectionMode','single','PromptString','Select a file describing the topology');
    
    % recursively call ft_read_headshape
    tmp1 = ft_read_headshape(coordfiles{selcoord});
    tmp2 = ft_read_headshape(topofiles{seltopo});
    
    % quick and dirty sanity check to see whether the indexing of the
    % points in the topology matches the number of points
    if max(tmp2.tri(:))~=size(tmp1.pos,1)
      error('there''s a mismatch between the number of points used in the topology, and described by the coordinates');
    end
    
    shape.pos = tmp1.pos;
    shape.tri = tmp2.tri;
    
  case 'neuromag_mex'
    [co,ki,nu] = hpipoints(filename);
    fid = co(:,find(ki==1))';
    
    [junk, NZ] = max(fid(:,2));
    [junk, L]  = min(fid(:,1));
    [junk, R]  = max(fid(:,1));
    
    shape.fid.pos = fid([NZ L R], :);
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
    shape.pos=[src(1).rr(inuse1,:); src(2).rr(inuse2,:)];
    
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
    if isfield(src(1), 'use_tri_area')
      shape.area = [src(1).use_tri_area(:); src(2).use_tri_area(:)];
    end
    if isfield(src(1), 'use_tri_nn')
      shape.nn = [src(1).use_tri_nn; src(2).use_tri_nn];
    end
    shape.orig.pos = [src(1).rr; src(2).rr];
    shape.orig.tri = [src(1).tris; src(2).tris + src(1).np];
    shape.orig.inuse = [src(1).inuse src(2).inuse]';
    shape.orig.nn    = [src(1).nn; src(2).nn];
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
    
  case {'neuromag_fif' 'neuromag_mne'}
    
    orig = read_neuromag_hc(filename);
    switch coordsys
      case 'head'
        fidN=1;
        pntN=1;
        for i=1:size(orig.head.pos,1)
          if strcmp(orig.head.label{i}, 'LPA') || strcmp(orig.head.label{i}, 'Nasion') || strcmp(orig.head.label{i}, 'RPA')
            shape.fid.pos(fidN,1:3) = orig.head.pos(i,:);
            shape.fid.label{fidN} = orig.head.label{i};
            fidN = fidN + 1;
          else
            shape.pos(pntN,1:3) = orig.head.pos(i,:);
            shape.label{pntN} = orig.head.label{i};
            pntN = pntN + 1;
          end
        end
        shape.coordsys = orig.head.coordsys;
      case 'dewar'
        fidN=1;
        pntN=1;
        for i=1:size(orig.dewar.pos,1)
          if strcmp(orig.dewar.label{i}, 'LPA') || strcmp(orig.dewar.label{i}, 'Nasion') || strcmp(orig.dewar.label{i}, 'RPA')
            shape.fid.pos(fidN,1:3) = orig.dewar.pos(i,:);
            shape.fid.label{fidN} = orig.dewar.label{i};
            fidN = fidN + 1;
          else
            shape.pos(pntN,1:3) = orig.dewar.pos(i,:);
            shape.label{pntN} = orig.dewar.label{i};
            pntN = pntN + 1;
          end
        end
        shape.coordsys = orig.dewar.coordsys;
      otherwise
        error('incorrect coordinates specified');
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
      shape.fid.pos = cat(1, marker(1:5).meg_pos);
      sw_ind = [3 1 2];
      shape.fid.pos(1:3,:)= shape.fid.pos(sw_ind, :);
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
    shape.fid.pos = [];
    shape.fid.label = {};
    while ind < nr_items
      if strcmp(in_str{ind},'MEG:x=')
        shape.fid.pos = [shape.fid.pos; str2num(strtok(in_str{ind+1},[',','['])) ...
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
    
    shape.fid.pos(1:3,:)= shape.fid.pos(sw_ind, :);
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
      shape.fid.pos = [
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
      shape.pos = zeros(siz(:)');
      for i=1:siz(1)
        line = fgetl(fid);
        shape.pos(i,:) = sscanf(line, '%f');
      end
    end
    
    fclose(fid);
    
  case 'ply'
    [vert, face] = read_ply(filename);
    shape.pos = [vert.x vert.y vert.z];
    if isfield(vert, 'red') && isfield(vert, 'green') && isfield(vert, 'blue')
      shape.color = double([vert.red vert.green vert.blue])/255;
    end
    switch size(face,2)
      case 3
        shape.tri = face;
      case 4
        shape.tet = face;
      case 8
        shape.hex = face;
    end
    
  case 'polhemus_fil'
    [shape.fid.pos, shape.pos, shape.fid.label] = read_polhemus_fil(filename, 0);
    
  case 'polhemus_pos'
    [shape.fid.pos, shape.pos, shape.fid.label] = read_ctf_pos(filename);
    
  case 'spmeeg_mat'
    tmp = load(filename);
    if isfield(tmp.D, 'fiducials') && ~isempty(tmp.D.fiducials)
      shape = tmp.D.fiducials;
    else
      error('no headshape found in SPM EEG file');
    end
    
  case 'matlab'
    % determine which variables are contained in the file
    tmp = load(filename);
    
    if isfield(tmp, 'shape')
      shape = tmp.shape;
    elseif isfield(tmp, 'headshape')
      shape = tmp.headshape;
    elseif isfield(tmp, 'bnd')
      % the variable in the file is most likely a precomputed triangulation of some sort
      shape = tmp.bnd;
    elseif isfield(tmp, 'mesh')
      % the variable in the file is most likely a precomputed triangulation of some sort
      shape = tmp.mesh;
    elseif isfield(tmp, 'elec')
      % represent the electrodes as headshape
      tmp.elec        = ft_datatype_sens(tmp.elec);
      shape.fid.pos   = tmp.elec.chanpos;
      shape.fid.label = tmp.elec.label;
    elseif isfield(tmp, 'Vertices')
      % this applies to BrainStorm cortical meshes
      shape.pos = tmp.Vertices;
      % copy some optional fields over with a new name
      shape = copyfields(tmp, shape, {'Faces', 'Curvature', 'SulciMap'});
      shape = renamefields(shape, {'Faces', 'Curvature', 'SulciMap'}, {'tri', 'curv', 'sulc'});
    else
      error('no headshape found in MATLAB file');
    end
    
  case {'freesurfer_triangle_binary', 'freesurfer_quadrangle'}
    % the freesurfer toolbox is required for this
    ft_hastoolbox('freesurfer', 1);
    
    [pos, tri] = read_surf(filename);
    
    if min(tri(:)) == 0
      % start counting from 1
      tri = tri + 1;
    end
    shape.pos = pos;
    shape.tri = tri;
    
    % for the left and right
    [path,name,ext] = fileparts(filename);
    
    if strcmp(ext, '.inflated') % does the shift only for inflated surface
      if strcmp(name, 'lh')
        % assume freesurfer inflated mesh in mm, mni space
        % move the mesh a bit to the left, to avoid overlap with the right
        % hemisphere
        shape.pos(:,1) = shape.pos(:,1) - max(shape.pos(:,1)) - 10;
        
      elseif strcmp(name, 'rh')
        % id.
        % move the mesh a bit to the right, to avoid overlap with the left
        % hemisphere
        shape.pos(:,1) = shape.pos(:,1) - min(shape.pos(:,1)) + 10;
      end
    end
    
    if exist(fullfile(path, [name,'.sulc']), 'file'), shape.sulc = read_curv(fullfile(path, [name,'.sulc'])); end
    if exist(fullfile(path, [name,'.curv']), 'file'), shape.curv = read_curv(fullfile(path, [name,'.curv'])); end
    if exist(fullfile(path, [name,'.area']), 'file'), shape.area = read_curv(fullfile(path, [name,'.area'])); end
    if exist(fullfile(path, [name,'.thickness']), 'file'), shape.thickness = read_curv(fullfile(path, [name,'.thickness'])); end
    
  case 'stl'
    [pos, tri, nrm] = read_stl(filename);
    shape.pos = pos;
    shape.tri = tri;
    
    case 'obj'
      ft_hastoolbox('wavefront', 1);
        % Implemented for structure.io .obj thus far 
        obj = read_wobj(filename);
        shape.pos     = obj.vertices;
        shape.pos     = shape.pos - repmat(sum(shape.pos)/length(shape.pos),[length(shape.pos),1]); %centering vertices
        shape.tri     = obj.objects(2).data.vertices;
        if (~isempty(image) && exist('hasimage','var'))
            texture = obj.vertices_texture;
            
            %Refines the mesh and textures to increase resolution of the
            %colormapping
            for i = 1:1
                [shape.pos, shape.tri,texture] = refine(shape.pos,shape.tri,'banks',texture);              
            end 
            picture     = imread(image);
            color = uint8(zeros(length(shape.pos),3));
            for i=1:length(shape.pos)
                color(i,1:3) = picture(floor((1-texture(i,2))*length(picture)),1+floor(texture(i,1)*length(picture)),1:3);
            end
        end
        
  case 'vtk'
    [pos, tri] = read_vtk(filename);
    shape.pos = pos;
    shape.tri = tri;
    
  case 'off'
    [pos, plc] = read_off(filename);
    shape.pos  = pos;
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
      shape.pos = hdr.orig.Var.headshapedata;
    else
      error('the NetMEG file "%s" does not contain headshape data', filename);
    end
    
  case 'vista'
    ft_hastoolbox('simbio', 1);
    [nodes,elements,labels] = read_vista_mesh(filename);
    shape.pos     = nodes;
    if size(elements,2)==8
      shape.hex     = elements;
    elseif size(elements,2)==4
      shape.tet = elements;
    else
      error('unknown elements format')
    end
    % representation of data is compatible with ft_datatype_parcellation
    shape.tissue = zeros(size(labels));
    numlabels = size(unique(labels),1);
    shape.tissuelabel = {};
    for i = 1:numlabels
      ulabel = unique(labels);
      shape.tissue(labels == ulabel(i)) = i;
      shape.tissuelabel{i} = num2str(ulabel(i));
    end
    
  case 'tet'
    % the toolbox from Gabriel Peyre has a function for this
    ft_hastoolbox('toolbox_graph', 1);
    [vertex, face] = read_tet(filename);
    %     'vertex' is a '3 x nb.vert' array specifying the position of the vertices.
    %     'face' is a '4 x nb.face' array specifying the connectivity of the tet mesh.
    shape.pos = vertex';
    shape.tet = face';
    
  case 'tetgen_ele'
    % reads in the tetgen format and rearranges according to FT conventions
    % tetgen files also return a 'faces' field, which is not used here
    [p, f, x] = fileparts(filename);
    filename = fullfile(p, f); % without the extension
    IMPORT = importdata([filename '.ele'],' ',1);
    shape.tet = IMPORT.data(:,2:5);
    if size(IMPORT.data,2)==6
      labels = IMPORT.data(:,6);
      % representation of tissue type is compatible with ft_datatype_parcellation
      numlabels = size(unique(labels),1);
      ulabel    = unique(labels);
      shape.tissue      = zeros(size(labels));
      shape.tissuelabel = {};
      for i = 1:numlabels
        shape.tissue(labels == ulabel(i)) = i;
        shape.tissuelabel{i} = num2str(ulabel(i));
      end
    end
    IMPORT = importdata([filename '.node'],' ',1);
    shape.pos = IMPORT.data(:,2:4);
    
  case 'brainsuite_dfs'
    % this requires the readdfs function from the BrainSuite MATLAB utilities
    ft_hastoolbox('brainsuite', 1);
    
    dfs = readdfs(filename);
    % these are expressed in MRI dimensions
    shape.pos  = dfs.vertices;
    shape.tri  = dfs.faces;
    shape.unit = 'unkown';
    
    % the filename is something like 2467264.right.mid.cortex.svreg.dfs
    % whereas the corresponding MRI is 2467264.nii and might be gzipped
    [p, f, x] = fileparts(filename);
    while ~isempty(x)
      [junk, f, x] = fileparts(f);
    end
    
    if exist(fullfile(p, [f '.nii']), 'file')
      fprintf('reading accompanying MRI file "%s"\n', fullfile(p, [f '.nii']));
      mri = ft_read_mri(fullfile(p, [f '.nii']));
      transform = eye(4);
      transform(1:3,4) = mri.transform(1:3,4); % only use the translation
      shape.pos  = ft_warp_apply(transform, shape.pos);
      shape.unit = mri.unit;
    elseif exist(fullfile(p, [f '.nii.gz']), 'file')
      fprintf('reading accompanying MRI file "%s"\n', fullfile(p, [f '.nii']));
      mri = ft_read_mri(fullfile(p, [f '.nii.gz']));
      transform = eye(4);
      transform(1:3,4) = mri.transform(1:3,4); % only use the translation
      shape.pos  = ft_warp_apply(transform, shape.pos);
      shape.unit = mri.unit;
    else
      warning('could not find accompanying MRI file, returning vertices in voxel coordinates');
    end
    
  case 'brainvisa_mesh'
    % this requires the loadmesh function from the BrainVISA MATLAB utilities
    ft_hastoolbox('brainvisa', 1);
    [shape.pos, shape.tri, shape.nrm] = loadmesh(filename);
    shape.tri = shape.tri + 1; % they should be 1-offset, not 0-offset
    shape.unit = 'unkown';
    
    if exist([filename '.minf'], 'file')
      minffid = fopen([filename '.minf']);
      hdr=fgetl(minffid);
      tfm_idx = strfind(hdr,'''transformations'':') + 21;
      transform = sscanf(hdr(tfm_idx:end),'%f,',[4 4])';
      fclose(minffid);
      if ~isempty(transform)
        shape.pos = ft_warp_apply(transform, shape.pos);
        shape = rmfield(shape, 'unit'); % it will be determined later on, based on the size
      end
    end
    
    if isempty(transform)
      % the transformation was not present in the minf file, try to get it from the MRI
      
      % the filename is something like subject01_Rwhite_inflated_4d.mesh
      % and it is accompanied by subject01.nii
      [p, f, x] = fileparts(filename);
      f = tokenize(f, '_');
      f = f{1};
      
      if exist(fullfile(p, [f '.nii']), 'file')
        fprintf('reading accompanying MRI file "%s"\n', fullfile(p, [f '.nii']));
        mri = ft_read_mri(fullfile(p, [f '.nii']));
        shape.pos  = ft_warp_apply(mri.transform, shape.pos);
        shape.unit = mri.unit;
        transform = true; % used for feedback
      elseif exist(fullfile(p, [f '.nii.gz']), 'file')
        fprintf('reading accompanying MRI file "%s"\n', fullfile(p, [f '.nii.gz']));
        mri = ft_read_mri(fullfile(p, [f '.nii.gz']));
        shape.pos  = ft_warp_apply(mri.transform, shape.pos);
        shape.unit = mri.unit;
        transform = true; % used for feedback
      end
    end
    
    if isempty(transform)
      warning('cound not determine the coordinate transformation, returning vertices in voxel coordinates');
    end
    
  case 'brainvoyager_srf'
    [pos, tri, srf] = read_bv_srf(filename);
    shape.pos = pos;
    shape.tri = tri;
    
    % FIXME add details from srf if possible
    % FIXME do transform
    % FIXME remove vertices that are not in a triangle
    % FIXME add unit
  
	case 'besa_sfp'
		[lab, pos] = read_besa_sfp(filename, 0);
		shape.pos = pos;
		
		% assume that all non-'headshape' points are fiducial markers
		hs = strmatch('headshape', lab);
		lab(hs) = [];
		pos(hs, :) = [];
		shape.fid.label = lab;
		shape.fid.pos = pos;
		
  case 'asa_elc'
    elec = ft_read_sens(filename);
    
    shape.fid.pos   = elec.chanpos;
    shape.fid.label = elec.label;
    
    npos = read_asa(filename, 'NumberHeadShapePoints=', '%d');
    if ~isempty(npos) && npos>0
      origunit = read_asa(filename, 'UnitHeadShapePoints', '%s', 1);
      pos  = read_asa(filename, 'HeadShapePoints', '%f', npos, ':');
      
      pos = ft_scalingfactor(origunit, 'mm')*pos;
      
      shape.pos = pos;
    end
    
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
          shape.fid.pos   = elec.chanpos;
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
        headmodel = ft_read_vol(filename);
        if ~ft_voltype(headmodel, 'bem')
          error('skin surface can only be extracted from boundary element model');
        else
          if ~isfield(headmodel, 'skin')
            headmodel.skin = find_outermost_boundary(headmodel.bnd);
          end
          shape.pos = headmodel.bnd(headmodel.skin).pos;
          shape.tri = headmodel.bnd(headmodel.skin).tri; % also return the triangulation
          success = 1;
        end
      catch
        success = false;
      end % try
    end
    
    if ~success
      error('unknown fileformat "%s" for head shape information', fileformat);
    end
end % switch fileformat

if isfield(shape, 'label')
  % ensure that it is a column
  shape.label = shape.label(:);
end

if isfield(shape, 'fid') && isfield(shape.fid, 'label')
  % ensure that it is a column
  shape.fid.label = shape.fid.label(:);
end

% this will add the units to the head shape and optionally convert
if ~isempty(unit)
  shape = ft_convert_units(shape, unit);
else
  try
    % ft_convert_units will fail for triangle-only gifties.
    shape = ft_convert_units(shape);
  catch
  end
end

% ensure that vertex positions are given in pos, not in pnt
shape = fixpos(shape);
% ensure that the numerical arrays are represented in double precision and not as integers
shape = ft_struct2double(shape);

if (~isempty(image) && exist('hasimage','var'))
    shape.color = color;
end    
    
