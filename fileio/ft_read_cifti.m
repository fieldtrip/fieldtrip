function source = ft_read_cifti(filename, varargin)

% FT_READ_CIFTI reads geometry and functional data or connectivity from a
% cifti file and returns a structure according to FT_DATATYPE_SOURCE.
%
% Use as
%   source = ft_read_cifti(filename, ...)
% where additional input arguments should be specified as key-value pairs
% and can include
%   readdata       = boolean, can be false or true
%
% See also FT_WRITE_CIFTI, READ_NIFTI2_HDR, WRITE_NIFTI2_HDR

% Copyright (C) 2013-2014, Robert Oostenveld
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

brainstructure  = ft_getopt(varargin, 'brainstructure'); % the default is determined further down
readdata        = ft_getopt(varargin, 'readdata', []);   % the default depends on file size, see below
cortexleft      = ft_getopt(varargin, 'cortexleft', {});
cortexright     = ft_getopt(varargin, 'cortexright', {});

% convert 'yes'/'no' into boolean
readdata = istrue(readdata);

% ensure that the external toolbox is present, this adds gifti/@xmltree
ft_hastoolbox('gifti', 1);

% read the header section
hdr = read_nifti2_hdr(filename);

% xml_offset = 540+12;
% xml_size   = hdr.vox_offset-xml_offset-8;

fid = fopen(filename, 'rb', hdr.endian);

% determine the file size, this is used to catch endian errors
fseek(fid, 0, 'eof');
filesize = ftell(fid);
fseek(fid, 0, 'bof');

% set the default for readdata
if isempty(readdata)
  if filesize>1e9
    warning('filesize>1GB, not reading data by default. Please specify the ''readdata'' option.');
    readdata = false;
  else
    readdata = true;
  end
end

fseek(fid, 540, 'bof');
hdrext = fread(fid, [1 4], 'int8');
if hdrext(1)~=1
  error('cifti requires a header extension');
end

% determine the size of the header extension
esize = fread(fid, 1, 'int32=>int32');
etype = fread(fid, 1, 'int32=>int32');

hdrsize = 540;
voxsize = filesize-hdr.vox_offset;
if esize>(filesize-hdrsize-voxsize)
  warning('the endianness of the header extension is inconsistent with the nifti-2 header');
  esize = swapbytes(esize);
  etype = swapbytes(etype);
end

if etype~=32 && etype~=swapbytes(int32(32)) % FIXME there is an endian problem
  error('the header extension type is not cifti');
end

% read the extension content, subtract the 8 bytes from esize and etype
xmldata = fread(fid, [1 esize-8], 'uint8=>char');

% the size of the extension must be an integer multiple of 16 bytes according to http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/extension.html
% consequently there might be some zero-padding at the end of the XML section
if any(xmldata==0)
  xmldata = xmldata(xmldata>0);
end

% write the xml section to a temporary file
xmlfile = 'debug.xml';
tmp = fopen(xmlfile, 'w');
fwrite(tmp, xmldata);
fclose(tmp);

% this requires the xmltree object from Guillaume Flandin
% see http://www.artefact.tk/software/matlab/xml/
% it is also included with the gifti toolbox
tree = xmltree(xmldata);

% convert from xmltree object to a generic MATLAB structure
cii = tree2struct(tree);

if readdata
  % read the voxel data section
  fseek(fid, hdr.vox_offset, 'bof');
  switch hdr.datatype
    case   2, [voxdata, nitemsread] = fread(fid, inf, 'uchar');   assert(nitemsread==prod(hdr.dim(2:end)), 'could not read all data');
    case   4, [voxdata, nitemsread] = fread(fid, inf, 'short');   assert(nitemsread==prod(hdr.dim(2:end)), 'could not read all data');
    case   8, [voxdata, nitemsread] = fread(fid, inf, 'int');     assert(nitemsread==prod(hdr.dim(2:end)), 'could not read all data');
    case  16, [voxdata, nitemsread] = fread(fid, inf, 'float');   assert(nitemsread==prod(hdr.dim(2:end)), 'could not read all data');
    case  64, [voxdata, nitemsread] = fread(fid, inf, 'double');  assert(nitemsread==prod(hdr.dim(2:end)), 'could not read all data');
    case 512, [voxdata, nitemsread] = fread(fid, inf, 'ushort');  assert(nitemsread==prod(hdr.dim(2:end)), 'could not read all data');
    case 768, [voxdata, nitemsread] = fread(fid, inf, 'uint');    assert(nitemsread==prod(hdr.dim(2:end)), 'could not read all data');
    otherwise, error('unsupported datatype');
  end
  % include the data in the output structure, note that the fieldname might be changed further down
  cii.data = squeeze(reshape(voxdata, hdr.dim(2:end)));
end
fclose(fid);

% include the nifti-2 header in the output structure
cii.hdr = hdr;

% convert to FieldTrip source representation, i.e. according to FT_DATATYPE_SOURCE and FT_DATATYPE_PARCELLATION
source = struct2source(cii);

% try to get the geometrical information from a corresponding gifti files
% the following assumes the HCP convention
[p, f, x] = fileparts(filename);
t = tokenize(f, '.');

subject  = 'unknown';
dataname = 'unknown';
geomodel = '';

if length(t)==2
  subject  = t{1};
  dataname = t{2};
elseif length(t)==3
  subject  = t{1};
  dataname = t{2};
  content  = t{3};
elseif length(t)==4
  subject  = t{1};
  dataname = t{2};
  geomodel = t{3};
  content  = t{4};
elseif length(t)==5
  subject  = t{1};
  dataname = [t{2} '.' t{3}];
  geomodel = t{4};
  content  = t{5};
else
  error('cannot parse file name');
end

% the surface anatomy is represented in an external file
% which can be difficult to match with the data

Lfilelist = {
  [subject '.L' '.midthickness'  '.' geomodel '.surf.gii']
  [subject '.L' '.pial'          '.' geomodel '.surf.gii']
  [subject '.L' '.white'         '.' geomodel '.surf.gii']
  [subject '.L' '.inflated'      '.' geomodel '.surf.gii']
  [subject '.L' '.very_inflated' '.' geomodel '.surf.gii']
  [subject '.L' '.sphere'        '.' geomodel '.surf.gii']
  [subject '.L' '.'              '.' geomodel '.surf.gii']
  [subject '.L' '.midthickness'               '.surf.gii']
  [subject '.L' '.pial'                       '.surf.gii']
  [subject '.L' '.white'                      '.surf.gii']
  [subject '.L' '.inflated'                   '.surf.gii']
  [subject '.L' '.very_inflated'              '.surf.gii']
  [subject '.L' '.sphere'                     '.surf.gii']
  [subject '.L'                               '.surf.gii']
  [subject '.CORTEX_LEFT'                     '.surf.gii']
  };

Rfilelist = {
  [subject '.R' '.midthickness'  '.' geomodel  '.surf.gii']
  [subject '.R' '.pial'          '.' geomodel  '.surf.gii']
  [subject '.R' '.white'         '.' geomodel  '.surf.gii']
  [subject '.R' '.inflated'      '.' geomodel  '.surf.gii']
  [subject '.R' '.very_inflated' '.' geomodel  '.surf.gii']
  [subject '.R' '.sphere'        '.' geomodel  '.surf.gii']
  [subject '.R' '.'              '.' geomodel  '.surf.gii']
  [subject '.R' '.midthickness'                '.surf.gii']
  [subject '.R' '.pial'                        '.surf.gii']
  [subject '.R' '.white'                       '.surf.gii']
  [subject '.R' '.inflated'                    '.surf.gii']
  [subject '.R' '.very_inflated'               '.surf.gii']
  [subject '.R' '.sphere'                      '.surf.gii']
  [subject '.R'                                '.surf.gii']
  [subject '.CORTEX_RIGHT'                     '.surf.gii']
  };

Bfilelist = {
  [subject '.midthickness'  '.' geomodel '.surf.gii']
  [subject '.pial'          '.' geomodel '.surf.gii']
  [subject '.white'         '.' geomodel '.surf.gii']
  [subject '.inflated'      '.' geomodel '.surf.gii']
  [subject '.very_inflated' '.' geomodel '.surf.gii']
  [subject '.sphere'        '.' geomodel '.surf.gii']
  [subject                  '.' geomodel '.surf.gii']
  [subject '.midthickness'               '.surf.gii']
  [subject '.pial'                       '.surf.gii']
  [subject '.white'                      '.surf.gii']
  [subject '.inflated'                   '.surf.gii']
  [subject '.very_inflated'              '.surf.gii']
  [subject '.sphere'                     '.surf.gii']
  [subject                               '.surf.gii']
  [subject '.CORTEX'                     '.surf.gii']
  };

Lfilelist = cat(1, cortexleft, Lfilelist);
Rfilelist = cat(1, cortexright, Rfilelist);

% set the default
if isempty(brainstructure) && isfield(source, 'brainstructure') && isfield(source, 'brainstructurelabel')
  brainstructure = 'brainstructure';
end

if isfield(source, brainstructure)
  % it contains information about anatomical structures, including cortical surfaces
  BrainStructure      = source.( brainstructure         );
  BrainStructurelabel = source.([brainstructure 'label']);
  
  if all(ismember({'CORTEX_LEFT', 'CORTEX_RIGHT'}, BrainStructurelabel))
    for i=1:length(Lfilelist)
      Lfilename = fullfile(p, Lfilelist{i});
      Rfilename = fullfile(p, Rfilelist{i});
      
      if exist(Lfilename, 'file') && exist(Rfilename, 'file')
        warning('reading left hemisphere geometry from %s',  Lfilename);
        meshL = ft_read_headshape(Lfilename, 'unit', 'mm'); % volume and surface should be in consistent units, gifti is defined in mm, wb_view also expects mm
        warning('reading right hemisphere geometry from %s',  Rfilename);
        meshR = ft_read_headshape(Rfilename, 'unit', 'mm'); % volume and surface should be in consistent units, gifti is defined in mm, wb_view also expects mm
        
        indexL = find(BrainStructure==find(strcmp(BrainStructurelabel, 'CORTEX_LEFT')));
        indexR = find(BrainStructure==find(strcmp(BrainStructurelabel, 'CORTEX_RIGHT')));
        
        source.pos(indexL,:) = meshL.pnt;
        source.pos(indexR,:) = meshR.pnt;
        
        source.tri = [
          indexL(meshL.tri)
          indexR(meshR.tri)
          ];
        
        break % only read a single pair of meshes
      end
    end
  elseif ismember({'CORTEX'}, BrainStructurelabel)
    for i=1:length(Bfilelist)
      Bfilename = fullfile(p, Bfilelist{i});
      
      if exist(Bfilename, 'file')
        warning('reading surface geometry from %s',  Bfilename);
        meshB   = ft_read_headshape(Bfilename, 'unit', 'mm'); % volume and surface should be in consistent units, gifti is defined in mm, wb_view also expects mm
        indexB  = find(BrainStructure==find(strcmp(BrainStructurelabel, 'CORTEX')));
        source.pos(indexB,:) = meshB.pnt;
        source.tri = indexB(meshB.tri);
      end
      
      break % only read a single mesh
    end
  end
end

if readdata
  if isfield(source, 'data')
    % rename the data field
    source.(fixname(dataname)) = source.data;
    source = rmfield(source, 'data');
  end
  
  % rename the datalabel field
  if isfield(source, 'datalabel')
    source.(fixname([dataname 'label'])) = source.datalabel;
    source = rmfield(source, 'datalabel');
  end
end

return % function ft_read_cifti

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the cifti XML section can be represented in various manners
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cifti = tree2struct(tree)

numericAttributeTypes = {'NumberOfMatrices', 'AppliesToMatrixDimension', 'IndexOffset', 'IndexCount', 'SurfaceNumberOfNodes', 'VolumeDimensions', 'SurfaceNumberOfVertices', 'SeriesStart', 'SeriesStep', 'NumberOfSeriesPoints', 'SeriesExponent', 'Vertices', 'MeterExponent'};

Cifti                  = struct(); % the parent of the XML tree, it only contains version info
Cifti.MatrixIndicesMap = struct(); % this is the interesting content

Volume    = [];
Parcel    = [];
NamedMap  = [];
Surface   = [];

attr = attributes(tree, 'get', 1);
if ~iscell(attr), attr = {attr}; end % treat one attribute just like multiple attributes
for j=1:numel(attr)
  if any(strcmp(attr{j}.key, numericAttributeTypes))
    Cifti.(attr{j}.key) = str2num(attr{j}.val);
  else
    Cifti.(attr{j}.key) = attr{j}.val;
  end
end

uid_MatrixIndicesMap = find(tree,'/CIFTI/Matrix/MatrixIndicesMap');
for i=1:length(uid_MatrixIndicesMap)
  map = branch(tree, uid_MatrixIndicesMap(i));
  
  % get the attributes of each map
  attr = attributes(map, 'get', 1);
  for j=1:numel(attr)
    if any(strcmp(attr{j}.key, numericAttributeTypes))
      MatrixIndicesMap(i).(attr{j}.key) = str2num(attr{j}.val);
    else
      MatrixIndicesMap(i).(attr{j}.key) = attr{j}.val;
    end
  end
  
  uid_Volume = find(tree,'/CIFTI/Matrix/MatrixIndicesMap/Volume');
  % the following will fail if there are multiple volumes
  if ~isempty(uid_Volume)
    volume = branch(tree, uid_Volume);
    attr = attributes(volume, 'get', 1); % there should only be one attribute here
    if ~iscell(attr), attr = {attr}; end % treat one attribute just like multiple attributes
    for j=1:numel(attr)
      if any(strcmp(attr{j}.key, numericAttributeTypes))
        Volume.(attr{j}.key) = str2num(attr{j}.val);
      else
        Volume.(attr{j}.key) = attr{j}.val;
      end
    end
    uid_Transform = find(volume,'/Volume/TransformationMatrixVoxelIndicesIJKtoXYZ');
    if ~isempty(uid_Transform)
      transform = branch(volume, uid_Transform);
      attr = attributes(transform, 'get', 1);
      if isstruct(attr), attr = {attr}; end % treat one attribute just like multiple attributes
      for j=1:numel(attr)
        if any(strcmp(attr{j}.key, numericAttributeTypes))
          Volume.(attr{j}.key) = str2num(attr{j}.val);
        else
          Volume.(attr{j}.key) = attr{j}.val;
        end
      end
      Volume.Transform = str2num(get(transform, 2, 'value'));
      Volume.Transform = reshape(Volume.Transform, [4 4])'; % it needs to be transposed
    end
  end
  
  uid_Surface = find(map, '/MatrixIndicesMap/Surface');
  if ~isempty(uid_Surface)
    for j=1:length(uid_Surface)
      surface = branch(map, uid_Surface(j));
      
      % get the attributes of each surface
      attr = attributes(surface, 'get', 1);
      if isstruct(attr), attr = {attr}; end % treat one attribute just like multiple attributes
      for k=1:numel(attr)
        if any(strcmp(attr{k}.key, numericAttributeTypes))
          Surface(j).(attr{k}.key) = str2num(attr{k}.val);
        else
          Surface(j).(attr{k}.key) = attr{k}.val;
        end
      end % for
      
    end
  end
  
  uid_Parcel = find(map, '/MatrixIndicesMap/Parcel');
  if ~isempty(uid_Parcel)
    for j=1:length(uid_Parcel)
      parcel = branch(map, uid_Parcel(j));
      
      % get the attributes of each parcel
      attr = attributes(parcel, 'get', 1);
      if isstruct(attr), attr = {attr}; end % treat one attribute just like multiple attributes
      for k=1:numel(attr)
        if any(strcmp(attr{k}.key, numericAttributeTypes))
          Parcel(j).(attr{k}.key) = str2num(attr{k}.val);
        else
          Parcel(j).(attr{k}.key) = attr{k}.val;
        end
      end % for
      
      uid_VoxelIndicesIJK = find(parcel, '/Parcel/VoxelIndicesIJK');
      if ~isempty(uid_VoxelIndicesIJK)
        Parcel(j).VoxelIndicesIJK = get(parcel, children(parcel, uid_VoxelIndicesIJK), 'value');
      else
        Parcel(j).VoxelIndicesIJK = [];
      end
      
      uid_Vertices = find(parcel, '/Parcel/Vertices');
      for k=1:length(uid_Vertices)
        vertices = branch(parcel, uid_Vertices(k));
        
        attr = attributes(vertices, 'get', 1);
        if isstruct(attr), attr = {attr}; end % treat one attribute just like multiple attributes
        for l=1:numel(attr)
          switch attr{l}.key
            case 'BrainStructure'
              Parcel(j).BrainStructure{k} = attr{l}.val;
          end
        end % for
        
        Parcel(j).Vertices{k} = get(vertices, children(vertices, find(vertices, 'Vertices')), 'value');
        
      end
    end
  end
  
  uid_NamedMap = find(map, '/MatrixIndicesMap/NamedMap');
  if ~isempty(uid_NamedMap)
    for j=1:length(uid_NamedMap)
      
      namedmap = branch(map, uid_NamedMap(j));
      NamedMap(j).MapName = get(namedmap, children(namedmap, find(namedmap, '/NamedMap/MapName')), 'value');
      
      uid_LabelTable = find(namedmap, '/NamedMap/LabelTable');
      for k=1:length(uid_LabelTable);
        labeltable = branch(namedmap, uid_LabelTable(k));
        uid_Label = find(labeltable, '/LabelTable/Label');
        for l=1:length(uid_Label)
          % there are also potentially intersting atributes here, but I don't know what to do with them
          NamedMap(j).LabelTable.Label{l} = get(labeltable, children(labeltable, uid_Label(l)), 'value');
          attr = attributes(branch(labeltable, uid_Label(l)), 'get', 1);
          if isstruct(attr), attr = {attr}; end % treat one attribute just like multiple attributes
          for m=1:numel(attr)
            switch attr{m}.key
              case 'Key'
                NamedMap(j).LabelTable.Key(l)   = str2double(attr{m}.val);
              case 'Red'
                NamedMap(j).LabelTable.Red(l)   = str2double(attr{m}.val);
              case 'Green'
                NamedMap(j).LabelTable.Green(l) = str2double(attr{m}.val);
              case 'Blue'
                NamedMap(j).LabelTable.Blue(l)  = str2double(attr{m}.val);
              case 'Alpha'
                NamedMap(j).LabelTable.Alpha(l) = str2double(attr{m}.val);
            end
          end
        end
      end % for
    end
  end
  
  uid_BrainModel = find(map, '/MatrixIndicesMap/BrainModel');
  for j=1:length(uid_BrainModel)
    brainmodel = branch(map, uid_BrainModel(j));
    
    % get the attributes of each model
    attr = attributes(brainmodel, 'get', 1);
    if isstruct(attr), attr = {attr}; end % treat one attribute just like multiple attributes
    for k=1:numel(attr)
      if any(strcmp(attr{k}.key, numericAttributeTypes))
        MatrixIndicesMap(i).BrainModel(j).(attr{k}.key) = str2num(attr{k}.val);
      else
        MatrixIndicesMap(i).BrainModel(j).(attr{k}.key) = attr{k}.val;
      end
    end % for
    
    switch MatrixIndicesMap(i).BrainModel(j).ModelType
      case 'CIFTI_MODEL_TYPE_SURFACE'
        switch Cifti.Version
          case {'1' '1.0'}
            uid = find(brainmodel, '/BrainModel/NodeIndices');
            try, MatrixIndicesMap(i).BrainModel(j).NodeIndices = str2num(get(brainmodel, children(brainmodel, uid), 'value')); end
          case {'2' '2.0'}
            uid = find(brainmodel, '/BrainModel/VertexIndices');
            try, MatrixIndicesMap(i).BrainModel(j).VertexIndices = str2num(get(brainmodel, children(brainmodel, uid), 'value')); end
          otherwise
            error('unsupported version');
        end % switch version
        
      case 'CIFTI_MODEL_TYPE_VOXELS'
        MatrixIndicesMap(i).BrainModel(j).VoxelIndicesIJK = str2num(get(brainmodel, children(brainmodel, find(brainmodel, '/BrainModel/VoxelIndicesIJK')), 'value'));
        
      otherwise
        error('unsupported ModelType');
    end % switch
  end % for each BrainModel
end % for each MatrixIndicesMap

% add the relevant sections to the main structure
Cifti.MatrixIndicesMap = MatrixIndicesMap;
Cifti.Volume           = Volume;
Cifti.NamedMap         = NamedMap;
Cifti.Surface          = Surface;
Cifti.Parcel           = Parcel;

return % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the cifti XML section can be represented in various manners
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function source = struct2source(Cifti)

dimord = cell(size(Cifti.MatrixIndicesMap));
hasbrainmodel = false;

% MatrixIndicesMap.IndicesMapToDataType can be
% CIFTI_INDEX_TYPE_BRAIN_MODELS The dimension represents one or more brain models.
% CIFTI_INDEX_TYPE_PARCELS      The dimension represents a parcellation scheme.
% CIFTI_INDEX_TYPE_SERIES       The dimension represents a series of regular samples.
% CIFTI_INDEX_TYPE_SCALARS      The dimension represents named scalar maps.
% CIFTI_INDEX_TYPE_LABELS       The dimension represents named label maps.

for i=1:length(Cifti.MatrixIndicesMap)
  switch Cifti.MatrixIndicesMap(i).IndicesMapToDataType
    case 'CIFTI_INDEX_TYPE_BRAIN_MODELS'
      dimord(Cifti.MatrixIndicesMap(i).AppliesToMatrixDimension+1) = {'pos'};
      
      IndexOffset         = [Cifti.MatrixIndicesMap(i).BrainModel(:).IndexOffset];
      IndexCount          = [Cifti.MatrixIndicesMap(i).BrainModel(:).IndexCount];
      ModelType           = nan(length(Cifti.MatrixIndicesMap(i).BrainModel),1);
      ModelTypelabel      = cell(length(Cifti.MatrixIndicesMap(i).BrainModel),1);
      BrainStructure      = nan(length(Cifti.MatrixIndicesMap(i).BrainModel),1);
      BrainStructurelabel = cell(length(Cifti.MatrixIndicesMap(i).BrainModel),1);
      
      tmp = cumsum([0 IndexCount]);
      if ~isequal(IndexOffset, tmp(1:end-1))
        % this happens in some of the example cifti1 files
        % and might be a bug in the actual format of the data in those files
        warning('inconsistency between IndexOffset and IndexCount');
      end
      
      % count the number of greynodes in all combined brain models
      geomCount = 0;
      Cifti.pos = nan(0,3);
      
      % concatenate all greynode positions
      for j=1:length(Cifti.MatrixIndicesMap(i).BrainModel)
        hasbrainmodel = true;
        
        switch Cifti.MatrixIndicesMap(i).BrainModel(j).ModelType
          case 'CIFTI_MODEL_TYPE_SURFACE'
            switch Cifti.Version
              case {'1' '1.0'}
                posbeg = geomCount + 1;
                posend = geomCount + Cifti.MatrixIndicesMap(i).BrainModel(j).SurfaceNumberOfNodes;
                geomCount = geomCount + Cifti.MatrixIndicesMap(i).BrainModel(j).SurfaceNumberOfNodes; % increment with the number of vertices in the (external) surface
                Cifti.pos(posbeg:posend,:) = nan;
                Cifti.dataIndex{j}         = IndexOffset(j) + (1:IndexCount(j));  % these are indices in the data
                Cifti.greynodeIndex{j}     = posbeg:posend;                       % these are indices in the greynodes (vertices or subcortical voxels)
                if isfield(Cifti.MatrixIndicesMap(i).BrainModel(j), 'NodeIndices')
                  % data is only present on a subset of vertices
                  Cifti.greynodeIndex{j} = Cifti.greynodeIndex{j}(Cifti.MatrixIndicesMap(i).BrainModel(j).NodeIndices+1);
                end
                
              case {'2' '2.0'}
                posbeg = geomCount + 1;
                posend = geomCount + Cifti.MatrixIndicesMap(i).BrainModel(j).SurfaceNumberOfVertices;
                geomCount = geomCount + Cifti.MatrixIndicesMap(i).BrainModel(j).SurfaceNumberOfVertices; % increment with the number of vertices in the (external) surface
                Cifti.pos(posbeg:posend,:) = nan;
                Cifti.dataIndex{j}         = IndexOffset(j) + (1:IndexCount(j));  % these are indices in the data
                Cifti.greynodeIndex{j}     = posbeg:posend;                       % these are indices in the greynodes (vertices or subcortical voxels)
                if isfield(Cifti.MatrixIndicesMap(i).BrainModel(j), 'VertexIndices')
                  % data is only present on a subset of vertices
                  Cifti.greynodeIndex{j} = Cifti.greynodeIndex{j}(Cifti.MatrixIndicesMap(i).BrainModel(j).VertexIndices+1);
                end
                
              otherwise
                error('unsupported version');
            end % switch version
            
            
          case 'CIFTI_MODEL_TYPE_VOXELS'
            posbeg = geomCount + 1;
            posend = geomCount + IndexCount(j);
            geomCount = geomCount + IndexCount(j); % increment with the number of vertices in the subcortical structure
            
            Cifti.pos(posbeg:posend,:) = reshape(Cifti.MatrixIndicesMap(i).BrainModel(j).VoxelIndicesIJK, 3, IndexCount(j))';
            Cifti.dataIndex{j}         = IndexOffset(j) + (1:IndexCount(j));  % these are indices in the data
            Cifti.greynodeIndex{j}     = posbeg:posend;                       % these are indices in the greynodes (vertices or subcortical voxels)
            
          otherwise
            error('unexpected ModelType');
        end % switch
        
        % perform a sanity check on the data and greynode indices
        assert(numel(Cifti.dataIndex{j})==numel(Cifti.greynodeIndex{j}));
        
        ModelType(posbeg:posend) = j; % indexed representation, see ft_datatype_parcellation
        ModelTypelabel{j} = Cifti.MatrixIndicesMap(i).BrainModel(j).ModelType;
        
        BrainStructure(posbeg:posend) = j; % indexed representation, see ft_datatype_parcellation
        BrainStructurelabel{j} = Cifti.MatrixIndicesMap(i).BrainModel(j).BrainStructure;
        
        if ~isempty(regexp(BrainStructurelabel{j}, '^CIFTI_STRUCTURE_', 'once'))
          BrainStructurelabel{j} = BrainStructurelabel{j}(17:end); % strip the first part
        end
        
      end % for all BrainModels
      
    case 'CIFTI_INDEX_TYPE_PARCELS'
      dimord(Cifti.MatrixIndicesMap(i).AppliesToMatrixDimension+1) = {'chan'};
      source.channel = {Cifti.Parcel(:).Name};
      
    case 'CIFTI_INDEX_TYPE_SERIES'
      % this only applies to cifti version 2
      switch Cifti.MatrixIndicesMap(i).SeriesUnit
        case 'SECOND'
          dimord(Cifti.MatrixIndicesMap(i).AppliesToMatrixDimension+1) = {'time'};
          Cifti.time = (((1:Cifti.MatrixIndicesMap(i).NumberOfSeriesPoints)-1) * Cifti.MatrixIndicesMap(i).SeriesStep + Cifti.MatrixIndicesMap(i).SeriesStart) * 10^Cifti.MatrixIndicesMap(i).SeriesExponent;
        case 'HZ'
          dimord(Cifti.MatrixIndicesMap(i).AppliesToMatrixDimension+1) = {'freq'};
          Cifti.freq = (((1:Cifti.MatrixIndicesMap(i).NumberOfSeriesPoints)-1) * Cifti.MatrixIndicesMap(i).SeriesStep + Cifti.MatrixIndicesMap(i).SeriesStart) * 10^Cifti.MatrixIndicesMap(i).SeriesExponent;
          % case 'METER'
          % case 'RADIAN'
        otherwise
          error('unsupported SeriesUnit');
      end % switch
      
    case 'CIFTI_INDEX_TYPE_SCALARS'
      dimord{Cifti.MatrixIndicesMap(i).AppliesToMatrixDimension+1} = []; % scalars are not explicitly represented
      for j=1:length(Cifti.NamedMap)
        Cifti.mapname{j} = fixname(Cifti.NamedMap(j).MapName);
      end
      
    case 'CIFTI_INDEX_TYPE_LABELS'
      dimord{Cifti.MatrixIndicesMap(i).AppliesToMatrixDimension+1} = []; % labels are not explicitly represented
      for j=1:length(Cifti.NamedMap)
        key = Cifti.NamedMap(j).LabelTable.Key;
        lab = Cifti.NamedMap(j).LabelTable.Label;
        sel = key>0;
        Cifti.labeltable{j}(key(sel)) = lab(sel);
        Cifti.mapname{j} = fixname(Cifti.NamedMap(j).MapName);
      end
      
    case 'CIFTI_INDEX_TYPE_TIME_POINTS'
      % this only applies to cifti-1, in cifti-2 this has been replaced by CIFTI_INDEX_TYPE_SERIES
      dimord(Cifti.MatrixIndicesMap(i).AppliesToMatrixDimension+1) = {'time'};
      switch Cifti.MatrixIndicesMap(i).TimeStepUnits
        case 'NIFTI_UNITS_SEC'
          Cifti.fsample = 1/str2double(Cifti.MatrixIndicesMap(i).TimeStep);
        otherwise
          % other units should be trivial to implement
          error('unsupported TimeStepUnits');
      end
      
    otherwise
      error('unsupported IndicesMapToDataType');
  end % switch
end

dimord = dimord(~cellfun(@isempty, dimord));
source.dimord = sprintf('%s_', dimord{:});
source.dimord(end) = [];

if hasbrainmodel
  source.pos                  = Cifti.pos;
  source.unit                 = 'mm';   % volume and surface should be in consistent units, gifti is defined in mm, wb_view also expects mm
  source.brainstructure       = BrainStructure;
  source.brainstructurelabel  = BrainStructurelabel;
  if ~isempty(Cifti.Volume)
    % this only applies to the voxel coordinates, not to surface vertices which are NaN
    source.pos        = ft_warp_apply(Cifti.Volume.Transform, source.pos+1);  % one offset
    source.pos        = source.pos .* (10^Cifti.Volume.MeterExponent);        % convert from native to meter
    source.pos        = source.pos .* (10^3);                                 % convert from meter to milimeter
    source.dim        = Cifti.Volume.VolumeDimensions;
    source.transform  = Cifti.Volume.Transform;
  end
  Ngreynodes = size(source.pos,1);
end

if isfield(Cifti, 'data')
  if hasbrainmodel
    % make the data consistent with the graynode positions
    dataIndex     = [Cifti.dataIndex{:}];
    greynodeIndex = [Cifti.greynodeIndex{:}];
  else
    % the data is defined on parcels
    Ngreynodes    = length(Cifti.Parcel);
    dataIndex     = 1:Ngreynodes;
    greynodeIndex = 1:Ngreynodes;
  end
  
  switch source.dimord
    case {'pos' 'chan'}
      [m, n] = size(Cifti.data);
      if m>n
        dat = nan(Ngreynodes,n);
        dat(greynodeIndex(dataIndex),1) = Cifti.data;
      else
        dat = nan(Ngreynodes,m);
        dat(greynodeIndex(dataIndex),:) = transpose(Cifti.data);
      end
    case {'pos_pos' 'chan_chan'}
      dat = nan(Ngreynodes,Ngreynodes);
      dat(greynodeIndex(dataIndex),greynodeIndex(dataIndex)) = Cifti.data;
    case {'pos_time' 'chan_time'}
      Ntime = size(Cifti.data,2);
      dat = nan(Ngreynodes,Ntime);
      dat(greynodeIndex(dataIndex),:) = Cifti.data;
    case 'time_pos'
      Ntime = size(Cifti.data,1);
      dat = nan(Ngreynodes,Ntime);
      dat(greynodeIndex(dataIndex),:) = transpose(Cifti.data);
      source.dimord = 'pos_time';
    case 'time_chan'
      Ntime = size(Cifti.data,1);
      dat = nan(Ngreynodes,Ntime);
      dat(greynodeIndex(dataIndex),:) = transpose(Cifti.data);
      source.dimord = 'chan_time';
    otherwise
      error('unsupported dimord');
  end % switch
  
  if isfield(Cifti, 'mapname') && length(Cifti.mapname)>1
    % use distict names if there are multiple scalars or labels
    for i=1:length(Cifti.mapname)
      fieldname = fixname(Cifti.mapname{i});
      source.(fieldname) = dat(:,i);
      if isfield(Cifti, 'labeltable')
        source.([fieldname 'label']) = Cifti.labeltable{i};
      end
    end
  else
    % the name of the data will be based on the filename
    source.data = dat;
  end
end % if data

source = copyfields(Cifti, source, {'hdr', 'time', 'freq'});

return % function
