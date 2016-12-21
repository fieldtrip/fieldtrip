function source = ft_read_cifti(filename, varargin)

% FT_READ_CIFTI read functional data or functional connectivity from a cifti-1 or
% cifti-2 file. The functional data can consist of a dense or a parcellated
% representation. The geometrical description of the brainordinates can consist of
% triangulated surfaces or voxels in a regular 3-D volumetric grid. If available,
% it also reads the geometrical description of the surfaces from the accompanying
% gifti files.
%
% Use as
%   data = ft_read_cifti(filename, ...)
%
% If the file contains a dense representation of functional data, the output data
% structure is organized according to the FT_DATATYPE_SOURCE or FT_DATATYPE_VOLUME
% definition.
%
% If the contains a parcellated representation of functional data, the output data
% structure is organized according to the FT_DATATYPE_TIMELOCK or FT_DATATYPE_FREQ
% definition. In addition, the description of the geometry wil be represented in a
% data.brainordinate field, which is organized according to the FT_DATATYPE_SOURCE
% or FT_DATATYPE_VOLUME definition.
%
% Any optional input arguments should come in key-value pairs and may include
%   'readdata'         = boolean, can be false or true (default depends on file size)
%   'readsurface'      = boolean, can be false or true (default = true)
%   'cortexleft'       = string, filename with left cortex (optional, default is automatic)
%   'cortexright'      = string, filename with right cortex (optional, default is automatic)
%   'hemisphereoffset' = number, amount in milimeter to move the hemispheres apart from each other (default = 0)
%   'mapname'          = string, 'field' to represent multiple maps separately, or 'array' to represent as array (default = 'field')
%   'debug'            = boolean, write a debug.xml file (default = false)
%
% See also FT_WRITE_CIFTI, FT_READ_MRI, FT_WRITE_MRI

% Known limitations and bugs
% - it will fail if multiple MatrixIndicesMap contain BrainModels
% - fibers (i.e. dfan and dfibersamp) are unsupported/untested
% - metadata is unsupported

% Copyright (C) 2013-2015, Robert Oostenveld
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

readdata         = ft_getopt(varargin, 'readdata', []);   % the default depends on file size, see below
readsurface      = ft_getopt(varargin, 'readsurface', true);
cortexleft       = ft_getopt(varargin, 'cortexleft', {});
cortexright      = ft_getopt(varargin, 'cortexright', {});
hemisphereoffset = ft_getopt(varargin, 'hemisphereoffset', 0); % in mm, move the two hemispheres apart from each other
debug            = ft_getopt(varargin, 'debug', false);
mapname          = ft_getopt(varargin, 'mapname', 'field');
dataformat       = ft_getopt(varargin, 'dataformat', []);

% convert 'yes'/'no' into boolean
readdata = istrue(readdata);

if ft_filetype(filename, 'compressed')
  % the file is compressed, unzip on the fly
  inflated = true;
  origfile = filename;
  filename = inflate_file(filename);
else
  inflated = false;
  origfile = filename;
end

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
    warning('Not reading data by default in case filesize>1GB. Please specify the ''readdata'' option.');
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

if debug
  try
    % write the xml section to a temporary file
    xmlfile = 'debug.xml';
    tmp = fopen(xmlfile, 'w');
    fwrite(tmp, xmldata);
    fclose(tmp);
  end
end

% ensure that the external toolbox is present, this adds gifti/@xmltree
ft_hastoolbox('gifti', 1);

% convert the character data to an xmltree object
tree = xmltree(xmldata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert the xmltree object to a generic structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numericAttributeTypes = {'NumberOfMatrices', 'AppliesToMatrixDimension', 'IndexOffset', 'IndexCount', 'SurfaceNumberOfNodes', 'VolumeDimensions', 'SurfaceNumberOfVertices', 'SeriesStart', 'SeriesStep', 'NumberOfSeriesPoints', 'SeriesExponent', 'Vertices', 'MeterExponent'};

Cifti             = struct(); % the parent of the XML tree, it only contains version info
MatrixIndicesMap  = [];
Parcel            = [];
NamedMap          = [];
Surface           = [];
Volume            = [];
BrainModel        = struct([]); % this will be constructed both for dense and parcellated files

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
  
  switch Cifti.Version
    case {'1' '1.0'}
      uid_Volume = find(tree,'/CIFTI/Matrix/Volume');
    case {'2' '2.0'}
      uid_Volume = find(map,'/MatrixIndicesMap/Volume');
  end
  % the following will fail if there are multiple volumes
  if ~isempty(uid_Volume)
    switch Cifti.Version
      case {'1' '1.0'}
        volume = branch(tree, uid_Volume);
      case {'2' '2.0'}
        volume = branch(map, uid_Volume);
    end
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
  
  uid_Surface = find(map,'/MatrixIndicesMap/Surface');
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
  
  uid_Parcel = find(map,'/MatrixIndicesMap/Parcel');
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
        tmp = str2num(get(parcel, children(parcel, uid_VoxelIndicesIJK), 'value'));
        Parcel(j).VoxelIndicesIJK = reshape(tmp, 3, [])' + 1; % transpose, one offset
      else
        Parcel(j).VoxelIndicesIJK = [];
      end
      
      Parcel(j).BrainStructure = {};
      
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
        
        Parcel(j).Vertices{k} = str2num(get(vertices, children(vertices, find(vertices, 'Vertices')), 'value')) + 1; % one offset
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
      end % for each LabelTable
    end % for each NamedMap
  end % if NamedMap
  
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
    end % for each attr
    
    switch MatrixIndicesMap(i).BrainModel(j).ModelType
      case 'CIFTI_MODEL_TYPE_SURFACE'
        switch Cifti.Version
          case {'1' '1.0'}
            uid = find(brainmodel, '/BrainModel/NodeIndices');
            % use VertexIndices rather than NodeIndices
            try, MatrixIndicesMap(i).BrainModel(j).VertexIndices = str2num(get(brainmodel, children(brainmodel, uid), 'value')) + 1; end
            % also copy the global surface information to a higher level
            Surface(end+1).BrainStructure          = MatrixIndicesMap(i).BrainModel(j).BrainStructure;
            Surface(end  ).SurfaceNumberOfVertices = MatrixIndicesMap(i).BrainModel(j).SurfaceNumberOfNodes;
            
          case {'2' '2.0'}
            uid = find(brainmodel, '/BrainModel/VertexIndices');
            try, MatrixIndicesMap(i).BrainModel(j).VertexIndices = str2num(get(brainmodel, children(brainmodel, uid), 'value')) + 1; end
            % also copy the global surface information to a higher level
            Surface(end+1).BrainStructure          = MatrixIndicesMap(i).BrainModel(j).BrainStructure;
            Surface(end  ).SurfaceNumberOfVertices = MatrixIndicesMap(i).BrainModel(j).SurfaceNumberOfVertices;
            
          otherwise
            error('unsupported cifti version');
        end % switch version
        
        
      case 'CIFTI_MODEL_TYPE_VOXELS'
        tmp = str2num(get(brainmodel, children(brainmodel, find(brainmodel, '/BrainModel/VoxelIndicesIJK')), 'value'));
        MatrixIndicesMap(i).BrainModel(j).VoxelIndicesIJK = reshape(tmp, 3, [])' + 1; % transpose, one offset
        
      otherwise
        error('unsupported ModelType');
    end % switch ModelType
  end % for each BrainModel
end % for each MatrixIndicesMap

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the voxel data section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if readdata
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
  % hdr.dim(1) is the number of dimensions
  % hdr.dim(2) is reserved for the x-dimension
  % hdr.dim(3) is reserved for the y-dimension
  % hdr.dim(4) is reserved for the z-dimension
  % hdr.dim(5) is reserved for the time-dimension
  % hdr.dim(6:8) are used for CIFTI
  voxdata = reshape(voxdata, hdr.dim(6:end));
end
fclose(fid);

if inflated
  % compressed file has been unzipped on the fly, clean up
  delete(filename);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert to FieldTrip source representation, i.e. according to FT_DATATYPE_SOURCE and FT_DATATYPE_PARCELLATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dimord = cell(size(MatrixIndicesMap));

% MatrixIndicesMap.IndicesMapToDataType can be
% CIFTI_INDEX_TYPE_BRAIN_MODELS The dimension represents one or more brain models.
% CIFTI_INDEX_TYPE_PARCELS      The dimension represents a parcellation scheme.
% CIFTI_INDEX_TYPE_SERIES       The dimension represents a series of regular samples.
% CIFTI_INDEX_TYPE_SCALARS      The dimension represents named scalar maps.
% CIFTI_INDEX_TYPE_LABELS       The dimension represents named label maps.

for i=1:length(MatrixIndicesMap)
  switch MatrixIndicesMap(i).IndicesMapToDataType
    case 'CIFTI_INDEX_TYPE_BRAIN_MODELS'
      dimord(MatrixIndicesMap(i).AppliesToMatrixDimension+1) = {'pos'};
      
      % concatenate all into an array of brain models
      for j=1:length(MatrixIndicesMap(i).BrainModel)
        BrainModel = cat(1, BrainModel, MatrixIndicesMap(i).BrainModel(j));
      end
      
    case 'CIFTI_INDEX_TYPE_PARCELS'
      dimord(MatrixIndicesMap(i).AppliesToMatrixDimension+1) = {'chan'};
      
      % construct an array of brain models
      IndexOffset = 0;
      BrainModelParcelName = {};
      
      for j=1:numel(Parcel)
        
        for k=1:numel(Parcel(j).BrainStructure)
          sel = strcmp({Surface(:).BrainStructure}, Parcel(j).BrainStructure{k});
          BrainModel(end+1).IndexOffset             = IndexOffset;
          BrainModel(end  ).IndexCount              = numel(Parcel(j).Vertices{k});
          BrainModel(end  ).BrainStructure          = Parcel(j).BrainStructure{k};
          BrainModel(end  ).ModelType               = 'CIFTI_MODEL_TYPE_SURFACE';
          BrainModel(end  ).SurfaceNumberOfVertices = Surface(sel).SurfaceNumberOfVertices;
          BrainModel(end  ).VertexIndices           = Parcel(j).Vertices{k};
          BrainModel(end  ).VoxelIndicesIJK         = [];
          IndexOffset = IndexOffset + numel(Parcel(j).Vertices{k});
          BrainModelParcelName{end+1} = Parcel(j).Name;
        end % for each BrainStructure
        
        if ~isempty(Parcel(j).VoxelIndicesIJK)
          BrainModel(end+1).IndexOffset             = IndexOffset;
          BrainModel(end  ).IndexCount              = size(Parcel(j).VoxelIndicesIJK,1);
          BrainModel(end  ).BrainStructure          = 'CIFTI_STRUCTURE_INVALID';
          BrainModel(end  ).ModelType               = 'CIFTI_MODEL_TYPE_VOXELS';
          BrainModel(end  ).SurfaceNumberOfVertices = [];
          BrainModel(end  ).VertexIndices           = [];
          BrainModel(end  ).VoxelIndicesIJK         = Parcel(j).VoxelIndicesIJK;
          IndexOffset = IndexOffset + size(Parcel(j).VoxelIndicesIJK,1);
          BrainModelParcelName{end+1} = Parcel(j).Name;
        end
        
      end % for each Parcel
      
    case 'CIFTI_INDEX_TYPE_SERIES'
      % this only applies to cifti version 2
      switch MatrixIndicesMap(i).SeriesUnit
        case 'SECOND'
          dimord(MatrixIndicesMap(i).AppliesToMatrixDimension+1) = {'time'};
          Cifti.time = (((1:MatrixIndicesMap(i).NumberOfSeriesPoints)-1) * MatrixIndicesMap(i).SeriesStep + MatrixIndicesMap(i).SeriesStart) * 10^MatrixIndicesMap(i).SeriesExponent;
        case 'HZ'
          dimord(MatrixIndicesMap(i).AppliesToMatrixDimension+1) = {'freq'};
          Cifti.freq = (((1:MatrixIndicesMap(i).NumberOfSeriesPoints)-1) * MatrixIndicesMap(i).SeriesStep + MatrixIndicesMap(i).SeriesStart) * 10^MatrixIndicesMap(i).SeriesExponent;
          % case 'METER'
          % case 'RADIAN'
        otherwise
          error('unsupported SeriesUnit');
      end % switch SeriesUnit
      
    case 'CIFTI_INDEX_TYPE_SCALARS'
      dimord{MatrixIndicesMap(i).AppliesToMatrixDimension+1} = []; % scalars are not explicitly represented
      for j=1:length(NamedMap)
        Cifti.mapname{j} = fixname(NamedMap(j).MapName);
      end
      
    case 'CIFTI_INDEX_TYPE_LABELS'
      dimord{MatrixIndicesMap(i).AppliesToMatrixDimension+1} = []; % labels are not explicitly represented
      for j=1:length(NamedMap)
        key = NamedMap(j).LabelTable.Key;
        lab = NamedMap(j).LabelTable.Label;
        sel = key>0;
        Cifti.labeltable{j}(key(sel)) = lab(sel);
        Cifti.mapname{j} = fixname(NamedMap(j).MapName);
      end
      
    case 'CIFTI_INDEX_TYPE_TIME_POINTS'
      % this only applies to cifti-1, in cifti-2 this has been replaced by CIFTI_INDEX_TYPE_SERIES
      dimord(MatrixIndicesMap(i).AppliesToMatrixDimension+1) = {'time'};
      switch MatrixIndicesMap(i).TimeStepUnits
        case 'NIFTI_UNITS_SEC'
          Cifti.fsample = 1/str2double(MatrixIndicesMap(i).TimeStep);
        otherwise
          % other units should be trivial to implement
          error('unsupported TimeStepUnits "%s"', MatrixIndicesMap(i).TimeStepUnits);
      end % switch TimeStepUnits
      
    otherwise
      error('unsupported IndicesMapToDataType');
  end % switch IndicesMapToDataType
end % for each MatrixIndicesMap

dimord = dimord(~cellfun(@isempty, dimord));
source.dimord = sprintf('%s_', dimord{:});
source.dimord(end) = [];

if ~isempty(BrainModel)
  % the BrainModel is constructed both for dense and parcellated data
  
  dataIndex           = cell(size(BrainModel));
  greynodeIndex       = cell(size(BrainModel));
  brainstructureIndex = cell(size(BrainModel));
  surfaceIndex        = nan(size(BrainModel)); % remains nan if it maps into the volume
  
  geombeg = [BrainModel.IndexOffset]+1;
  geomend = geombeg + [BrainModel.IndexCount] - 1;
  
  for i=1:numel(BrainModel)
    dataIndex{i} = geombeg(i):geomend(i);
    switch BrainModel(i).ModelType
      case 'CIFTI_MODEL_TYPE_SURFACE'
        
        switch Cifti.Version
          case {'1' '1.0'}
            try
              greynodeIndex{i}     = BrainModel(i).VertexIndices;
            catch
              greynodeIndex{i}     = 1:BrainModel(i).SurfaceNumberOfNodes;
            end
            brainstructureIndex{i} = 1:BrainModel(i).SurfaceNumberOfNodes;
          case {'2', '2.0'}
            greynodeIndex{i}       = BrainModel(i).VertexIndices;
            brainstructureIndex{i} = 1:BrainModel(i).SurfaceNumberOfVertices;
          otherwise
            error('unsupported cifti version');
        end % switch
        
        sel = find(strcmp({Surface(:).BrainStructure}, BrainModel(i).BrainStructure));
        assert(numel(sel)==1);
        surfaceIndex(i) = sel;
        
      case 'CIFTI_MODEL_TYPE_VOXELS'
        greynodeIndex{i}       = 1:BrainModel(i).IndexCount;
        brainstructureIndex{i} = 1:BrainModel(i).IndexCount;
        surfaceIndex(i) = nan; % does not map onto surface
        
      otherwise
        error('unsupported ModelType "%s"', BrainModel(i).ModelType);
    end
  end % for each BrainModel
  
  pos      = zeros(0,3);
  posIndex = zeros(0,1);
  
  % concatenate all vertices of all surfaces, including vertices that do not have data
  for i=1:numel(Surface)
    pos      = cat(1, pos, nan(Surface(i).SurfaceNumberOfVertices, 3));
    posIndex = cat(1, posIndex, i*ones(Surface(i).SurfaceNumberOfVertices, 1));
  end
  
  % it would be possible to represent all voxels, but for efficiency we only include voxel positions with data
  if ~isempty(Volume) && any(isnan(surfaceIndex))
    tmp       = ft_warp_apply(Volume.Transform, cat(1, BrainModel(isnan(surfaceIndex)).VoxelIndicesIJK));
    pos       = cat(1, pos, tmp);
    posIndex  = cat(1, posIndex, nan(size(tmp,1),1));
  end
  
  % the surface vertices come before the voxels
  if ~isempty(Surface)
    voxeloffset = sum([Surface.SurfaceNumberOfVertices]);
  else
    voxeloffset = 0;
  end
  
  greynodeOffset = nan(size(BrainModel));
  for i=1:numel(BrainModel)
    if strcmp(BrainModel(i).ModelType, 'CIFTI_MODEL_TYPE_SURFACE')
      sel = find(strcmp({Surface(:).BrainStructure}, BrainModel(i).BrainStructure));
      greynodeOffset(i) = find(posIndex==sel, 1, 'first') - 1;
    else
      sel = strcmp({BrainModel(1:i-1).ModelType}, 'CIFTI_MODEL_TYPE_VOXELS');
      greynodeOffset(i) = voxeloffset + sum([BrainModel(sel).IndexCount]);
    end
    % shift the greynodes to become consistent with the voxel data
    greynodeIndex{i} = greynodeIndex{i} + greynodeOffset(i);
    % shift the brainstructures to become consistent with the voxel data
    brainstructureIndex{i} = brainstructureIndex{i} + greynodeOffset(i);
  end
  
end % if BrainModel

brainordinate.brainstructure      = zeros(size(pos,1),1);
brainordinate.brainstructurelabel = {};

for i=1:numel(BrainModel)
  indx = find(strcmp(brainordinate.brainstructurelabel, BrainModel(i).BrainStructure));
  if isempty(indx)
    indx = length(brainordinate.brainstructurelabel)+1;
  end
  brainordinate.brainstructure(brainstructureIndex{i}) = indx;
  brainordinate.brainstructurelabel{indx} = BrainModel(i).BrainStructure;
end

list1 = {
  'CIFTI_STRUCTURE_CORTEX_LEFT'
  'CIFTI_STRUCTURE_CORTEX_RIGHT'
  'CIFTI_STRUCTURE_CEREBELLUM'
  'CIFTI_STRUCTURE_ACCUMBENS_LEFT'
  'CIFTI_STRUCTURE_ACCUMBENS_RIGHT'
  'CIFTI_STRUCTURE_ALL_GREY_MATTER'
  'CIFTI_STRUCTURE_ALL_WHITE_MATTER'
  'CIFTI_STRUCTURE_AMYGDALA_LEFT'
  'CIFTI_STRUCTURE_AMYGDALA_RIGHT'
  'CIFTI_STRUCTURE_BRAIN_STEM'
  'CIFTI_STRUCTURE_CAUDATE_LEFT'
  'CIFTI_STRUCTURE_CAUDATE_RIGHT'
  'CIFTI_STRUCTURE_CEREBELLAR_WHITE_MATTER_LEFT'
  'CIFTI_STRUCTURE_CEREBELLAR_WHITE_MATTER_RIGHT'
  'CIFTI_STRUCTURE_CEREBELLUM_LEFT'
  'CIFTI_STRUCTURE_CEREBELLUM_RIGHT'
  'CIFTI_STRUCTURE_CEREBRAL_WHITE_MATTER_LEFT'
  'CIFTI_STRUCTURE_CEREBRAL_WHITE_MATTER_RIGHT'
  'CIFTI_STRUCTURE_CORTEX'
  'CIFTI_STRUCTURE_DIENCEPHALON_VENTRAL_LEFT'
  'CIFTI_STRUCTURE_DIENCEPHALON_VENTRAL_RIGHT'
  'CIFTI_STRUCTURE_HIPPOCAMPUS_LEFT'
  'CIFTI_STRUCTURE_HIPPOCAMPUS_RIGHT'
  'CIFTI_STRUCTURE_INVALID'
  'CIFTI_STRUCTURE_OTHER'
  'CIFTI_STRUCTURE_OTHER_GREY_MATTER'
  'CIFTI_STRUCTURE_OTHER_WHITE_MATTER'
  'CIFTI_STRUCTURE_PALLIDUM_LEFT'
  'CIFTI_STRUCTURE_PALLIDUM_RIGHT'
  'CIFTI_STRUCTURE_PUTAMEN_LEFT'
  'CIFTI_STRUCTURE_PUTAMEN_RIGHT'
  'CIFTI_STRUCTURE_THALAMUS_LEFT'
  'CIFTI_STRUCTURE_THALAMUS_RIGHT'
  };

list2 = {
  'CORTEX_LEFT'
  'CORTEX_RIGHT'
  'CEREBELLUM'
  'ACCUMBENS_LEFT'
  'ACCUMBENS_RIGHT'
  'ALL_GREY_MATTER'
  'ALL_WHITE_MATTER'
  'AMYGDALA_LEFT'
  'AMYGDALA_RIGHT'
  'BRAIN_STEM'
  'CAUDATE_LEFT'
  'CAUDATE_RIGHT'
  'CEREBELLAR_WHITE_MATTER_LEFT'
  'CEREBELLAR_WHITE_MATTER_RIGHT'
  'CEREBELLUM_LEFT'
  'CEREBELLUM_RIGHT'
  'CEREBRAL_WHITE_MATTER_LEFT'
  'CEREBRAL_WHITE_MATTER_RIGHT'
  'CORTEX'
  'DIENCEPHALON_VENTRAL_LEFT'
  'DIENCEPHALON_VENTRAL_RIGHT'
  'HIPPOCAMPUS_LEFT'
  'HIPPOCAMPUS_RIGHT'
  'INVALID'
  'OTHER'
  'OTHER_GREY_MATTER'
  'OTHER_WHITE_MATTER'
  'PALLIDUM_LEFT'
  'PALLIDUM_RIGHT'
  'PUTAMEN_LEFT'
  'PUTAMEN_RIGHT'
  'THALAMUS_LEFT'
  'THALAMUS_RIGHT'
  };

% replace the long name with the short name, i.e remove 'CIFTI_STRUCTURE_' where applicable
[dum, indx1, indx2] = intersect(brainordinate.brainstructurelabel, list1);
brainordinate.brainstructurelabel(indx1) = list2(indx2);

if ~isempty(Parcel)
  brainordinate.parcellation = zeros(size(pos,1),1);
  brainordinate.parcellationlabel = {};
  for i=1:numel(Parcel)
    brainordinate.parcellationlabel{i} = Parcel(i).Name;
    sel = find(strcmp(BrainModelParcelName, Parcel(i).Name));
    for j=1:numel(sel)
      brainordinate.parcellation(greynodeIndex{sel(j)}) = i; % FIXME should this be greynodeIndex or brainstructureIndex?
    end
  end
end

if readdata
  if isempty(Parcel)
    % the data is dense, make it consistent with the graynode positions
    dataIndex     = [dataIndex{:}];
    greynodeIndex = [greynodeIndex{:}];
    Ngreynodes    = size(pos,1);
  else
    % the data is defined on parcels
    Ngreynodes    = length(Parcel);
    dataIndex     = 1:Ngreynodes;
    greynodeIndex = 1:Ngreynodes;
  end
  
  switch source.dimord
    % the following representations are directly consistent with FieldTrip
    case {'pos' 'chan'}
      [m, n] = size(voxdata);
      if m>n
        dat = nan(Ngreynodes,n);
        dat(greynodeIndex(dataIndex),:) = voxdata;
      else
        dat = nan(Ngreynodes,m);
        dat(greynodeIndex(dataIndex),:) = transpose(voxdata);
      end
    case {'pos_time' 'chan_time'}
      Ntime = size(voxdata,2);
      dat = nan(Ngreynodes,Ntime);
      dat(greynodeIndex(dataIndex),:) = voxdata;
    case {'pos_pos' 'chan_chan'}
      dat = nan(Ngreynodes,Ngreynodes);
      dat(greynodeIndex(dataIndex),greynodeIndex(dataIndex)) = voxdata;
    case {'pos_pos_time' 'chan_chan_time'}
      Ntime = size(voxdata,3);
      dat = nan(Ngreynodes,Ngreynodes,Ntime);
      dat(greynodeIndex(dataIndex),greynodeIndex(dataIndex),:) = voxdata;
      
      % the following representations need to be transposed to be consistent with FieldTrip
    case 'time_pos'
      Ntime = size(voxdata,1);
      dat = nan(Ngreynodes,Ntime);
      dat(greynodeIndex(dataIndex),:) = transpose(voxdata);
      source.dimord = 'pos_time';
    case 'time_chan'
      Ntime = size(voxdata,1);
      dat = nan(Ngreynodes,Ntime);
      dat(greynodeIndex(dataIndex),:) = transpose(voxdata);
      source.dimord = 'chan_time';
    case 'time_pos_pos'
      Ntime = size(voxdata,1);
      dat = nan(Ngreynodes,Ngreynodes,Ntime);
      dat(greynodeIndex(dataIndex),greynodeIndex(dataIndex),:) = permute(voxdata, [2 3 1]);
      source.dimord = 'pos_pos_time';
    case 'time_chan_chan'
      Ntime = size(voxdata,1);
      dat = nan(Ngreynodes,Ngreynodes,Ntime);
      dat(greynodeIndex(dataIndex),greynodeIndex(dataIndex),:) = permute(voxdata, [2 3 1]);
      source.dimord = 'chan_chan_time';
      
    otherwise
      error('unsupported dimord %s', source.dimord);
  end % switch
  
  if isfield(Cifti, 'mapname') && isfield(Cifti, 'labeltable') && strcmp(mapname, 'array')
    allthesame = true;
    for i=2:length(Cifti.labeltable)
      allthesame = allthesame && isequal(Cifti.labeltable{1}, Cifti.labeltable{i});
    end
    if allthesame
      warning('using the same labels for all maps in the array');
      source.datalabel = Cifti.labeltable{1};
      Cifti = rmfield(Cifti, 'labeltable');
    else
      error('multiple maps cannot be represented as array in the presence of different labeltables');
    end
  end
  
  if isfield(Cifti, 'mapname') && (length(Cifti.mapname)>1 || isfield(Cifti, 'labeltable'))
    switch mapname
      case 'field'
        % use distict names if there are multiple scalars or labels
        for i=1:length(Cifti.mapname)
          fieldname = Cifti.mapname{i};
          if isfield(Cifti, 'labeltable')
            if length(fieldname)>58
              % truncate it, needed to be able to append 'label' to the end
              fieldname = fieldname(1:58);
              % append 'label' to the end
              source.([fieldname 'label']) = Cifti.labeltable{i};
            end
          end
          source.(fieldname) = dat(:,i);
        end
      case 'array'
        source.mapname = {NamedMap.MapName}; % keep the original names, not the field names
        source.data    = dat;
        source.dimord  = [source.dimord '_mapname'];
      otherwise
        error('incorrect specification of mapname "%s"', mapname);
    end % switch mapname
  else
    % the name of the data will be based on the filename
    source.data = dat;
  end
end % if readdata

source = copyfields(Cifti, source, {'time', 'freq'});
source.hdr = hdr;
source.unit = 'mm'; % per definition

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try to get the geometrical information from the corresponding gifti files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use the filename prior to decompression
filename = origfile;

[p, f, x] = fileparts(filename);
t = tokenize(f, '.');

subject  = 'unknown';
dataname = 'unknown';
geomodel = '';

% the following assumes HCP/WorkBench/Caret file naming conventions
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
  warning('cannot parse file name');
end

if readsurface
  % construct a list of possible file names for the surface geometry
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
  
  % assume that the surface files are in the same directory as the cifti file
  for i=1:numel(Lfilelist)
    Lfilelist{i} = fullfile(p, Lfilelist{i});
  end
  for i=1:numel(Rfilelist)
    Rfilelist{i} = fullfile(p, Rfilelist{i});
  end
  for i=1:numel(Bfilelist)
    Bfilelist{i} = fullfile(p, Bfilelist{i});
  end
  
  Lfilelist = cat(1, cortexleft, Lfilelist);
  Rfilelist = cat(1, cortexright, Rfilelist);
  
  tri = zeros(0,3);
  for i=1:length(Surface)
    
    switch Surface(i).BrainStructure
      case 'CIFTI_STRUCTURE_CORTEX_LEFT'
        for j=1:length(Lfilelist)
          if exist(Lfilelist{j}, 'file')
            fprintf('reading CORTEX_LEFT surface from %s\n', Lfilelist{j});
            mesh = ft_read_headshape(Lfilelist{j}, 'unit', 'mm'); % volume and surface should be in consistent units, gifti is defined in mm, wb_view also expects mm
            mesh.pos(:,1) = mesh.pos(:,1) - hemisphereoffset;
            pos(posIndex==i,:) = mesh.pos;
            tri = cat(1, tri, mesh.tri + find(posIndex==i, 1, 'first') - 1);
            break
          end
        end % for each Lfilelist
        
      case 'CIFTI_STRUCTURE_CORTEX_RIGHT'
        for j=1:length(Rfilelist)
          if exist(Rfilelist{j}, 'file')
            fprintf('reading CORTEX_RIGHT surface from %s\n', Rfilelist{j});
            mesh = ft_read_headshape(Rfilelist{j}, 'unit', 'mm'); % volume and surface should be in consistent units, gifti is defined in mm, wb_view also expects mm
            mesh.pos(:,1) = mesh.pos(:,1) + hemisphereoffset;
            pos(posIndex==i,:) = mesh.pos;
            tri = cat(1, tri, mesh.tri + find(posIndex==i, 1, 'first') - 1);
            break
          end
        end % for each Rfilelist
        
      otherwise
        for j=1:length(Bfilelist)
          if exist(Bfilelist{j}, 'file')
            fprintf('reading %s surface from %s\n', Surface(i).BrainStructure(17:end), Bfilelist{j});
            mesh = ft_read_headshape(Bfilelist{j}, 'unit', 'mm'); % volume and surface should be in consistent units, gifti is defined in mm, wb_view also expects mm
            pos(posIndex==i,:) = mesh.pos;
            tri = cat(1, tri, mesh.tri + find(posIndex==i, 1, 'first') - 1);
            break
          end
        end % for each Bfilelist
        
    end % switch BrainStructure
  end
else
  tri = [];
end % if readsurface

% add the vertex and voxel positions
brainordinate.pos = pos;
if ~isempty(tri)
  % add the surface triangulations
  brainordinate.tri = tri;
end

if ~isempty(Volume)
  brainordinate.dim        = Volume.VolumeDimensions;
  brainordinate.transform  = Volume.Transform;
end

if isempty(Parcel)
  % copy the geometrical description of the brainordinates into the main structure
  source = copyfields(brainordinate, source, fieldnames(brainordinate));
else
  % it is a parcellated source structure, i.e. represented by one channel per parcel
  % copy the geometrical description of the brainordinates into a sub-structure
  source.brainordinate = brainordinate;
  source.label = {Parcel(:).Name}';
end

% haslabeltable = false;
% if ~isempty(NamedMap)
%   % the following assumes a single NamedMap
%   if isfield(NamedMap, 'LabelTable')
%     % use the key-label combination in the label table
%     haslabeltable    = true;
%     key              = NamedMap.LabelTable.Key;
%     source.datalabel = NamedMap.LabelTable.Label(:);
%   end
% end

if readdata
  if isfield(source, 'data')
    % rename the data field
    source.(fixname(dataname)) = source.data;
    source = rmfield(source, 'data');
%     % adopt FT convention for parcel-to-label mapping
%     if haslabeltable
%       tempdata = nan+zeros(size(source.data));
%       for k = 1:numel(key)
%         tempdata(source.data==key(k)) = k;
%       end
%       source.data = tempdata;
%     end
%     source = rmfield(source, 'data');
  end
  
  % rename the datalabel field
  if isfield(source, 'datalabel')
    source.(fixname([dataname 'label'])) = source.datalabel;
    source = rmfield(source, 'datalabel');
  end
  
end % if readdata
