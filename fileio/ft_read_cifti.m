function source = ft_read_cifti(filename, varargin)

% FT_READ_CIFTI reads geometry and functional data or connectivity from a
% cifti file and returns a structure according to FT_DATATYPE_SOURCE.
%
% Use as
%   source = ft_read_cifti(filename, ...)
% where additional input arguments should be specified as key-value pairs
% and can include
%   readdata         = boolean, can be false or true
%   cortexleft       = string, filename with left cortex (optional, default is automatic)
%   cortexright      = string, filename with right cortex (optional, default is automatic)
%   hemisphereoffset = number, amount in milimeter to move the two hemispheres apart from each other (default = 0)
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

readdata         = ft_getopt(varargin, 'readdata', []);   % the default depends on file size, see below
cortexleft       = ft_getopt(varargin, 'cortexleft', {});
cortexright      = ft_getopt(varargin, 'cortexright', {});
hemisphereoffset = ft_getopt(varargin, 'hemisphereoffset', 0); % in mm, move the two hemispheres apart from each other

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert from xmltree object to a generic structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numericAttributeTypes = {'NumberOfMatrices', 'AppliesToMatrixDimension', 'IndexOffset', 'IndexCount', 'SurfaceNumberOfNodes', 'VolumeDimensions', 'SurfaceNumberOfVertices', 'SeriesStart', 'SeriesStep', 'NumberOfSeriesPoints', 'SeriesExponent', 'Vertices', 'MeterExponent'};

Cifti             = struct(); % the parent of the XML tree, it only contains version info
MatrixIndicesMap  = [];
Parcel            = [];
NamedMap          = [];
Surface           = [];
Volume            = [];
BrainModel        = struct([]);

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
  
  uid_Volume = find(map,'/MatrixIndicesMap/Volume');
  % the following will fail if there are multiple volumes
  if ~isempty(uid_Volume)
    volume = branch(map, uid_Volume);
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
        Parcel(j).VoxelIndicesIJK = reshape(tmp, 3, [])' + 1; % one offset
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
            try, MatrixIndicesMap(i).BrainModel(j).NodeIndices = str2num(get(brainmodel, children(brainmodel, uid), 'value')) + 1; end
            try, MatrixIndicesMap(i).BrainModel(j).VertexIndices = MatrixIndicesMap(i).BrainModel(j).NodeIndices; end % convert from cifti-1
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
            error('unsupported version');
        end % switch version
        
        
      case 'CIFTI_MODEL_TYPE_VOXELS'
        tmp = str2num(get(brainmodel, children(brainmodel, find(brainmodel, '/BrainModel/VoxelIndicesIJK')), 'value'));
        MatrixIndicesMap(i).BrainModel(j).VoxelIndicesIJK = reshape(tmp, 3, [])' + 1;
        
      otherwise
        error('unsupported ModelType');
    end % switch
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
  % include the data in the output structure, note that the fieldname might be changed further down
  data = squeeze(reshape(voxdata, hdr.dim(2:end)));
end
fclose(fid);

% include the nifti-2 header in the output structure
hdr = hdr;

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
      
      %       IndexOffset         = [MatrixIndicesMap(i).BrainModel(:).IndexOffset];
      %       IndexCount          = [MatrixIndicesMap(i).BrainModel(:).IndexCount];
      %       ModelType           = nan(length(MatrixIndicesMap(i).BrainModel),1);
      %       ModelTypelabel      = cell(length(MatrixIndicesMap(i).BrainModel),1);
      %       BrainStructure      = nan(length(MatrixIndicesMap(i).BrainModel),1);
      %       BrainStructurelabel = cell(length(MatrixIndicesMap(i).BrainModel),1);
      
      %       tmp = cumsum([0 IndexCount]);
      %       if ~isequal(IndexOffset, tmp(1:end-1))
      %         % this happens in some of the example cifti1 files
      %         % and might be a bug in the actual format of the data in those files
      %         warning('inconsistency between IndexOffset and IndexCount');
      %       end
      
      % concatenate all brain models, this is consistent with parcels
      for j=1:length(MatrixIndicesMap(i).BrainModel)
        % FIXME
        % %         if strcmp(MatrixIndicesMap(i).BrainModel(end).ModelType, 'CIFTI_MODEL_TYPE_SURFACE') && ~isfield(MatrixIndicesMap(i).BrainModel(end), 'VertexIndices')
        % %           if isfield(MatrixIndicesMap(i).BrainModel(end), 'NodeIndices')
        % %             MatrixIndicesMap(i).BrainModel(end).VertexIndices = MatrixIndicesMap(i).BrainModel(end).NodeIndices; % convert from cifti-1
        % %           elseif (MatrixIndicesMap(i).BrainModel(end).IndexCount==MatrixIndicesMap(i).BrainModel(end).SurfaceNumberOfNodes)
        % %             MatrixIndicesMap(i).BrainModel(end).VertexIndices = 1:MatrixIndicesMap(i).BrainModel(end).IndexCount; % assume it spans the whole surface
        % %           end
        % %         end
        BrainModel = cat(1, BrainModel, MatrixIndicesMap(i).BrainModel(j));
      end % for all BrainModels
      
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
        end % for each brain structure
        
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
        
      end % for each parcel
      
      
      %       source.label = {Parcel(:).Name};
      %
      %       tmp1 = {};
      %       tmp2 = [];
      %       for j=1:numel(Parcel)
      %         tmp1{end+1} = Parcel(j).Name;
      %         tmp2(end+1) = sum(cellfun(@numel, Parcel(j).Vertices));
      %       end
      %       begpos = [0 cumsum(tmp2(1:end-1))] + 1;
      %       endpos = begpos + tmp2 - 1;
      %       for j=1:length(tmp1)
      %         sel = find(strcmp(tmp1, tmp1{j}), 1, 'first');
      %         Parcellation(begpos(j):endpos(j)) = sel;
      %         Parcellationlabel{sel} = tmp1{sel};
      %       end
      %
      %       tmp1 = {};
      %       tmp2 = [];
      %       for j=1:numel(Parcel)
      %         for k=1:numel(Parcel(j).BrainStructure)
      %           tmp1{end+1} = Parcel(j).BrainStructure{k};
      %           tmp2(end+1) = numel(Parcel(j).Vertices{k});
      %         end
      %       end
      %       begpos = [0 cumsum(tmp2(1:end-1))] + 1;
      %       endpos = begpos + tmp2 - 1;
      %       for j=1:length(tmp1)
      %         sel = find(strcmp(tmp1, tmp1{j}), 1, 'first');
      %         BrainStructure(begpos(j):endpos(j)) = sel;
      %         BrainStructurelabel{sel} = tmp1{sel}(17:end); % strip the first part, i.e. CIFTI_STRUCTURE_
      %       end
      %
      %       Cifti.dataIndex = {};
      %       Cifti.greynodeIndex = {};
      %       for j=1:max(Parcellation)
      %         for k=1:max(BrainStructure)
      %           Cifti.dataIndex{end+1}     = find(Parcellation==j & BrainStructure==k);
      %           Cifti.greynodeIndex{end+1} = Parcel(j).Vertices{k} + 1; % convert from zero to one-offset
      %         end
      %       end
      %
      %       hasbrainmodel = true;
      %       Cifti.pos = nan(numel(BrainStructure),3);
      
      
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
      end % switch
      
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
          error('unsupported TimeStepUnits');
      end
      
    otherwise
      error('unsupported IndicesMapToDataType');
  end % switch
end % fo length MatrixIndicesMap

dimord = dimord(~cellfun(@isempty, dimord));
source.dimord = sprintf('%s_', dimord{:});
source.dimord(end) = [];

if ~isempty(BrainModel)
  % this applies both to dense and parcellated formats
  
  dataIndex     = cell(size(BrainModel));
  greynodeIndex = cell(size(BrainModel));
  surfaceIndex  = nan(size(BrainModel)); % remains nan if it maps into the volume
  
  geombeg = [BrainModel.IndexOffset]+1;
  geomend = geombeg + [BrainModel.IndexCount] - 1;
  
  for i=1:numel(BrainModel)
    dataIndex{i} = geombeg(i):geomend(i);
    switch BrainModel(i).ModelType
      case 'CIFTI_MODEL_TYPE_SURFACE'
        try
          greynodeIndex{i} = BrainModel(i).VertexIndices;
        catch
          greynodeIndex{i} = 1:BrainModel(i).SurfaceNumberOfNodes;
        end
        
        sel = find(strcmp({Surface(:).BrainStructure}, BrainModel(i).BrainStructure));
        assert(numel(sel)==1);
        surfaceIndex(i) = sel;
        
      case 'CIFTI_MODEL_TYPE_VOXELS'
        greynodeIndex{i} = 1:BrainModel(i).IndexCount;
        surfaceIndex(i) = nan; % does not map onto surface
        
      otherwise
        error('unsupported ModelType "%s"', BrainModel(i).ModelType);
    end
  end
  
  
  %           case 'CIFTI_MODEL_TYPE_SURFACE'
  %             posbeg = geomCount + 1;
  %             posend = geomCount + MatrixIndicesMap(i).BrainModel(j).SurfaceNumberOfVertices;
  %             geomCount = geomCount + MatrixIndicesMap(i).BrainModel(j).SurfaceNumberOfVertices; % increment with the number of vertices in the (external) surface
  %             Cifti.pos(posbeg:posend,:) = nan;
  %             Cifti.dataIndex{j}         = IndexOffset(j) + (1:IndexCount(j));  % these are indices in the data
  %             Cifti.greynodeIndex{j}     = posbeg:posend;                       % these are indices in the greynodes (vertices or subcortical voxels)
  %             if isfield(MatrixIndicesMap(i).BrainModel(j), 'VertexIndices')
  %               % data is only present on a subset of vertices
  %               Cifti.greynodeIndex{j} = Cifti.greynodeIndex{j}(MatrixIndicesMap(i).BrainModel(j).VertexIndices+1);
  %             end
  %
  %           case 'CIFTI_MODEL_TYPE_VOXELS'
  %             posbeg = geomCount + 1;
  %             posend = geomCount + IndexCount(j);
  %             geomCount = geomCount + IndexCount(j); % increment with the number of vertices in the subcortical structure
  %
  %             Cifti.pos(posbeg:posend,:) = reshape(MatrixIndicesMap(i).BrainModel(j).VoxelIndicesIJK, 3, IndexCount(j))';
  %             Cifti.dataIndex{j}         = IndexOffset(j) + (1:IndexCount(j));  % these are indices in the data
  %             Cifti.greynodeIndex{j}     = posbeg:posend;                       % these are indices in the greynodes (vertices or subcortical voxels)
  %
  
  
  %
  %   ModelType(posbeg:posend) = j; % indexed representation, see ft_datatype_parcellation
  %   ModelTypelabel{j} = MatrixIndicesMap(i).BrainModel(j).ModelType;
  %
  %   BrainStructure(posbeg:posend) = j; % indexed representation, see ft_datatype_parcellation
  %   BrainStructurelabel{j} = MatrixIndicesMap(i).BrainModel(j).BrainStructure;
  %
  %   if ~isempty(regexp(BrainStructurelabel{j}, '^CIFTI_STRUCTURE_', 'once'))
  %     BrainStructurelabel{j} = BrainStructurelabel{j}(17:end); % strip the first part, i.e. CIFTI_STRUCTURE_
  %   end
  %
  %   brainordinate.pos                  = Cifti.pos;
  %   brainordinate.unit                 = 'mm';   % volume and surface should be in consistent units, gifti is defined in mm, wb_view also expects mm
  %   if ~isempty(Cifti.Volume)
  %     % this only applies to the voxel coordinates, not to surface vertices which are NaN
  %     brainordinate.pos        = ft_warp_apply(Cifti.Volume.Transform, brainordinate.pos+1);  % one offset
  %     brainordinate.pos        = brainordinate.pos .* (10^Cifti.Volume.MeterExponent);        % convert from native to meter
  %     brainordinate.pos        = brainordinate.pos .* (10^3);                                 % convert from meter to milimeter
  %     brainordinate.dim        = Cifti.Volume.VolumeDimensions;
  %     brainordinate.transform  = Cifti.Volume.Transform;
  %   end
  %   Ngreynodes = size(brainordinate.pos,1);
  %   if any(strcmp(dimord, 'pos'))
  %     % it is a dense source structure, i.e. represented with a position for each brainordinate
  %     brainordinate.brainstructure      = BrainStructure;
  %     brainordinate.brainstructurelabel = BrainStructurelabel;
  %     % copy all brainordinate information over to the main level
  %     source = copyfields(brainordinate, source, fieldnames(brainordinate));
  %   elseif any(strcmp(dimord, 'chan'))
  %     % it is a parcellated source structure, i.e. represented by one channel per parcel
  %     brainordinate.brainstructure      = BrainStructure;
  %     brainordinate.brainstructurelabel = BrainStructurelabel;
  %     brainordinate.parcellation        = Parcellation;
  %     brainordinate.parcellationlabel   = Parcellationlabel;
  %     % keep the mapping between positions and greynodes
  %     brainordinate.posindex = cat(2, Cifti.dataIndex{:});
  %     brainordinate.srfindex = cat(2, Cifti.greynodeIndex{:});
  %     % store the brainordinate information in a sub-structure
  %     source.brainordinate = brainordinate;
  %   end
  
  
  pos      = zeros(0,3);
  posIndex = zeros(0,1);
  
  % represent the position of all vertices of all surfaces, not only the ones used
  % in the brain models
  for i=1:numel(Surface)
    pos      = cat(1, pos, nan(Surface(i).SurfaceNumberOfVertices, 3));
    posIndex = cat(1, posIndex, i*ones(Surface(i).SurfaceNumberOfVertices, 1));
  end
  
  % it would be possible to also represent all voxels, but for efficiency we only
  % represent the ones that have data assigned to them
  if ~isempty(Volume)
    tmp       = ft_warp_apply(Volume.Transform, cat(1, BrainModel(isnan(surfaceIndex)).VoxelIndicesIJK));
    pos       = cat(1, pos, tmp);
    posIndex  = cat(1, posIndex, nan(size(tmp,1),1));
  end
  
  % the surfaces come before the voxels
  if ~isempty(Surface)
    voxeloffset = sum([Surface.SurfaceNumberOfVertices]);
  else
    voxeloffset  = 0;
  end
  
  for i=1:numel(BrainModel)
    if strcmp(BrainModel(i).ModelType, 'CIFTI_MODEL_TYPE_SURFACE')
      sel = find(strcmp({Surface(:).BrainStructure}, BrainModel(i).BrainStructure));
      greynodeOffset(i) = find(posIndex==sel, 1, 'first') - 1;
    else
      sel = strcmp({BrainModel(1:i-1).ModelType}, 'CIFTI_MODEL_TYPE_VOXELS');
      greynodeOffset(i) = voxeloffset + sum([BrainModel(sel).IndexCount]);
    end
    greynodeIndex{i} = greynodeIndex{i} + greynodeOffset(i);
  end
  
end % if brainmodel

brainordinate.brainstructure = zeros(size(pos,1),1);
brainordinate.brainstructurelabel = {};

for i=1:numel(BrainModel)
  indx = find(strcmp(brainordinate.brainstructurelabel, BrainModel(i).BrainStructure));
  if isempty(indx)
    indx = length(brainordinate.brainstructurelabel)+1;
  end
  brainordinate.brainstructurelabel{indx} = BrainModel(i).BrainStructure;
  brainordinate.brainstructure(greynodeIndex{i}) = indx;
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

% replace the long name with the short name
[dum, indx1, indx2] = intersect(brainordinate.brainstructurelabel, list1);
brainordinate.brainstructurelabel(indx1) = list2(indx2);

if ~isempty(Parcel)
  brainordinate.parcellation = zeros(size(pos,1),1);
  brainordinate.parcellationlabel = {};
  for i=1:numel(Parcel)
    brainordinate.parcellationlabel{i} = Parcel(i).Name;
    sel = find(strcmp(BrainModelParcelName, Parcel(i).Name));
    for j=1:numel(sel)
      brainordinate.parcellation(greynodeIndex{sel(j)}) = i;
    end
  end
end

if readdata
  if isempty(Parcel)
    % the data is dense, make it consistent with the graynode positions
    dataIndex     = [dataIndex{:}];
    greynodeIndex = [greynodeIndex{:}];
    Ngreynodes    = numel(greynodeIndex);
  else
    % the data is defined on parcels
    Ngreynodes    = length(Parcel);
    dataIndex     = 1:Ngreynodes;
    greynodeIndex = 1:Ngreynodes;
  end
  
  switch source.dimord
    case {'pos' 'chan'}
      [m, n] = size(data);
      if m>n
        dat = nan(Ngreynodes,n);
        dat(greynodeIndex(dataIndex),1) = data;
      else
        dat = nan(Ngreynodes,m);
        dat(greynodeIndex(dataIndex),:) = transpose(data);
      end
    case {'pos_pos' 'chan_chan'}
      dat = nan(Ngreynodes,Ngreynodes);
      dat(greynodeIndex(dataIndex),greynodeIndex(dataIndex)) = data;
    case {'pos_time' 'chan_time'}
      Ntime = size(data,2);
      dat = nan(Ngreynodes,Ntime);
      dat(greynodeIndex(dataIndex),:) = data;
    case 'time_pos'
      Ntime = size(data,1);
      dat = nan(Ngreynodes,Ntime);
      dat(greynodeIndex(dataIndex),:) = transpose(data);
      source.dimord = 'pos_time';
    case 'time_chan'
      Ntime = size(data,1);
      dat = nan(Ngreynodes,Ntime);
      dat(greynodeIndex(dataIndex),:) = transpose(data);
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

source = copyfields(Cifti, source, {'time', 'freq'});
source.hdr = hdr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try to get the geometrical information from the corresponding gifti files
% the following assumes HCP/WorkBench conventions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


tri = zeros(0,3);
for i=1:length(Surface)
  
  switch Surface(i).BrainStructure
    case 'CIFTI_STRUCTURE_CORTEX_LEFT'
      
      for j=1:length(Lfilelist)
        if exist(Lfilelist{j}, 'file')
          warning('reading CORTEX_LEFT from %s', Lfilelist{j});
          mesh = ft_read_headshape(Lfilelist{j}, 'unit', 'mm'); % volume and surface should be in consistent units, gifti is defined in mm, wb_view also expects mm
          pos(posIndex==i,:) = mesh.pnt;
          tri = cat(1, tri, mesh.tri + find(posIndex==i, 1, 'first') - 1);
          
          break
        end
      end % for
      
    case 'CIFTI_STRUCTURE_CORTEX_RIGHT'
      
      for j=1:length(Rfilelist)
        if exist(Rfilelist{j}, 'file')
          warning('reading CORTEX_RIGHT from %s', Rfilelist{j});
          mesh = ft_read_headshape(Rfilelist{j}, 'unit', 'mm'); % volume and surface should be in consistent units, gifti is defined in mm, wb_view also expects mm
          pos(posIndex==i,:) = mesh.pnt;
          tri = cat(1, tri, mesh.tri + find(posIndex==i, 1, 'first') - 1);
          break
        end
      end % for
      
    otherwise
      keyboard
  end
end

% if isfield(source, 'brainordinate')
%   % it is a parcellated source structure, i.e. represented by one channel per parcel
%   brainordinate = source.brainordinate;
% else
%   % it is a dense source structure, i.e. represented with a position for each brainordinate
%   brainordinate = source;
% end
%
% % update the geometrical description of the brainordinates
% if isfield(brainordinate, 'brainstructure')
%   % it contains information about anatomical structures, including cortical surfaces
%   BrainStructure      = brainordinate.brainstructure;
%   BrainStructurelabel = brainordinate.brainstructurelabel;
%
%   if all(ismember({'CORTEX_LEFT', 'CORTEX_RIGHT'}, BrainStructurelabel))
%     for i=1:length(Lfilelist)
%       Rfilename = fullfile(p, Rfilelist{i});
%
%       if  && exist(Rfilename, 'file')
%         warning('reading left hemisphere geometry from %s',  Lfilename);
%         meshL =
%         warning('reading right hemisphere geometry from %s',  Rfilename);
%         meshR = ft_read_headshape(Rfilename, 'unit', 'mm'); % volume and surface should be in consistent units, gifti is defined in mm, wb_view also expects mm
%
%         % move the two hemispheres apart from each other
%         meshL.pnt(:,1) = meshL.pnt(:,1) - hemisphereoffset;
%         meshR.pnt(:,1) = meshR.pnt(:,1) + hemisphereoffset;
%
%         if isfield(brainordinate, 'parcellation')
%           % in case of a model with CIFTI_MODEL_TYPE_PARCELLATION, the size of the geometry is not specified in the cifti file
%           % here we extend the positions to all brainordinates, including the ones that do not have data on them
%
%           % FIXME there are these two that I should have used
%           % <Surface BrainStructure="CIFTI_STRUCTURE_CORTEX_LEFT" SurfaceNumberOfVertices="32492"/>
%           % <Surface BrainStructure="CIFTI_STRUCTURE_CORTEX_RIGHT" SurfaceNumberOfVertices="32492"/>
%
%           % concatenate the vertices of both hemispheres
%           brainordinate.pos = [
%             meshL.pnt
%             meshR.pnt
%             ];
%
%           % concatenate the triangles of both hemispheres
%           brainordinate.tri = [
%             meshL.tri
%             meshR.tri + size(meshL.pnt,1); % the vertices of the right hemisphere have an offset
%             ];
%
%           indexL = find(brainordinate.brainstructure==find(strcmp(brainordinate.brainstructurelabel, 'CORTEX_LEFT')));
%           indexR = find(brainordinate.brainstructure==find(strcmp(brainordinate.brainstructurelabel, 'CORTEX_RIGHT')));
%           brainordinate.srfindex(indexR) = brainordinate.srfindex(indexR) + size(meshL.pnt,1); % the vertices of the right hemisphere have an offset
%
%           % extend the brainstructure and parcellation to all vertices in the concatenated hemispheres
%           tmp1 = zeros(1, size(brainordinate.pos,1));
%           tmp2 = zeros(1, size(brainordinate.pos,1));
%           tmp1(brainordinate.srfindex) = brainordinate.brainstructure(brainordinate.posindex);
%           tmp2(brainordinate.srfindex) = brainordinate.parcellation(brainordinate.posindex);
%           brainordinate.brainstructure = tmp1;
%           brainordinate.parcellation   = tmp2;
%           brainordinate = removefields(brainordinate, {'posindex', 'srfindex'}); % these are not needed any more
%
%         else
%           % in case of a model with CIFTI_MODEL_TYPE_SURFACE or CIFTI_MODEL_TYPE_VOLUME, the size of the geometry is specified in the cifti file
%           % the positions already represents all brainordinates, including the ones that do not have data on them
%
%           % insert the positions of the hemispheres on the correct place
%           indexL = find(brainordinate.brainstructure==find(strcmp(brainordinate.brainstructurelabel, 'CORTEX_LEFT')));
%           indexR = find(brainordinate.brainstructure==find(strcmp(brainordinate.brainstructurelabel, 'CORTEX_RIGHT')));
%
%           brainordinate.pos(indexL,:) = meshL.pnt;
%           brainordinate.pos(indexR,:) = meshR.pnt;
%
%           brainordinate.tri = [
%             indexL(meshL.tri)
%             indexR(meshR.tri)
%             ];
%
%         end % if isfield parcellation
%
%         break % only read a single pair of meshes
%       end
%     end
%
%   elseif ismember({'CORTEX'}, brainordinate.brainstructurelabel)
%     for i=1:length(Bfilelist)
%       Bfilename = fullfile(p, Bfilelist{i});
%
%       if exist(Bfilename, 'file')
%         warning('reading surface geometry from %s',  Bfilename);
%         meshB   = ft_read_headshape(Bfilename, 'unit', 'mm'); % volume and surface should be in consistent units, gifti is defined in mm, wb_view also expects mm
%         indexB  = find(brainordinate.brainstructure==find(strcmp(brainordinate.brainstructurelabel, 'CORTEX')));
%         brainordinate.pos(indexB,:) = meshB.pnt;
%         brainordinate.tri = indexB(meshB.tri);
%       end
%
%       break % only read a single mesh
%     end
%   end
% end

try,
  brainordinate.pos = pos;
end

try,
  if ~isempty(tri)
    brainordinate.tri = tri;
  end
end


if ~isempty(Volume)
  brainordinate.dim        = Volume.VolumeDimensions;
  brainordinate.transform  = Volume.Transform;
end

% copy the updated brainordinates back into the source structure
if isempty(Parcel)
  source = copyfields(brainordinate, source, fieldnames(brainordinate));
else
  % it is a parcellated source structure, i.e. represented by one channel per parcel
  source.brainordinate = brainordinate;
  source.label = {Parcel(:).Name};
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
