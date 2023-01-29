classdef ProbeClass < matlab.mixin.Copyable
    
    % SNIRF-spec class properties
    properties
        wavelengths
        wavelengthsEmission
        sourcePos2D
        detectorPos2D
        landmarkPos2D
        sourcePos3D
        detectorPos3D
        frequencies
        timeDelays
        timeDelayWidths
        momentOrder
        correlationTimeDelays
        correlationTimeDelayWidths
        sourceLabels
        detectorLabels
        landmarkLabels
    end
    
    % Non-SNIRF class properties
    properties (Access = private)
        filename
        fileformat
    end
    
    
    methods
        
        % -------------------------------------------------------
        function obj = ProbeClass(varargin)
            % Set class properties not part of the SNIRF format
            obj.fileformat = 'hdf5';
            
            % Set SNIRF fomat properties
            if nargin>0
                if isstruct(varargin{1})
                    SD = varargin{1};
                    obj.wavelengths = SD.Lambda;
                    obj.wavelengthsEmission  = [];
                    if size(SD.SrcPos, 2) == 3 & SD.SrcPos(1, 3) ~= 0
                        obj.sourcePos3D  = SD.SrcPos;
                        obj.detectorPos3D  = SD.DetPos;
                        obj.sourcePos2D  = [];
                        obj.detectorPos2D  = [];
                    else
                        obj.sourcePos2D  = SD.SrcPos;
                        obj.detectorPos2D  = SD.DetPos;
                        obj.sourcePos3D  = [];
                        obj.detectorPos3D  = [];
                    end
                    obj.frequencies  = 1;
                    obj.timeDelays  = 0;
                    obj.timeDelayWidths  = 0;
                    obj.momentOrder = [];
                    obj.correlationTimeDelays = 0;
                    obj.correlationTimeDelayWidths = 0;
                    for ii=1:size(SD.SrcPos)
                        obj.sourceLabels{ii} = ['S',num2str(ii)];
                    end
                    for ii=1:size(SD.DetPos)
                        obj.detectorLabels{ii} = ['D',num2str(ii)];
                    end
                elseif ischar(varargin{1})
                    obj.filename = varargin{1};
                    obj.Load(varargin{1});
                end
            else
                obj.wavelengths          = [];
                obj.wavelengthsEmission  = [];
                obj.sourcePos2D  = [];
                obj.detectorPos2D  = [];
                obj.sourcePos3D  = [];
                obj.detectorPos3D  = [];
                obj.frequencies  = 1;
                obj.timeDelays  = 0;
                obj.timeDelayWidths  = 0;
                obj.momentOrder = [];
                obj.correlationTimeDelays = 0;
                obj.correlationTimeDelayWidths = 0;
                obj.sourceLabels = {};
                obj.detectorLabels = {};
            end
        end

        
        
        % -------------------------------------------------------
        function ForwardCompatibility(obj)
            if size(obj.sourcePos2D,2)<3
                obj.sourcePos2D       = [obj.sourcePos2D, zeros(size(obj.sourcePos2D,1), 1)];
            end
            if size(obj.detectorPos2D,2)<3
                obj.detectorPos2D     = [obj.detectorPos2D, zeros(size(obj.detectorPos2D,1), 1)];
            end
        end

        
        
        % -------------------------------------------------------
        function BackwardCompatibility(obj)
            if isempty(obj.sourcePos2D)
                obj.sourcePos2D   = HDF5_DatasetLoad(gid, 'sourcePos', [], '2D');
            end
            if isempty(obj.detectorPos2D)
                obj.detectorPos2D = HDF5_DatasetLoad(gid, 'detectorPos', [], '2D');
            end
            if isempty(obj.landmarkPos2D)
                obj.landmarkPos2D = HDF5_DatasetLoad(gid, 'landmarkPos', [], '2D');
            end
        end

        
        
        % -------------------------------------------------------
        function err = LoadHdf5(obj, fileobj, location)
            err = 0;
            
            % Arg 1
            if ~exist('fileobj','var') || (ischar(fileobj) && ~exist(fileobj,'file'))
                fileobj = '';
            end
                        
            % Arg 2
            if ~exist('location', 'var') || isempty(location)
                location = '/nirs/probe';
            elseif location(1)~='/'
                location = ['/',location];
            end
              
            % Error checking            
            if ~isempty(fileobj) && ischar(fileobj)
                obj.filename = fileobj;
            elseif isempty(fileobj)
                fileobj = obj.filename;
            end 
            if isempty(fileobj)
               err = -1;
               return;
            end
            
            try
                % Open group
                [gid, fid] = HDF5_GroupOpen(fileobj, location);
                
                % Load datasets
                obj.wavelengths               = HDF5_DatasetLoad(gid, 'wavelengths');
                obj.wavelengthsEmission       = HDF5_DatasetLoad(gid, 'wavelengthsEmission');
                obj.sourcePos2D               = HDF5_DatasetLoad(gid, 'sourcePos2D', [], '2D');
                obj.detectorPos2D             = HDF5_DatasetLoad(gid, 'detectorPos2D', [], '2D');
                obj.landmarkPos2D             = HDF5_DatasetLoad(gid, 'landmarkPos2D', [], '2D');
                obj.sourcePos3D               = HDF5_DatasetLoad(gid, 'sourcePos3D', [], '2D');
                obj.detectorPos3D             = HDF5_DatasetLoad(gid, 'detectorPos3D', [], '2D');
                obj.frequencies               = HDF5_DatasetLoad(gid, 'frequencies');
                obj.timeDelays                 = HDF5_DatasetLoad(gid, 'timeDelays');
                obj.timeDelayWidths            = HDF5_DatasetLoad(gid, 'timeDelayWidths');
                obj.momentOrder               = HDF5_DatasetLoad(gid, 'momentOrder');
                obj.correlationTimeDelays      = HDF5_DatasetLoad(gid, 'correlationTimeDelays');
                obj.correlationTimeDelayWidths = HDF5_DatasetLoad(gid, 'correlationTimeDelayWidths');
                obj.sourceLabels              = HDF5_DatasetLoad(gid, 'sourceLabels', obj.sourceLabels);
                obj.detectorLabels            = HDF5_DatasetLoad(gid, 'detectorLabels', obj.detectorLabels);
                obj.landmarkLabels            = HDF5_DatasetLoad(gid, 'landmarkLabels', obj.landmarkLabels);
                                
                % Close group
                HDF5_GroupClose(fileobj, gid, fid);
                
                assert(obj.IsValid())
                
            catch 
                err=-1;
                return;
            end
            
            % Call method to change future current and future versions of
            % SNIRF data to Homer3 compatible structure
            obj.ForwardCompatibility();
            
        end

        
        
        % -------------------------------------------------------
        function SaveHdf5(obj, fileobj, location)
            % Arg 1
            if ~exist('fileobj', 'var') || isempty(fileobj)
                error('Unable to save file. No file name given.')
            end
            
            % Arg 2
            if ~exist('location', 'var') || isempty(location)
                location = '/nirs/probe';
            elseif location(1)~='/'
                location = ['/',location];
            end
            
            if ~exist(fileobj, 'file')
                fid = H5F.create(fileobj, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
                H5F.close(fid);
            end     
            hdf5write_safe(fileobj, [location, '/wavelengths'], obj.wavelengths);
            hdf5write_safe(fileobj, [location, '/wavelengthsEmission'], obj.wavelengthsEmission);
            if ~isempty(obj.sourcePos2D)
                hdf5write_safe(fileobj, [location, '/sourcePos2D'], obj.sourcePos2D(:,1:2), 'rw:2D');
                hdf5write_safe(fileobj, [location, '/detectorPos2D'], obj.detectorPos2D(:,1:2), 'rw:2D');
            end
            if ~isempty(obj.sourcePos3D)
                hdf5write_safe(fileobj, [location, '/sourcePos3D'], obj.sourcePos3D(:,1:3), 'rw:2D');
                hdf5write_safe(fileobj, [location, '/detectorPos3D'], obj.detectorPos3D(:,1:3), 'rw:2D');
            end
            hdf5write_safe(fileobj, [location, '/frequencies'], obj.frequencies);
            hdf5write_safe(fileobj, [location, '/timeDelays'], obj.timeDelays);
            hdf5write_safe(fileobj, [location, '/timeDelayWidths'], obj.timeDelayWidths);
            hdf5write_safe(fileobj, [location, '/momentOrder'], obj.momentOrder);
            hdf5write_safe(fileobj, [location, '/correlationTimeDelays'], obj.correlationTimeDelays);
            hdf5write_safe(fileobj, [location, '/correlationTimeDelayWidths'], obj.correlationTimeDelayWidths);
            hdf5write_safe(fileobj, [location, '/sourceLabels'], obj.sourceLabels);
            hdf5write_safe(fileobj, [location, '/detectorLabels'], obj.detectorLabels);
        end
        
        
        
        % ---------------------------------------------------------
        function wls = GetWls(obj)
            wls = obj.wavelengths;
        end
        
        
        
        % ---------------------------------------------------------
        function srcpos = GetSrcPos(obj)
            srcpos = obj.sourcePos2D;
        end
        
        
        % ---------------------------------------------------------
        function detpos = GetDetPos(obj)
            detpos = obj.detectorPos2D;
        end
        
        
        % -------------------------------------------------------
        function B = eq(obj, obj2)
            B = false;
            if ~all(obj.wavelengths(:)==obj2.wavelengths(:))
                return;
            end
            if ~all(obj.wavelengthsEmission(:)==obj2.wavelengthsEmission(:))
                return;
            end
            if ~all(obj.sourcePos2D(:)==obj2.sourcePos2D(:))
                return;
            end
            if ~all(obj.detectorPos2D(:)==obj2.detectorPos2D(:))
                return;
            end
            if ~all(obj.frequencies(:)==obj2.frequencies(:))
                return;
            end
            if ~all(obj.timeDelays(:)==obj2.timeDelays(:))
                return;
            end
            if ~all(obj.timeDelayWidths(:)==obj2.timeDelayWidths(:))
                return;
            end
            if ~all(obj.momentOrder(:)==obj2.momentOrder(:))
                return;
            end
            if ~all(obj.correlationTimeDelays(:)==obj2.correlationTimeDelays(:))
                return;
            end
            if ~all(obj.correlationTimeDelayWidths(:)==obj2.correlationTimeDelayWidths(:))
                return;
            end
            if length(obj.sourceLabels)~=length(obj2.sourceLabels)
                return;
            end
            for ii=1:length(obj.sourceLabels)
                if ~strcmp(obj.sourceLabels{ii}, obj2.sourceLabels{ii})
                    return;
                end
            end
            if length(obj.detectorLabels)~=length(obj2.detectorLabels)
                return;
            end
            for ii=1:length(obj.detectorLabels)
                if ~strcmp(obj.detectorLabels{ii}, obj2.detectorLabels{ii})
                    return;
                end
            end
            B = true;
        end
        
        
        % ----------------------------------------------------------------------------------
        function nbytes = MemoryRequired(obj)
            nbytes = 0;
            fields = properties(obj);
            for ii=1:length(fields)
                nbytes = nbytes + eval(sprintf('sizeof(obj.%s)', fields{ii}));
            end
        end
        
        
        % ----------------------------------------------------------------------------------
        function b = IsEmpty(obj)
            b = true;
            if isempty(obj.wavelengths)
                return;
            end
            if isempty(obj.sourcePos2D) && isempty(obj.detectorPos2D) && ...
                    isempty(obj.sourcePos3D) && isempty(obj.detectorPos3D) 
                return;
            end
            b = false;
        end

        
        % ----------------------------------------------------------------------------------
        function b = IsValid(obj)
            b = false;
            if obj.IsEmpty()
                return;
            end
            if iscolumn(obj.sourcePos2D)
                return;
            end
            if length(obj.sourcePos2D)>4
                if size(obj.sourcePos2D,2) > size(obj.sourcePos2D,1)
                    return;
                end
            end
            if iscolumn(obj.detectorPos2D)
                return;
            end
            if length(obj.detectorPos2D)>4
                if size(obj.detectorPos2D,2) > size(obj.detectorPos2D,1)
                    return;
                end
            end
            b = true;
        end
        
    end
    
end
