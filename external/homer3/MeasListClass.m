classdef MeasListClass < matlab.mixin.Copyable
    
    % SNIRF-spec class properties
    properties
        sourceIndex
        detectorIndex
        wavelengthIndex
        dataType
        dataTypeLabel
        dataTypeIndex   % Used for condition when dataType=99999 ("Processed") and dataTypeLabel='HRF...'
        sourcePower
        detectorGain
        moduleIndex
    end
    
    % Non-SNIRF class properties
    properties
        filename
        fileformat
    end
    
    
    methods

        function obj = MeasListClass(varargin)
            %
            %  Syntax:
            %     obj = MeasListClass()
            %     obj = MeasListClass(ml)
            %     obj = MeasListClass(sourceIndex, detectorIndex, wavelengthIndex)
            %     obj = MeasListClass(sourceIndex, detectorIndex, dataType)
            %     obj = MeasListClass(sourceIndex, detectorIndex, dataType, dataTypeLabel)
            %     obj = MeasListClass(sourceIndex, detectorIndex, dataType, dataTypeLabel, condition)
            %     
            %  Inputs:
            %     ml             - When there's one argument, ml is the measurent list, which 
            %                      can be either a nirs style matrix or a MeasListClass object.
            %     sourceIndex    - When there are more than 2 arguments, ...
            %     detectorIndex  - When there are more than 2 arguments, ...
            %     dataType       - When there are more than 2 arguments, ...
            %     dataTypeLabel  - When there are more than 2 arguments, ...
            %     dataTypeIndex  - When there are more than 2 arguments, ...
            %
            %  Example:
            %
            
            % Fields which are part of the SNIRF spec which are loaded and saved 
            % from/to SNIRF files
            obj.sourceIndex      = 0;
            obj.detectorIndex    = 0;
            obj.wavelengthIndex  = 0;
            obj.dataType         = 0;
            obj.dataTypeLabel    = '';
            obj.dataTypeIndex    = 0;
            obj.sourcePower      = 0;
            obj.detectorGain     = 0;
            obj.moduleIndex      = 0;
            
            dataTypeValues = DataTypeValues();

            if nargin==1 && isa(varargin{1}, 'MeasListClass')
                obj                  = varargin{1}.copy();                    % shallow copy ok because MeasListClass has no handle properties 
            elseif nargin==1 
                obj.sourceIndex      = varargin{1}(1);
                obj.detectorIndex    = varargin{1}(2);
                obj.wavelengthIndex  = varargin{1}(4);
                obj.dataType         = dataTypeValues.Raw.CW.Amplitude;
            elseif nargin==3
                obj.sourceIndex      = varargin{1};
                obj.detectorIndex    = varargin{2};
                obj.dataType         = varargin{3};
            elseif nargin==4
                obj.sourceIndex      = varargin{1};
                obj.detectorIndex    = varargin{2};
                obj.dataType         = varargin{3};
                obj.dataTypeLabel    = varargin{4};
            elseif nargin==5
                obj.sourceIndex      = varargin{1};
                obj.detectorIndex    = varargin{2};
                obj.dataType         = varargin{3};
                obj.dataTypeLabel    = varargin{4};
                obj.dataTypeIndex    = varargin{5};
            end
            
            % Set base class properties not part of the SNIRF format
            obj.fileformat = 'hdf5';

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
                location = '/nirs/data1/measurementList1';
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
                obj.sourceIndex     = HDF5_DatasetLoad(gid, 'sourceIndex');
                obj.detectorIndex   = HDF5_DatasetLoad(gid, 'detectorIndex');
                obj.wavelengthIndex = HDF5_DatasetLoad(gid, 'wavelengthIndex');
                obj.dataType        = HDF5_DatasetLoad(gid, 'dataType');
                obj.dataTypeLabel   = HDF5_DatasetLoad(gid, 'dataTypeLabel', obj.dataTypeLabel);
                obj.detectorIndex   = HDF5_DatasetLoad(gid, 'detectorIndex');
                obj.sourcePower     = HDF5_DatasetLoad(gid, 'sourcePower');
                obj.sourcePower     = HDF5_DatasetLoad(gid, 'sourcePower');
                obj.moduleIndex     = HDF5_DatasetLoad(gid, 'moduleIndex');
                
                HDF5_GroupClose(fileobj, gid, fid);
            catch ME
                err = -1;
                return
            end

        end

        
        % -------------------------------------------------------
        function SaveHdf5(obj, fileobj, location)

            % Arg 1
            if ~exist('fileobj', 'var') || isempty(fileobj)
                error('Unable to save file. No file name given.')
            end
            
            % Arg 2
            if ~exist('location', 'var') || isempty(location)
                location = '/nirs/data1/measurementList1';
            elseif location(1)~='/'
                location = ['/',location];
            end
            
            if ~exist(fileobj, 'file')
                fid = H5F.create(fileobj, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
                H5F.close(fid);
            end
            
            hdf5write_safe(fileobj, [location, '/sourceIndex'], obj.sourceIndex);
            hdf5write_safe(fileobj, [location, '/detectorIndex'], obj.detectorIndex);
            hdf5write_safe(fileobj, [location, '/wavelengthIndex'], obj.wavelengthIndex);
            hdf5write_safe(fileobj, [location, '/dataType'], obj.dataType);
            hdf5write_safe(fileobj, [location, '/dataTypeLabel'], obj.dataTypeLabel);
            hdf5write_safe(fileobj, [location, '/dataTypeIndex'], obj.dataTypeIndex);
            hdf5write_safe(fileobj, [location, '/sourcePower'], obj.sourcePower);
            hdf5write_safe(fileobj, [location, '/detectorGain'], obj.detectorGain);
            hdf5write_safe(fileobj, [location, '/moduleIndex'], obj.moduleIndex);
        end

                
        % ---------------------------------------------------------
        function idx = GetSourceIndex(obj)
            idx = obj.sourceIndex;
        end
        
        
        % ---------------------------------------------------------
        function idx = GetDetectorIndex(obj)
            idx = obj.detectorIndex;
        end
        
        
        % ---------------------------------------------------------
        function idx = GetWavelengthIndex(obj)
            idx = obj.wavelengthIndex;
        end
        
        
        % ---------------------------------------------------------
        function SetWavelengthIndex(obj, val)
            obj.wavelengthIndex = val;
        end
        
        
        % ---------------------------------------------------------
        function SetDataType(obj, dataType, dataTypeLabel)
            obj.dataType = dataType;
            obj.dataTypeLabel = dataTypeLabel;
        end
        
        
        % ---------------------------------------------------------
        function SetDataTypeLabel(obj, dataTypeLabel)
            obj.dataTypeLabel = dataTypeLabel;
        end
        
        
        % ---------------------------------------------------------
        function [dataType, dataTypeLabel] = GetDataType(obj)
            dataType = obj.dataType;
            dataTypeLabel = obj.dataTypeLabel;
        end
        
        
        % ---------------------------------------------------------
        function dataTypeLabel = GetDataTypeLabel(obj)
            dataTypeLabel = obj.dataTypeLabel;
        end
        
        
        % ---------------------------------------------------------
        function SetCondition(obj, val)
            obj.dataTypeIndex = val;
        end
        
        
        % ---------------------------------------------------------
        function val = GetCondition(obj)
            val = obj.dataTypeIndex;
        end
        
        
        % -------------------------------------------------------
        function b = IsEmpty(obj)
            b = false;
            if obj.sourceIndex==0 && obj.detectorIndex==0
                b = true;
            end
        end
        
        % -------------------------------------------------------
        function B = eq(obj, obj2)
            B = false;       
            if obj.sourceIndex~=obj2.sourceIndex
                return;
            end
            if obj.detectorIndex~=obj2.detectorIndex
                return;
            end
            if obj.wavelengthIndex~=obj2.wavelengthIndex
                return;
            end
            if obj.dataType~=obj2.dataType
                return;
            end
            if ~strcmp(obj.dataTypeLabel, obj2.dataTypeLabel)
                return;
            end
            if obj.dataTypeIndex~=obj2.dataTypeIndex
                return;
            end
            if obj.sourcePower~=obj2.sourcePower
                return;
            end
            if obj.detectorGain~=obj2.detectorGain
                return;
            end
            if obj.moduleIndex~=obj2.moduleIndex
                return;
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
        
    end
    
end

