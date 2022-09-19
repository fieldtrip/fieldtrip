classdef DataClass < matlab.mixin.Copyable
    
    % SNIRF-spec class properties
    properties
        dataTimeSeries
        time
        measurementList
    end
    
    % Non-SNIRF class properties
    properties (Access = private)
        filename
        fileformat
    end
    
    methods
        
        % -------------------------------------------------------
        function obj = DataClass(varargin)
            %
            %  Syntax:
            %     obj = DataClass()
            %     obj = DataClass(filename)
            %     obj = DataClass(data)
            %     obj = DataClass(d,t,ml)
            %
            %  Input:
            %     filename - When there's one argument and it is a char,
            %                then it's interepreted as a filename path
            %     data - When there's one argument and it it is not a char string, 
            %            it can be either a DataClass or NirsClass object.
            %     d    - When there are three arguments, d is the data time course matrix
            %     t    - When there are three arguments, t is the data time vector
            %     ml   - When there are three arguments, ml is the measurent list, which
            %            can be either a nirs style matrix or a MeasListClass object
            %
            %  Examples:
            %     
            %     % Example 1 - Create DataClass object and initialize it with SNIRF data variable 
            %     %             from file neuro_run01.snirf
            %
            %     data = DataClass('c:/users/public/subjects/subj1/neuro_run01.snirf')   
            %
            %    
            %     % Example 2 - Create DataClass object and initialize it with time course data and time vectors 
            %     %             from the .nirs file ./s1/neuro_run01.nirs
            %
            %     nirs = NirsClass('./s1/neuro_run01.nirs')
            %     data = DataClass(nirs.d, nirs.t)
            % 
            obj.fileformat = 'hdf5';
            
            % Set SNIRF fomat properties
            obj.measurementList = MeasListClass().empty();
            
            if nargin==0
                return;
            elseif nargin==1
                if isa(varargin{1}, 'DataClass')
                    obj.Copy(varargin{1});
                elseif isa(varargin{1}, 'NirsClass')
                    obj.dataTimeSeries = varargin{1}.d;
                    obj.time = varargin{1}.t;
                    for ii=1:size(varargin{1}.ml,1)
                        obj.measurementList(end+1) = MeasListClass(varargin{1}.ml(ii,:));
                    end
                elseif isa(varargin{1}, 'char')
                    obj.filename = varargin{1};
                    obj.Load();
                end
            elseif nargin==3
                if ~all(isreal(varargin{1}(:)))
                    return;
                end
                if ~all(isreal(varargin{2}(:)))
                    return;
                end
                if ~isa(varargin{3}, 'MeasListClass') && ~all(iswholenum(varargin{3}(:)))
                    return;
                end
                if isa(varargin{3}, 'MeasListClass')
                    obj.dataTimeSeries = varargin{1};
                    obj.time = varargin{2};
                    for ii=1:length(varargin{3})
                        obj.measurementList(end+1) = MeasListClass(varargin{3}(ii));
                    end
                else
                    obj.dataTimeSeries = varargin{1};
                    obj.time = varargin{2};
                    for ii=1:size(varargin{3},1)
                        obj.measurementList(end+1) = MeasListClass(varargin{3}(ii,:));
                    end
                end
            else
                obj.dataTimeSeries = double([]);
                obj.time = double([]);
            end
        end
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load/Save from/to file methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        % -------------------------------------------------------
        function err = LoadHdf5(obj, fileobj, location)
            err = 0;
            
            % Arg 1
            if ~exist('fileobj','var') || (ischar(fileobj) && ~exist(fileobj,'file'))
                fileobj = '';
            end
                      
            % Arg 2
            if ~exist('location', 'var') || isempty(location)
                location = '/nirs/data1';
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
                
                obj.dataTimeSeries  = HDF5_DatasetLoad(gid, 'dataTimeSeries');
                obj.time            = HDF5_DatasetLoad(gid, 'time');
                
                ii=1;
                while 1
                    if ii > length(obj.measurementList)
                        obj.measurementList(ii) = MeasListClass;
                    end
                    if obj.measurementList(ii).LoadHdf5(fileobj, [location, '/measurementList', num2str(ii)]) < 0
                        obj.measurementList(ii).delete();
                        obj.measurementList(ii) = [];
                        if ii==1
                            err=-1;
                        end
                        break;
                    end
                    ii=ii+1;
                end
                
                % Close group
                HDF5_GroupClose(fileobj, gid, fid);
            catch ME
                err = -1;
            end
        end
        
        
        % -------------------------------------------------------
        function SaveHdf5(obj, fileobj, location)
            if ~exist('fileobj', 'var') || isempty(fileobj)
                error('Unable to save file. No file name given.')
            end
            
            % Arg 2
            if ~exist('location', 'var') || isempty(location)
                location = '/nirs/data1';
            elseif location(1)~='/'
                location = ['/',location];
            end
            
            if ~exist(fileobj, 'file')
                fid = H5F.create(fileobj, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
                H5F.close(fid);
            end
            
            hdf5write_safe(fileobj, [location, '/dataTimeSeries'], obj.dataTimeSeries);
            hdf5write_safe(fileobj, [location, '/time'], obj.time);
            
            for ii=1:length(obj.measurementList)
                obj.measurementList(ii).SaveHdf5(fileobj, [location, '/measurementList', num2str(ii)]);
            end
        end
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set/Get properties methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        % -------------------------------------------------------
        function b = IsEmpty(obj)
            b = false;
            if isempty(obj)
                b = true;
            end
            if isempty(obj.dataTimeSeries) || isempty(obj.time) || isempty(obj.measurementList)
                b = true;
            end
            if isempty(obj.measurementList)
                return;
            end
        end
        
        
        % ---------------------------------------------------------
        function val = GetT(obj)
            val = obj.time;
        end
        
        
        % ---------------------------------------------------------
        function wls = GetWls(obj)
            wls = obj.measurementList(1).GetWls();
        end
        
        
        % ---------------------------------------------------------
        function ml = GetMeasList(obj)
            % Preallocate for speed 
            ml = ones(length(obj.measurementList), 4);
            
            % Convert obj.measurementList to matrix
            for ii = 1:length(obj.measurementList)
                % If this data contains block average then only get the measurements for first condition. That will
                % contain all the measurement channels
                if obj.measurementList(ii).GetCondition()>1
                    break;
                end
                ml(ii,:) = [obj.measurementList(ii).GetSourceIndex(), obj.measurementList(ii).GetDetectorIndex(), 1, obj.measurementList(ii).GetWavelengthIndex()];
            end
            
            % Remove unused rows that were pre-allocated
            ml(ii+1:end,:) = [];

            % Sort according to wavelength
            ml = sortrows(ml,4);
        end
        
        
        % ---------------------------------------------------------
        function ml = GetMeasListSrcDetPairs(obj)
            ml = zeros(0, 2);
            jj=1;
            for ii=1:length(obj.measurementList)
                if isempty(find(ml(:,1)==obj.measurementList(ii).GetSourceIndex() & ml(:,2)==obj.measurementList(ii).GetDetectorIndex()))
                    ml(jj,:) = [obj.measurementList(ii).GetSourceIndex(), obj.measurementList(ii).GetDetectorIndex()];
                    jj=jj+1;
                end
            end
        end
        
        
        % ---------------------------------------------------------
        function idxs = GetMeasurementListIdxs(obj, CondIdxs)
            % Get all the measurementList array idxs matching the
            % conditions in the CondNames argument
            idxs = zeros(1,length(obj.measurementList));
            kk=1;
            for iCh = 1:length(obj.measurementList)
                for iCond = 1:length(CondIdxs)
                    if sum(obj.measurementList(iCh).dataTypeIndex == CondIdxs)
                        idxs(kk) = iCh;
                        kk=kk+1;
                        break;
                    end
                end
            end
            idxs(idxs==0) = [];
        end
        
        
        % ---------------------------------------------------------
        function t = GetTime(obj)
            t = obj.time;
        end
        
        
        % ---------------------------------------------------------
        function d = GetDataTimeSeries(obj, options)
            d = [];
            if ~exist('options','var') || isempty(options)
                options = '';
            end
            if isempty(obj.dataTimeSeries)
                return;
            end
            if ~strcmp(options, 'reshape')
                d = obj.dataTimeSeries;
                return
            end
            
            % Get information for each ch in d matrix
            dataTypeLabels = {};
            srcDetPairs = zeros(0,2);
            conditions = [];
            wavelengths = [];
            hh=1; jj=1; kk=1; ll=1;
            for ii=1:length(obj.measurementList)
                if ~ismember(obj.measurementList(ii).GetDataTypeLabel(), dataTypeLabels)
                    dataTypeLabels{hh} = obj.measurementList(ii).GetDataTypeLabel();
                    hh=hh+1;
                end
                if isempty(find(srcDetPairs(:,1)==obj.measurementList(ii).GetSourceIndex() & srcDetPairs(:,2)==obj.measurementList(ii).GetDetectorIndex()))
                    srcDetPairs(jj,:) = [obj.measurementList(ii).GetSourceIndex(), obj.measurementList(ii).GetDetectorIndex()];
                    jj=jj+1;
                end
                if ~ismember(obj.measurementList(ii).GetCondition(), conditions)
                    conditions(kk) = obj.measurementList(ii).GetCondition();
                    kk=kk+1;
                end
                if ~ismember(obj.measurementList(ii).GetWavelengthIndex(), wavelengths)
                    wavelengths(ll) = obj.measurementList(ii).GetWavelengthIndex();
                    ll=ll+1;
                end
            end
            dim1 = length(obj.dataTimeSeries(:,1));
            if all(wavelengths(:)~=0) && all(conditions(:)==0)
                dim2 = length(wavelengths(:)) * size(srcDetPairs,1);
                dim3 = 1;
                dim4 = 1;
            elseif all(wavelengths(:)~=0) && all(conditions(:)~=0)
                dim2 = length(wavelengths(:)) * size(srcDetPairs,1);
                dim3 = length(conditions(:));
                dim4 = 1;
            elseif all(wavelengths(:)==0) && all(conditions(:)==0)
                dim2 = length(dataTypeLabels);
                dim3 = size(srcDetPairs,1);
                dim4 = 1;
            elseif all(wavelengths(:)==0) && all(conditions(:)~=0)
                dim2 = length(dataTypeLabels);
                dim3 = size(srcDetPairs,1);
                dim4 = length(conditions(:));
            end
            d = reshape(obj.dataTimeSeries, dim1, dim2, dim3, dim4);
        end
        
        
        % ---------------------------------------------------------
        function SetDataTimeSeries(obj, val)
            if ~exist('val','var')
                return;
            end
            obj.dataTimeSeries = val;
        end
        
        
        % ---------------------------------------------------------
        function SetTime(obj, val, datacheck)
            if ~exist('val','var')
                return;
            end
            if ~exist('datacheck','var')
                datacheck = false;
            end
            if isempty(obj.dataTimeSeries) && datacheck==true
                obj.time = [];
                return;
            end
            obj.time = val;
        end
        
        
        % ---------------------------------------------------------
        function SetDataType(obj, dataType, dataTypeLabel, chIdxs)
            if ~exist('dataType','var') ||  isempty(dataType)
                return;
            end
            if ~exist('dataTypeLabel','var') ||  isempty(dataTypeLabel)
                dataTypeLabel = '';
            end
            if ~exist('chIdxs','var') || isempty(chIdxs)
                chIdxs = 1:length(obj.measurementList);
            end
            for ii=chIdxs
                obj.measurementList(ii).SetDataType(dataType, dataTypeLabel);
            end
        end
        
        
        % ---------------------------------------------------------
        function SetDataTypeDod(obj)
            vals = DataTypeValues();
            for ii=1:length(obj.measurementList)
                obj.measurementList(ii).SetDataType(vals.Processed, 'dOD');
            end
        end
        
        
        
        % ---------------------------------------------------------
        function SetMl(obj, val)
            obj.measurementList = val.copy();      % shallow copy ok because MeasListClass has no handle properties
        end
        
        
        % ---------------------------------------------------------
        function val = GetMl(obj)
            val = obj.measurementList;
        end
        
        
        % ---------------------------------------------------------
        function val = GetDataType(obj)
            val = zeros(length(obj.measurementList),1);
            for ii=1:length(obj.measurementList)
                val(ii) = obj.measurementList(ii).GetDataType();
            end
        end
        
        
        % ---------------------------------------------------------
        function val = GetDataTypeLabel(obj, ch_idx)
            if ~exist('ch_idx','var')
                ch_idx = 1:length(obj.measurementList);
            end
            val = repmat({''}, length(ch_idx),1);
            for ii=ch_idx
                val{ii} = obj.measurementList(ii).GetDataTypeLabel();
            end
            val = unique(val);
        end
        
        
        % ---------------------------------------------------------
        function val = GetCondition(obj, ch_idx)
            if ~exist('ch_idx','var')
                ch_idx = 1:length(obj.measurementList);
            end
            val = zeros(length(ch_idx),1);
            for ii=ch_idx
                val(ii) = obj.measurementList(ii).GetCondition();
            end
        end
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Adding/deleting data methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        % ---------------------------------------------------------
        function AddChannelHb(obj, isrc, idet, iHb, icond)
            if ~exist('isrc','var') || isempty(isrc)
                return;
            end
            if ~exist('idet','var') || isempty(idet)
                return;
            end
            if ~exist('iHb','var') || isempty(iHb)
                iHb = 1;
            end
            if ~exist('icond','var') || isempty(icond)
                icond = 0;
            end
            switch(iHb)
                case 1
                    AddChannelHbO(obj, isrc, idet, icond);
                case 2
                    AddChannelHbR(obj, isrc, idet, icond);
                case 3
                    AddChannelHbT(obj, isrc, idet, icond);
            end
            
        end
        
        
        % ---------------------------------------------------------
        function AddChannelHbO(obj, isrc, idet, icond)
            if ~exist('isrc','var') || isempty(isrc)
                return;
            end
            if ~exist('idet','var') || isempty(idet)
                return;
            end
            if ~exist('icond','var') || isempty(icond)
                icond = 0;
            end
            vals = DataTypeValues();
            if icond==0
                obj.measurementList(end+1) = MeasListClass(isrc, idet, vals.Processed, 'HbO');
            else
                obj.measurementList(end+1) = MeasListClass(isrc, idet, vals.Processed, 'HRF HbO', icond);
            end
        end
        
        
        % ---------------------------------------------------------
        function AddChannelHbR(obj, isrc, idet, icond)
            if ~exist('isrc','var') || isempty(isrc)
                return;
            end
            if ~exist('idet','var') || isempty(idet)
                return;
            end
            if ~exist('icond','var') || isempty(icond)
                icond = 0;
            end
            vals = DataTypeValues();
            if icond==0
                obj.measurementList(end+1) = MeasListClass(isrc, idet, vals.Processed, 'HbR');
            else
                obj.measurementList(end+1) = MeasListClass(isrc, idet, vals.Processed, 'HRF HbR', icond);
            end
        end
        
        
        % ---------------------------------------------------------
        function AddChannelHbT(obj, isrc, idet, icond)
            if ~exist('isrc','var') || isempty(isrc)
                return;
            end
            if ~exist('idet','var') || isempty(idet)
                return;
            end
            if ~exist('icond','var') || isempty(icond)
                icond = 0;
            end
            vals = DataTypeValues();
            if icond==0
                obj.measurementList(end+1) = MeasListClass(isrc, idet, vals.Processed, 'HbT');
            else
                obj.measurementList(end+1) = MeasListClass(isrc, idet, vals.Processed, 'HRF HbT', icond);
            end
        end
        
        
        % ---------------------------------------------------------
        function AddChannelDod(obj, isrc, idet, wl, icond)
            if ~exist('isrc','var') || isempty(isrc)
                return;
            end
            if ~exist('idet','var') || isempty(idet)
                return;
            end
            if ~exist('icond','var') || isempty(icond)
                icond = 0;
            end
            vals = DataTypeValues();
            if icond==0
                obj.measurementList(end+1) = MeasListClass(isrc, idet, vals.Processed, 'dOD');
            else
                obj.measurementList(end+1) = MeasListClass(isrc, idet, vals.Processed, 'HRF dOD', icond);
            end
            obj.measurementList(end).SetWavelengthIndex(wl);
        end
        
        
        % ---------------------------------------------------------
        function AppendDataTimeSeries(obj, y)
            obj.dataTimeSeries(:, end+1:end+size(y(:,:),2)) = y(:,:);
        end
        
        
        % ---------------------------------------------------------
        function TruncateTpts(obj, n)
            obj.dataTimeSeries(end-n+1:end, :) = [];
            obj.time(end-n+1:end) = [];
        end
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        % ----------------------------------------------------------------------------------
        function Copy(obj, obj2)
            if isempty(obj)
                obj = DataClass();
            end
            if isempty(obj2)
                obj = DataClass();
                return;
            end
            if ~isa(obj2, 'DataClass')
                return;
            end
            for ii=1:length(obj2.measurementList)
                obj.measurementList(ii) = obj2.measurementList(ii).copy();      % shallow copy ok because MeasListClass has no handle properties
            end
            obj.dataTimeSeries = obj2.dataTimeSeries;
            obj.time = obj2.time;
        end
        
        
        % -------------------------------------------------------
        function B = eq(obj, obj2)
            B = false;
            if length(obj.dataTimeSeries(:)) ~= length(obj2.dataTimeSeries(:))
                return;
            end
            if ndims(obj.dataTimeSeries) ~= ndims(obj2.dataTimeSeries)
                return;
            end
            if ~all(size(obj.dataTimeSeries)==size(obj2.dataTimeSeries))
                return;
            end
            if ~all(obj.dataTimeSeries(:)==obj2.dataTimeSeries(:))
                return;
            end
            if ~all(obj.time(:)==obj2.time(:))
                return;
            end
            if length(obj.measurementList)~=length(obj2.measurementList)
                return;
            end
            for ii=1:length(obj.measurementList)
                if ~(obj.measurementList(ii) == obj2.measurementList(ii))
                    return;
                end
            end
            B = true;
        end
        
        
        % ----------------------------------------------------------------------------------
        function nbytes = MemoryRequired(obj)
            nbytes = 0;
            if isempty(obj)
                return
            end
            nbytes = sizeof(obj.dataTimeSeries) + sizeof(obj.time);
            for ii=1:length(obj.measurementList)
                nbytes = nbytes + obj.measurementList(ii).MemoryRequired();
            end
        end
        
        
    end
end

