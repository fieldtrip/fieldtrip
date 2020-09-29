classdef SnirfClass < matlab.mixin.Copyable
    
    properties
        formatVersion
        metaDataTags
        data
        stim
        probe
        aux
    end
    
    properties (Access = public)
        fid
        gid
        location
        nirsdatanum
        nirs_tb
        stim0
        filename        
        fileformat
        options
    end
    
    methods
        
        % -------------------------------------------------------
        function obj = SnirfClass(varargin)
            %
            % Syntax:
            %   obj = SnirfClass()
            %   obj = SnirfClass(filename);
            %   obj = SnirfClass(filename, nirsdatanum);
            %   obj = SnirfClass(filename, nirsdatanum, options);
            %   obj = SnirfClass(filename, options);
            %   obj = SnirfClass(dotnirs);
            %   obj = SnirfClass(dotnirs, numdatabllocks);
            %   obj = SnirfClass(data, stim);
            %   obj = SnirfClass(data, stim, probe);
            %   obj = SnirfClass(data, stim, probe, aux);
            %   obj = SnirfClass(d, t, SD, aux, s);
            %   obj = SnirfClass(d, t, SD, aux, s, CondNames);
            %
            %   Also for debugging/simulation of time bases
            %
            %   obj = SnirfClass(dotnirs, tfactors);
            %
            % Example 1:
            %
            %   % Save .nirs file in SNIRF format
            %   snirf1 = SnirfClass(load('neuro_run01.nirs','-mat'));
            %   snirf1.Save('neuro_run01.snirf');
            %   snirf1.Info()
            %
            %   % Check that the file was saved correctly
            %   snirf2 = SnirfClass();
            %   snirf2.Load('neuro_run01.snirf');
            %   snirf2.Info()
            %
            % Example 2:
            %
            %   Nirs2Snirf('Simple_Probe1.nirs');
            %   obj = SnirfClass('Simple_Probe1.snirf');
            %
            %   Here's some of the output:
            %
            %   obj(1).data ====>
            %
            %       DataClass with properties:
            %
            %           data: [1200x8 double]
            %           time: [1200x1 double]
            %           measurementList: [1x8 MeasListClass]
            %
            
            % Initialize properties from SNIRF spec
            obj.Initialize()
            
            % Set class properties NOT part of the SNIRF format
            obj.fileformat = 'hdf5';
            obj.location = '/nirs';
            obj.nirsdatanum = 1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Between 1 and 4 arguments covers the following syntax variants
            %
            % obj = SnirfClass(filename);
            % obj = SnirfClass(filename, nirsdatanum);
            % obj = SnirfClass(filename, nirsdatanum, options);
            % obj = SnirfClass(filename, options);
            % obj = SnirfClass(dotnirs);
            % obj = SnirfClass(dotnirs, numdatabllocks);
            % obj = SnirfClass(data, stim);
            % obj = SnirfClass(data, stim, probe);
            % obj = SnirfClass(data, stim, probe, aux);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin>0 && nargin<5
                
                % obj = SnirfClass(filename);
                if isa(varargin{1}, 'SnirfClass')
                    obj.Copy(varargin{1});
                    return;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % obj = SnirfClass(filename, nirsdatanum);
                % obj = SnirfClass(filename, nirsdatanum, options);
                % obj = SnirfClass(filename, options);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ischar(varargin{1})
                    obj.filename = varargin{1};
                    if nargin>1
                        % obj = SnirfClass(filename, nirsdatanum);
                        if isnumeric(varargin{2})
                            obj.nirsdatanum = varargin{2};
                            
                            % obj = SnirfClass(filename, nirsdatanum, options);
                            if nargin>2
                                obj.options = varargin{3};
                            end
                            
                            % obj = SnirfClass(filename, options);
                        elseif ischar(varargin{2})
                            obj.options = varargin{2};
                            
                        end
                    end
                    
                    % Conditional loading of snirf file data
                    if strcmpi(obj.options, 'memory')
                        obj.Load(varargin{1});
                    end
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % obj = SnirfClass(dotnirs);
                    % obj = SnirfClass(dotnirs, numdatabllocks);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                elseif isstruct(varargin{1}) || isa(varargin{1}, 'NirsClass')
                    
                    % obj = SnirfClass(dotnirs);
                    tfactors = 1;    % Debug simulation parameter
                    
                    % obj = SnirfClass(dotnirs, numdatabllocks);
                    if nargin==2
                        tfactors = varargin{2};
                    end
                    dotnirs = varargin{1};
                    obj.GenSimulatedTimeBases(dotnirs, tfactors);
                    for ii=1:length(tfactors)
                        obj.data(ii) = DataClass(obj.nirs_tb(ii).d, obj.nirs_tb(ii).t(:), obj.nirs_tb(ii).SD.MeasList);
                    end
                    
                    for ii=1:size(dotnirs.s,2)
                        if isfield(dotnirs, 'CondNames')
                            obj.stim(ii) = StimClass(dotnirs.s(:,ii), dotnirs.t(:), dotnirs.CondNames{ii});
                        else
                            obj.stim(ii) = StimClass(dotnirs.s(:,ii), dotnirs.t(:), num2str(ii));
                        end
                    end
                    obj.probe      = ProbeClass(dotnirs.SD);
                    for ii=1:size(dotnirs.aux,2)
                        obj.aux(ii) = AuxClass(dotnirs.aux(:,ii), dotnirs.t(:), sprintf('aux%d',ii));
                    end
                    
                    % Add metadatatags
                    obj.metaDataTags   = MetaDataTagsClass();
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % obj = SnirfClass(data, stim);
                    % obj = SnirfClass(data, stim, probe);
                    % obj = SnirfClass(data, stim, probe, aux);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                elseif isa(varargin{1}, 'DataClass')
                    
                    % obj = SnirfClass(data, stim);
                    data = varargin{1};
                    obj.SetData(data);
                    stim = varargin{2};
                    obj.SetStim(stim);
                    
                    % obj = SnirfClass(data, stim, probe);
                    if nargin>2
                        probe = varargin{3};
                        obj.SetSd(probe);
                    end
                    
                    % obj = SnirfClass(data, stim, probe, aux);
                    if nargin>3
                        aux = varargin{4};
                        obj.SetAux(aux);
                    end
                    
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Between 5 and 6 arguments covers the following syntax variants
                %
                % obj = SnirfClass(d, t, SD, aux, s);
                % obj = SnirfClass(d, t, SD, aux, s, CondNames);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif nargin>4
                
                % obj = SnirfClass(d, t, SD, aux, s);
                d         = varargin{1};
                t         = varargin{2}(:);
                SD        = varargin{3};
                aux       = varargin{4};
                s         = varargin{5};
                CondNames = {};
                
                % obj = SnirfClass(d, t, SD, aux, s, CondNames);
                if nargin>5
                    CondNames = varargin{6};
                end
                
                obj.data(1) = DataClass(d, t(:), SD.MeasList);
                for ii=1:size(s,2)
                    if nargin==5
                        condition = num2str(ii);
                    else
                        condition = CondNames{ii};
                    end
                    obj.stim(ii) = StimClass(s(:,ii), t(:), condition);
                end
                obj.probe      = ProbeClass(SD);
                for ii=1:size(aux,2)
                    obj.aux(ii) = AuxClass(aux, t(:), sprintf('aux%d',ii));
                end
                
                % Add metadatatags
                obj.metaDataTags   = MetaDataTagsClass();
                
            end
            
        end
        
        
        
        % -------------------------------------------------------
        function Initialize(obj)
            obj.formatVersion = '1.0';
            obj.metaDataTags   = MetaDataTagsClass().empty();
            obj.data           = DataClass().empty();
            obj.stim           = StimClass().empty();
            obj.probe          = ProbeClass().empty();
            obj.aux            = AuxClass().empty();
            
            obj.stim0          = StimClass().empty();
            obj.options        = 'memory';
        end
        
        
        % -------------------------------------------------------
        function err = Copy(obj, obj2)
            err=0;
            if ~isa(obj2, 'SnirfClass')
                err=1;
                return;
            end
            obj.formatVersion = obj2.formatVersion;
            obj.metaDataTags  = CopyHandles(obj2.metaDataTags);
            obj.data          = CopyHandles(obj2.data);
            obj.stim          = CopyHandles(obj2.stim);
            obj.probe         = CopyHandles(obj2.probe);
            obj.aux           = CopyHandles(obj2.aux);
            
            try
                obj.stim0     = CopyHandles(obj2.stim0);
            catch
            end
        end
        
        
        
        % -------------------------------------------------------
        function objnew = CopyMutable(obj, options)
            if nargin==1
                options = '';
            end
            
            % If we're working off the snirf file instead of loading everything into memory
            % then we have to load stim here from file before accessing it.
            if strcmpi(obj.options, 'file')
                obj.LoadStim(obj.filename);
            end
            
            % Generate new instance of SnirfClass
            objnew = SnirfClass();
            
            objnew.filename = obj.filename;
            
            % Copy mutable properties to new object instance;
            objnew.stim = CopyHandles(obj.stim);
            
            if strcmp(options, 'extended')
                t = obj.GetTimeCombined();
                objnew.data = DataClass([],t,[]);
            end
        end
        
        
        
        % -------------------------------------------------------
        function SortStims(obj)
            if isempty(obj.stim)
                return;
            end
            temp = CopyHandles(obj.stim);
            delete(obj.stim);
            names = cell(length(temp),1);
            for ii=1:length(temp)
                names{ii} = temp(ii).name;
            end
            [~,idx] = sort(names);
            obj.stim = temp(idx).copy;
        end
        
        
        
        % -------------------------------------------------------
        function err = SetLocation(obj)
            err = 0;
            gid1 = HDF5_GroupOpen(obj.fid, sprintf('%s%d', obj.location, obj.nirsdatanum));
            gid2 = HDF5_GroupOpen(obj.fid, obj.location);
            
            if gid1.double > 0
                obj.location = sprintf('%s%d', obj.location, obj.nirsdatanum);
                return;
            elseif gid2.double > 0
                return;
            end
            err = -1;
        end
        
        
        
        
        % -------------------------------------------------------
        function err = LoadFormatVersion(obj)
            err = 0;
            formatVersionFile = HDF5_DatasetLoad(obj.gid, 'formatVersion'); %#ok<*PROPLC>
            formatVersionFile = str2double(formatVersionFile);
            formatVersionCurr = str2double(obj.formatVersion);
            if formatVersionFile < formatVersionCurr
                fprintf('Warning: Current SNIRF version is %0.1f. Cannot load older version (%0.1f) file. Backward compatibility not yet implemented ...\n', formatVersionCurr, formatVersionFile)
                err = -2;
                return
            end
        end
        
        
        
        % -------------------------------------------------------
        function err = LoadMetaDataTags(obj, fileobj)
            err = 0;
            obj.metaDataTags = MetaDataTagsClass();
            if obj.metaDataTags.LoadHdf5(fileobj, [obj.location, '/metaDataTags']) < 0
                err = -1;
            end
        end
        
        
        
        % -------------------------------------------------------
        function err = LoadData(obj, fileobj)
            err = 0;
            ii=1;
            while 1
                if ii > length(obj.data)
                    obj.data(ii) = DataClass;
                end
                if obj.data(ii).LoadHdf5(fileobj, [obj.location, '/data', num2str(ii)]) < 0
                    obj.data(ii).delete();
                    obj.data(ii) = [];
                    if ii==1
                        err = -1;
                    end
                    break;
                end
                ii=ii+1;
            end
        end
        
        
        % -------------------------------------------------------
        function err = LoadStim(obj, fileobj)
            err = 0;
            
            ii=1;
            while 1
                if ii > length(obj.stim)
                    obj.stim(ii) = StimClass;
                end
                if obj.stim(ii).LoadHdf5(fileobj, [obj.location, '/stim', num2str(ii)]) < 0
                    obj.stim(ii).delete();
                    obj.stim(ii) = [];
                    if ii==1
                        err = -1;
                    end
                    break;
                end
                ii=ii+1;
            end
            obj.SortStims();
            
            % Load original, unedited stims, if they exist
            ii=1;
            while 1
                if ii > length(obj.stim0)
                    obj.stim0(ii) = StimClass;
                end
                if obj.stim0(ii).LoadHdf5(fileobj, [obj.location, '/stim0', num2str(ii)]) < 0
                    obj.stim0(ii).delete();
                    obj.stim0(ii) = [];
                    break;
                end
                ii=ii+1;
            end
            
        end
        
        
        
        % -------------------------------------------------------
        function err = LoadProbe(obj, fileobj, ~)
            obj.probe = ProbeClass();
            err = obj.probe.LoadHdf5(fileobj, [obj.location, '/probe']);
        end
        
        
        
        % -------------------------------------------------------
        function err = LoadAux(obj, fileobj)
            err = 0;
            ii=1;
            while 1
                if ii > length(obj.aux)
                    obj.aux(ii) = AuxClass;
                end
                if obj.aux(ii).LoadHdf5(fileobj, [obj.location, '/aux', num2str(ii)]) < 0
                    obj.aux(ii).delete();
                    obj.aux(ii) = [];
                    if ii==1
                        err = -1;
                    end
                    break;
                end
                ii=ii+1;
            end
        end
        
        
        
        % -------------------------------------------------------
        function err = LoadHdf5(obj, fileobj, ~)
            err = 0;
            
            % Arg 1
            if ~exist('fileobj','var') || ~exist(fileobj,'file')
                fileobj = '';
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
            
            % Don't reload if not empty
            if ~obj.IsEmpty()
                return;
            end
            
            
            %%%%%%%%%%%% Ready to load from file
            
            try
                
                % Open group
                [obj.gid, obj.fid] = HDF5_GroupOpen(fileobj, '/');
                
                if obj.SetLocation() < 0
                    err = -1;
                    return
                end
                
                %%%% Load formatVersion
                if obj.LoadFormatVersion() < 0
                    err = -2;
                end
                
                %%%% Load metaDataTags
                if obj.LoadMetaDataTags(obj.fid) < 0
                    err = -3;
                end
                
                %%%% Load data
                if obj.LoadData(obj.fid) < 0
                    err = -4;
                end
                
                %%%% Load stim
                if obj.LoadStim(obj.fid)
                    err = -5;
                end
                
                %%%% Load probe
                if obj.LoadProbe(obj.fid)
                    err = -6;
                end
                
                %%%% Load aux
                if obj.LoadAux(obj.fid)
                    err = -7;
                end
                
                
                % Close group
                HDF5_GroupClose(fileobj, obj.gid, obj.fid);
                
            catch
                
                err = -1;
                
            end
            
            if obj.fid>0
                H5F.close(obj.fid);
            end
            
        end
        
        
        % -------------------------------------------------------
        function err = Load(obj, fileobj)
            err = LoadHdf5(obj, fileobj);
        end
        
        
        
        % -------------------------------------------------------
        function SaveMetaDataTags(obj, fileobj)
            obj.metaDataTags.SaveHdf5(fileobj, [obj.location, '/metaDataTags']);
        end
        
        
        
        % -------------------------------------------------------
        function SaveData(obj, fileobj)
            for ii=1:length(obj.data)
                obj.data(ii).SaveHdf5(fileobj, [obj.location, '/data', num2str(ii)]);
            end
        end
        
        
        % -------------------------------------------------------
        function SaveStim(obj, fileobj)
            for ii=1:length(obj.stim)
                obj.stim(ii).SaveHdf5(fileobj, [obj.location, '/stim', num2str(ii)]);
            end
            if isempty(obj.stim0)
                obj.stim0 = obj.stim.copy();
                for ii=1:length(obj.stim0)
                    obj.stim0(ii).SaveHdf5(fileobj, [obj.location, '/stim0', num2str(ii)]);
                end
            end
        end
        
        
        % -------------------------------------------------------
        function SaveProbe(obj, fileobj)
            obj.probe.SaveHdf5(fileobj, [obj.location, '/probe']);
        end
        
        
        % -------------------------------------------------------
        function SaveAux(obj, fileobj)
            for ii=1:length(obj.aux)
                obj.aux(ii).SaveHdf5(fileobj, [obj.location, '/aux', num2str(ii)]);
            end
        end
        
        
        % -------------------------------------------------------
        function SaveHdf5(obj, fileobj, ~)
            % Arg 1
            if ~exist('fileobj','var') || isempty(fileobj)
                error('Unable to save file. No file name given.')
            end
            
            % Args
            if exist(fileobj, 'file')
                delete(fileobj);
            end
            obj.fid = H5F.create(fileobj, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
            H5F.close(obj.fid);
            
            %%%%% Save this object's properties
            
            % Save formatVersion
            if isempty(obj.formatVersion)
                obj.formatVersion = '1.1';
            end
            hdf5write_safe(fileobj, '/formatVersion', obj.formatVersion);
            
            % Save metaDataTags
            obj.SaveMetaDataTags(fileobj);
            
            % Save data
            obj.SaveData(fileobj);
            
            % Save stim
            obj.SaveStim(fileobj);
            
            % Save sd
            obj.SaveProbe(fileobj);
            
            % Save aux
            obj.SaveAux(fileobj);
        end
        
        
        % -------------------------------------------------------
        function Save(obj, fileobj)
            SaveHdf5(obj, fileobj);
        end

        
                
        % -------------------------------------------------------
        function [stimFromFile, changes] = UpdateStim(obj, fileobj)
            flags = zeros(length(obj.stim), 1);
            
            % Load stim from file and update it
            snirfFile = SnirfClass();
            snirfFile.LoadStim(fileobj);
            
            % Update stims from file with edited stims
            for ii = 1:length(obj.stim)
                for jj = 1:length(snirfFile.stim)
                    if strcmp(obj.stim(ii).GetName(), snirfFile.stim(jj).GetName())
                        if obj.stim(ii) ~= snirfFile.stim(jj)
                            snirfFile.stim(jj).Copy(obj.stim(ii));
                        end
                        flags(ii) = 1;
                        break;
                    end
                end
                if ~flags(ii)
                    % We have new stimulus condition added
                    if ~obj.stim(ii).IsEmpty()
                        snirfFile.stim(jj+1) = StimClass(obj.stim(ii));
                    end
                end
            end
            
            % If stims were edited then update snirf file with new stims
            changes = sum(flags);
            if changes
                snirfFile.SaveStim(fileobj);
            end
            stimFromFile = snirfFile.stim;
        end
        
        
        
        % -------------------------------------------------------
        function changes = StimChangesMade(obj)
            
            flags = zeros(length(obj.stim), 1);
            
            % Load stims from file
            snirf = SnirfClass();
            snirf.LoadStim(obj.filename);
            stimFromFile = snirf.stim;
            
            % Update stims from file with edited stims
            for ii = 1:length(obj.stim)
                for jj = 1:length(stimFromFile)
                    if strcmp(obj.stim(ii).GetName(), stimFromFile(jj).GetName())
                        if obj.stim(ii) ~= stimFromFile(jj)
                            flags(ii) = 1;
                        else
                            flags(ii) = -1;
                        end
                        break;
                    end
                end
                if flags(ii)==0
                    % We have new stimulus condition added
                    if ~obj.stim(ii).IsEmpty()
                        flags(ii) = 1;
                    end
                end
            end
            flags(flags ~= 1) = 0;
            changes = sum(flags)>0;
        end
        
        
        
        % -------------------------------------------------------
        function b = DataModified(obj)
            b = obj.StimChangesMade();
        end
        
        
        
        % -------------------------------------------------------
        function err = SaveMutable(obj, fileobj)
            if isempty(obj)
                return
            end
            
            % Arg 1
            if ~exist('fileobj','var') || ~exist(fileobj,'file')
                fileobj = '';
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
            
            % Update original stims and save back to file
            obj.UpdateStim(fileobj);
        end
        
        
        
        
        % -------------------------------------------------------
        function B = eq(obj, obj2)
            B = false;
            if ~strcmp(obj.formatVersion, obj2.formatVersion)
                return;
            end
            if length(obj.data)~=length(obj2.data)
                return;
            end
            for ii=1:length(obj.data)
                if ~(obj.data(ii)==obj2.data(ii))
                    return;
                end
            end
            if length(obj.stim)~=length(obj2.stim)
                return;
            end
            for ii=1:length(obj.stim)
                flag = false;
                for jj=1:length(obj2.stim)
                    if obj.stim(ii)==obj2.stim(jj)
                        flag = true;
                        break;
                    end
                end
                if flag==false
                    return;
                end
            end
            if ~(obj.probe==obj2.probe)
                return;
            end
            if length(obj.aux)~=length(obj2.aux)
                return;
            end
            for ii=1:length(obj.aux)
                if ~(obj.aux(ii)==obj2.aux(ii))
                    return;
                end
            end
            if length(obj.metaDataTags)~=length(obj2.metaDataTags)
                return;
            end
            for ii=1:length(obj.metaDataTags)
                if ~(obj.metaDataTags(ii)==obj2.metaDataTags(ii))
                    return;
                end
            end
            B = true;
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Basic methods to Set/Get native variable
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        % ---------------------------------------------------------
        function val = GetFormatVersion(obj)
            val = obj.formatVersion;
        end
        
        % ---------------------------------------------------------
        function val = GetFormatVersionString(obj)
            val = sprintf('SNIRF v%s', obj.formatVersion);
        end
        
        % ---------------------------------------------------------
        function SetData(obj, val)
            obj.data = CopyHandles(val);
        end
        
        % ---------------------------------------------------------
        function val = GetData(obj)
            val = obj.data;
        end
        
        % ---------------------------------------------------------
        function SetStim(obj, val)
            obj.stim = CopyHandles(val);
        end
        
        % ---------------------------------------------------------
        function val = GetStim(obj)
            val = obj.stim;
        end
        
        % ---------------------------------------------------------
        function SetSd(obj, val)
            obj.probe = CopyHandles(val);
        end
        
        % ---------------------------------------------------------
        function val = GetSd(obj)
            val = obj.probe;
        end
        
        % ---------------------------------------------------------
        function SetAux(obj, val)
            obj.aux = CopyHandles(val);
        end
        
        % ---------------------------------------------------------
        function val = GetAux(obj)
            val = obj.aux;
        end
        
        % ---------------------------------------------------------
        function SetMetaDataTags(obj, val)
            obj.metaDataTags = val;
        end
        
        % ---------------------------------------------------------
        function val = GetMetaDataTags(obj)
            val = obj.metaDataTags;
            if isempty(obj)
                return;
            end
            if isempty(obj.metaDataTags)
                return;
            end
            val = obj.metaDataTags.Get();
        end
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods that must be implemented as a child class of AcqDataClass
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        % ---------------------------------------------------------
        function t = GetTime(obj, iBlk)
            t = [];
            if ~exist('iBlk','var') || isempty(iBlk)
                iBlk=1;
            end
            if iBlk>length(obj.data)
                return;
            end
            t = obj.data(iBlk).GetTime();
        end
        
        
        % ---------------------------------------------------------
        function datamat = GetDataTimeSeries(obj, options, iBlk)
            datamat = [];
            if ~exist('options','var')
                options = '';
            end
            if ~exist('iBlk','var') || isempty(iBlk)
                iBlk=1;
            end
            if iBlk>length(obj.data)
                return;
            end
            datamat = obj.data(iBlk).GetDataTimeSeries(options);
        end
        
        
        % ---------------------------------------------------------
        function datamat = GetAuxDataMatrix(obj)
            datamat = [];
            if isempty(obj.aux)
                return;
            end
            for ii=1:length(obj.aux)
                datamat(:,ii) = obj.aux(ii).GetDataTimeSeries();
            end
        end
        
        
        % ---------------------------------------------------------
        function names = GetAuxNames(obj)
            names = {};
            for ii=1:length(obj.aux)
                names{ii} = obj.aux(ii).GetName();
            end
        end
        
        
        % ---------------------------------------------------------
        function aux = GetAuxiliary(obj)
            aux = struct('names',{{}}, 'data', obj.aux);
            aux.names = obj.GetAuxNames();
            aux.data = obj.GetAuxDataMatrix();
        end
        
        
        % ---------------------------------------------------------
        function ml = GetMeasList(obj, iBlk)
            ml = [];
            if ~exist('iBlk','var') || isempty(iBlk)
                iBlk=1;
            end
            if iBlk>length(obj.data)
                return;
            end
            ml = obj.data(iBlk).GetMeasList();
        end
        
        
        % ---------------------------------------------------------
        function wls = GetWls(obj)
            wls = obj.probe.GetWls();
        end
        
        
        % ---------------------------------------------------------
        function SetStims_MatInput(obj, s, t)
            if nargin<2
                return
            end
            if isempty(t)
                return;
            end
            for ii=1:size(s,2)
                tidxs = find(s(:,ii)~=0);
                for jj=1:length(tidxs)
                    if ~obj.stim(ii).Exists(t(tidxs(jj)))
                        obj.stim(ii).AddStims(t(tidxs(jj)));
                    else
                        obj.stim(ii).EditValue(t(tidxs(jj)), s(tidxs(jj),ii));
                    end
                end
            end
        end
        
        
        % ---------------------------------------------------------
        function s = GetStims(obj, t)
            s = zeros(length(t), length(obj.stim));
            for ii=1:length(obj.stim)
                [ts, v] = obj.stim(ii).GetStim();
                [~, k] = nearest_point(t, ts);
                if isempty(k)
                    continue;
                end
                s(k,ii) = v;
            end
        end
        
        
        % ----------------------------------------------------------------------------------
        function SetConditions(obj, CondNames)
            if nargin==1
                return;
            end
            CondNamesLocal = unique({obj.stim.name});
            stimnew = StimClass().empty;
            for ii=1:length(CondNames)
                k = find(strcmp(CondNamesLocal, CondNames{ii}));
                if ~isempty(k)
                    stimnew(ii) = StimClass(obj.stim(k));
                else
                    stimnew(ii) = StimClass(CondNames{ii});
                end
            end
            obj.stim = stimnew;
        end
        
        
        % ---------------------------------------------------------
        function CondNames = GetConditions(obj)
            CondNames = cell(1,length(obj.stim));
            for ii=1:length(obj.stim)
                CondNames{ii} = obj.stim(ii).GetName();
            end
        end
        
        
        
        % ---------------------------------------------------------
        function SD = GetSDG(obj)
            SD = [];
            if isempty(obj)
                return;
            end
            if isempty(obj.probe)
                return;
            end
            SD.Lambda = obj.probe.GetWls();
            SD.SrcPos = obj.probe.GetSrcPos();
            SD.DetPos = obj.probe.GetDetPos();
        end
        
        
        % ---------------------------------------------------------
        function srcpos = GetSrcPos(obj)
            srcpos = obj.probe.GetSrcPos();
        end
        
        
        % ---------------------------------------------------------
        function detpos = GetDetPos(obj)
            detpos = obj.probe.GetDetPos();
        end
        
        
        % ----------------------------------------------------------------------------------
        function n = GetDataBlocksNum(obj)
            n = length(obj.data);
        end
        
        
        % ----------------------------------------------------------------------------------
        function [iDataBlks, ich] = GetDataBlocksIdxs(obj, ich0)
            iDataBlks=[];
            ich={};
            if nargin==1
                ich0=[];
            end
            if isempty(ich0)
                iDataBlks=1:length(obj.data);
                return;
            end
            
            % Get channel matrix for whole probe
            mlAll = [];
            nDataBlks = length(obj.data);
            for iBlk = 1:nDataBlks
                mlAll = [mlAll; obj.GetMeasList(iBlk)];
            end
            
            iSrc = mlAll(ich0,1);
            iDet = mlAll(ich0,2);
            
            % Now search block by block for the selecdted channels
            ich = cell(nDataBlks,1);
            for iBlk=1:nDataBlks
                ml = obj.GetMeasList(iBlk);
                for ii=1:length(ich0)
                    k = find(ml(:,1)==iSrc(ii) & ml(:,2)==iDet(ii));
                    if ~isempty(k)
                        iDataBlks = [iDataBlks; iBlk];
                        ich{iBlk} = [ich{iBlk}, k(1)];
                    end
                end
            end
            
            % Important: make sure iDataBlks is row vector (: + transpose does that) .
            % For some reason a for-loop traversing through empty column vector doesn't work properly
            iDataBlks = sort(unique(iDataBlks(:)'));
            
        end
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pubic interface for .nirs processing stream
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        % ----------------------------------------------------------------------------------
        function d = Get_d(obj, iBlk)
            d = [];
            if ~exist('iBlk','var') || isempty(iBlk)
                iBlk = 1;
            end
            if iBlk>length(obj.data)
                return;
            end
            d = obj.data(iBlk).GetDataTimeSeries();
        end
        
        
        % ----------------------------------------------------------------------------------
        function t = Get_t(obj, iBlk)
            t = [];
            if ~exist('iBlk','var') || isempty(iBlk)
                iBlk = 1;
            end
            if iBlk>length(obj.data)
                return;
            end
            t = obj.data(iBlk).GetTime();
        end
        
        
        % ----------------------------------------------------------------------------------
        function SD = Get_SD(obj, iBlk)
            SD = [];
            if isempty(obj.probe)
                return;
            end
            if ~exist('iBlk','var') || isempty(iBlk)
                iBlk = 1;
            end
            if iBlk>length(obj.data)
                return;
            end
            SD.Lambda   = obj.probe.GetWls();
            SD.SrcPos   = obj.probe.GetSrcPos();
            SD.DetPos   = obj.probe.GetDetPos();
            SD.MeasList = obj.data(iBlk).GetMeasList();
            SD.MeasListAct = ones(size(SD.MeasList,1),1);
        end
        
        
        
        % ----------------------------------------------------------------------------------
        function s = Get_s(obj, iBlk)
            s = [];
            if ~exist('iBlk','var') || isempty(iBlk)
                iBlk = 1;
            end
            if iBlk>length(obj.data)
                return;
            end
            t = obj.data(iBlk).GetTime();
            s = zeros(length(t), length(obj.stim));
            for ii=1:length(obj.stim)
                [ts, v] = obj.stim(ii).GetStim();
                [~, k] = nearest_point(t, ts);
                if isempty(k)
                    continue;
                end
                s(k,ii) = v;
            end
        end
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % All other public methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        % ----------------------------------------------------------------------------------
        function AddStims(obj, tPts, condition)
            % Try to find existing condition to which to add stims.
            for ii=1:length(obj.stim)
                if strcmp(condition, obj.stim(ii).GetName())
                    obj.stim(ii).AddStims(tPts);
                    return;
                end
            end
            
            % Otherwise we have a new condition to which to add the stims.
            obj.stim(end+1) = StimClass(tPts, condition);
            obj.SortStims();
        end
        
        
        % ----------------------------------------------------------------------------------
        function DeleteStims(obj, tPts, condition)
            % Find all stims for any conditions which match the time points.
            for ii=1:length(obj.stim)
                obj.stim(ii).DeleteStims(tPts);
            end
        end
        
        
        % ----------------------------------------------------------------------------------
        function ToggleStims(obj, tPts, condition)
            % Find all stims for any conditions which match the time points.
            for ii=1:length(obj.stim)
                obj.stim(ii).ToggleStims(tPts);
            end
        end
        
        
        % ----------------------------------------------------------------------------------
        function MoveStims(obj, tPts, condition)
            if ~exist('tPts','var') || isempty(tPts)
                return;
            end
            if ~exist('condition','var') || isempty(condition)
                return;
            end
            
            % Find the destination condition to move stims (among the time pts in tPts)
            % to
            j = [];
            for ii=1:length(obj.stim)
                if strcmp(condition, obj.stim(ii).GetName())
                    j=ii;
                    break;
                end
            end
            
            % If no destination condition found among existing conditions,
            % then create a new condition to move stims to
            if isempty(j)
                j = length(obj.stim)+1;
                
                % Otherwise we have a new condition to which to add the stims.
                obj.stim(j) = StimClass([], condition);
                obj.SortStims();
                
                % Recalculate j after sort
                for ii=1:length(obj.stim)
                    if strcmp(condition, obj.stim(ii).GetName())
                        j=ii;
                        break;
                    end
                end
            end
            
            % Find all stims for any conditions which match the time points.
            for ii=1:length(tPts)
                for kk=1:length(obj.stim)
                    d = obj.stim(kk).GetData();
                    if isempty(d)
                        continue;
                    end
                    k = find(d(:,1)==tPts(ii));
                    if ~isempty(k)
                        if kk==j
                            continue;
                        end
                        
                        % If stim at time point tPts(ii) exists in stim
                        % condition kk, then move stim from obj.stim(kk) to
                        % obj.stim(j)
                        obj.stim(j).AddStims(tPts(ii), d(k(1),2), d(k(1),3));
                        
                        % After moving stim from obj.stim(kk) to
                        % obj.stim(j), delete it from obj.stim(kk)
                        d(k(1),:)=[];
                        obj.stim(kk).SetData(d);
                        
                        % Move on to next time point
                        break;
                    end
                end
            end
        end
        
        
        % ----------------------------------------------------------------------------------
        function SetStimTpts(obj, icond, tpts)
            obj.stim(icond).SetTpts(tpts);
        end
        
        
        % ----------------------------------------------------------------------------------
        function tpts = GetStimTpts(obj, icond)
            if icond>length(obj.stim)
                tpts = [];
                return;
            end
            tpts = obj.stim(icond).GetTpts();
        end
        
        
        % ----------------------------------------------------------------------------------
        function SetStimDuration(obj, icond, duration)
            obj.stim(icond).SetDuration(duration);
        end
        
        
        % ----------------------------------------------------------------------------------
        function duration = GetStimDuration(obj, icond)
            if icond>length(obj.stim)
                duration = [];
                return;
            end
            duration = obj.stim(icond).GetDuration();
        end
        
        
        % ----------------------------------------------------------------------------------
        function SetStimValues(obj, icond, vals)
            obj.stim(icond).SetValues(vals);
        end
        
        
        
        % ----------------------------------------------------------------------------------
        function vals = GetStimValues(obj, icond)
            if icond>length(obj.stim)
                vals = [];
                return;
            end
            vals = obj.stim(icond).GetValues();
        end
        
        
        % ----------------------------------------------------------------------------------
        function RenameCondition(obj, oldname, newname)
            if ~exist('oldname','var') || ~ischar(oldname)
                return;
            end
            if ~exist('newname','var')  || ~ischar(newname)
                return;
            end
            k=[];
            for ii=1:length(obj.stim)
                if strcmp(obj.stim(ii).GetName(), oldname)
                    k = ii;
                    break;
                end
            end
            if isempty(k)
                return;
            end
            obj.stim(k).SetName(newname);
            obj.SortStims();
        end
        
        
        % ----------------------------------------------------------------------------------
        function b = IsEmpty(obj)
            b = true;
            if isempty(obj)
                return;
            end
            if isempty(obj.data) || (isa(obj.data, 'DataClass') && obj.data(1).IsEmpty())
                return;
            end
            if isempty(obj.data) || (isa(obj.probe, 'ProbeClass') && obj.probe.IsEmpty())
                return;
            end
            b = false;
        end
        
        
        % ----------------------------------------------------------------------------------
        function nbytes = MemoryRequired(obj)
            nbytes = 0;
            nbytes = nbytes + sizeof(obj.formatVersion);
            for ii=1:length(obj.metaDataTags)
                nbytes = nbytes + obj.metaDataTags(ii).MemoryRequired();
            end
            for ii=1:length(obj.data)
                nbytes = nbytes + obj.data(ii).MemoryRequired();
            end
            for ii=1:length(obj.stim)
                nbytes = nbytes + obj.stim(ii).MemoryRequired();
            end
            if ~isempty(obj.probe)
                nbytes = nbytes + obj.probe.MemoryRequired();
            end
            for ii=1:length(obj.aux)
                nbytes = nbytes + obj.aux(ii).MemoryRequired();
            end
        end
        
        
        
        % ----------------------------------------------------------------------------------
        function Info(obj)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Read formatVersion
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('\n');
            fprintf('    FormatVersion:\n');
            fv = obj.GetFormatVersion();
            fprintf('        Format version: %s\n', fv);
            fprintf('\n');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Load meta data tags from file and extract the tag names and values for display
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('    MetaDataTags:\n');
            tags = obj.GetMetaDataTags();
            for ii=1:length(tags)
                fprintf('        Tag #%d: {''%s'', ''%s''}\n', ii, tags(ii).key, tags(ii).value);
            end
            fprintf('\n');
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Load data from file and extract .nirs-style d and ml matrices for display
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('    Data (.nirs-style display):\n');
            for ii=1:length(obj.data)
                
                % Display data matrix dimensions and data type
                d = obj.data(ii).GetDataTimeSeries();
                pretty_print_struct(d, 8, 1);
                
                % Display meas list dimensions and data type
                ml = obj.data(ii).GetMeasList();
                pretty_print_struct(ml, 8, 1);
                
            end
            fprintf('\n');
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Load probe and extract .nirs-style SD structure
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('    Probe (.nirs-style display):\n');
            SD = obj.GetSDG();
            pretty_print_struct(SD, 8, 1);
            fprintf('\n');
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Load stim from file and extract it for display
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('    Stim (.snirf-style display):\n');
            for ii=1:length(obj.stim)
                fprintf('        stim(%d): {name = ''%s'', data = [', ii, obj.stim(ii).name);
                for jj=1:size(obj.stim(ii).data,1)
                    if jj==size(obj.stim(ii).data,1)
                        fprintf('%0.1f', obj.stim(ii).data(jj,1));
                    else
                        fprintf('%0.1f, ', obj.stim(ii).data(jj,1));
                    end
                end
                fprintf(']}\n');
            end
            fprintf('\n');
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Load aux from file and extract nirs-style data for display
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('    Aux (.nirs-style display):\n');
            auxl = obj.GetAuxiliary();
            pretty_print_struct(auxl, 8, 1);
            fprintf('\n');
            
        end
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Private methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = private)
        
        function GenSimulatedTimeBases(obj, nirs, tfactors)
            obj.nirs_tb = struct('SD', struct('Lambda',nirs.SD.Lambda, 'MeasList',[], 'SrcPOos',[], 'DetPOos',[]), ...
                't',[], ...
                'd',[] ...
                );
            
            if length(tfactors)==1 && tfactors==1
                obj.nirs_tb = nirs;
                return;
            end
            
            % a) Subdivide data and measurement list among time bases
            % b) Resample time and data
            % c) Throw out source and detectors that don't belong in each time base
            nCh = size(nirs.SD.MeasList,1)/length(nirs.SD.Lambda);
            nTimeBases = length(tfactors);
            baseSize = round(nCh/nTimeBases);
            irows = 1:baseSize:nCh;
            
            % Assign channels for time bases
            obj.nirs_tb = repmat(obj.nirs_tb, nTimeBases,1);
            for iWl=1:length(nirs.SD.Lambda)
                iBase = 1;
                for ii=irows+(iWl-1)*nCh
                    istart = ii;
                    if ii+baseSize-1 <= iWl*nCh
                        iend = ii+baseSize-1;
                    else
                        iend = iWl*nCh;
                    end
                    nChAdd = iend-istart+1;
                    obj.nirs_tb(iBase).d(:,end+1:end+nChAdd) = nirs.d(:,istart:iend);
                    obj.nirs_tb(iBase).SD.MeasList(end+1:end+nChAdd,:) = nirs.SD.MeasList(istart:iend,:);
                    if iBase<nTimeBases
                        iBase = iBase+1;
                    end
                end
            end
            
            % Resample data time and throw out optodes that don't belong in each time base
            for iBase=1:length(obj.nirs_tb)
                % Resample time
                [n,d] = rat(tfactors(iBase));
                obj.nirs_tb(iBase).t = resample(nirs.t, n, d);
                
                % Resample data
                obj.nirs_tb(iBase).d = resample(obj.nirs_tb(iBase).d, n, d);
                
                % Throw out source and detectors that don't belong in each time base
                iSrc = unique(obj.nirs_tb(iBase).SD.MeasList(:,1));
                obj.nirs_tb(iBase).SD.SrcPos = nirs.SD.SrcPos(iSrc);
                iDet = unique(obj.nirs_tb(iBase).SD.MeasList(:,2));
                obj.nirs_tb(iBase).SD.DetPos = nirs.SD.DetPos(iDet);
            end
            
        end
        
    end
    
end


