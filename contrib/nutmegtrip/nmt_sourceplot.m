function nmt_sourceplot(cfg,functional)

% NMT_SOURCEPLOT
% plots functional source reconstruction data on slices or on
% a surface, optionally as an overlay on anatomical MRI data, where
% statistical data can be used to determine the opacity of the mask. Input
% data comes from FT_SOURCEANALYSIS, FT_SOURCEGRANDAVERAGE or statistical
% values from FT_SOURCESTATISTICS.
%
% Use as
%   ft_sourceplot(cfg, data)
% where the input data can contain an anatomical MRI, functional source
% reconstruction results and/or statistical data. Interpolation is not
% necessary for this function; performance is best [no "gaps" or overlaps
% displayed activation map] if the functional data consists of a uniform
% grid in the chosen coordinate space. If the functional data is evenly
% spaced in MNI coordinates, for example, the data is best plotted on the
% MNI brain or the subject's MNI-warped MRI.
%
%
% The configuration should contain:
%   cfg.funparameter  = string, field in data with the functional parameter of interest (default = [])
%   cfg.mripath = string, location of Nifti-format MRI
%   cfg.maskparameter = string, field in the data to be used for opacity masking of fun data (default = [])
%                        If values are between 0 and 1, zero is fully transparant and one is fully opaque.
%                        If values in the field are not between 0 and 1 they will be scaled depending on the values
%                        of cfg.opacitymap and cfg.opacitylim (see below)
%                        You can use masking in several ways, f.i.
%                        - use outcome of statistics to show only the significant values and mask the insignificant
%                          NB see also cfg.opacitymap and cfg.opacitylim below
%                        - use the functional data itself as mask, the highest value (and/or lowest when negative)
%                          will be opaque and the value closest to zero transparent
%                        - Make your own field in the data with values between 0 and 1 to control opacity directly
%
% The following parameters can be used in all methods:
%   cfg.atlas         = string, filename of atlas to use (default = []) see FT_READ_ATLAS
%                        for ROI masking (see "masking" below) or in "ortho-plotting" mode (see "ortho-plotting" below)
%
% The following parameters can be used for the functional data:
%   **TODO** cfg.funcolormap   = colormap for functional data, see COLORMAP (default = 'auto')
%                       'auto', depends structure funparameter, or on funcolorlim
%                         - funparameter: only positive values, or funcolorlim:'zeromax' -> 'hot'
%                         - funparameter: only negative values, or funcolorlim:'minzero' -> 'cool'
%                         - funparameter: both pos and neg values, or funcolorlim:'maxabs' -> 'default'
%                         - funcolorlim: [min max] if min & max pos-> 'hot', neg-> 'cool', both-> 'default'
%   **TODO** cfg.funcolorlim   = color range of the functional data (default = 'auto')
%                        [min max]
%                        'maxabs', from -max(abs(funparameter)) to +max(abs(funparameter))
%                        'zeromax', from 0 to max(funparameter)
%                        'minzero', from min(funparameter) to 0
%                        'auto', if funparameter values are all positive: 'zeromax',
%                          all negative: 'minzero', both possitive and negative: 'maxabs'
%
% The following parameters can be used for the masking data:
%   **TODO** cfg.opacitymap    = opacitymap for mask data, see ALPHAMAP (default = 'auto')
%                       'auto', depends structure maskparameter, or on opacitylim
%                         - maskparameter: only positive values, or opacitylim:'zeromax' -> 'rampup'
%                         - maskparameter: only negative values, or opacitylim:'minzero' -> 'rampdown'
%                         - maskparameter: both pos and neg values, or opacitylim:'maxabs' -> 'vdown'
%                         - opacitylim: [min max] if min & max pos-> 'rampup', neg-> 'rampdown', both-> 'vdown'
%                         - NB. to use p-values use 'rampdown' to get lowest p-values opaque and highest transparent
%   **TODO** cfg.opacitylim    = range of mask values to which opacitymap is scaled (default = 'auto')
%                        [min max]
%                        'maxabs', from -max(abs(maskparameter)) to +max(abs(maskparameter))
%                        'zeromax', from 0 to max(abs(maskparameter))
%                        'minzero', from min(abs(maskparameter)) to 0
%                        'auto', if maskparameter values are all positive: 'zeromax',
%                          all negative: 'minzero', both possitive and negative: 'maxabs'
%   **TODO** cfg.roi           = string or cell of strings, region(s) of interest from anatomical atlas (see cfg.atlas above)
%                        everything is masked except for ROI
%
% The following parameters apply for ortho-plotting
%   **TODO** cfg.location      = location of cut, (default = 'auto')
%                        'auto', 'center' if only anatomy, 'max' if functional data
%                        'min' and 'max' position of min/max funparameter
%                        'center' of the brain
%                        [x y z], coordinates in voxels or head, see cfg.locationcoordinates
%   **TODO** cfg.locationcoordinates = coordinate system used in cfg.location, 'head' or 'voxel' (default = 'head')
%                              'head', headcoordinates as mm or cm
%                              'voxel', voxelcoordinates as indices
%   **TODO** cfg.crosshair     = 'yes' or 'no' (default = 'yes')
%   **TODO** cfg.axis          = 'on' or 'off' (default = 'on')
%   **TODO** cfg.queryrange    = number, in atlas voxels (default 3)
%
%
% The following parameters apply for slice-plotting
%   **TODO** cfg.title         = string, title of the figure window
%
% **TODO**
% When cfg.method = 'surface', the functional data will be rendered onto a
% cortical mesh (can be an inflated mesh). If the input source data
% contains a tri-field, no interpolation is needed. If the input source
% data does not contain a tri-field (i.e. a description of a mesh), an
% interpolation is performed onto a specified surface. Note that the
% coordinate system in which the surface is defined should be the same as
% the coordinate system that is represented in source.pos.
%
% The following parameters apply to surface-plotting when an interpolation
% is required
%   cfg.surffile       = string, file that contains the surface (default = 'surface_white_both.mat')
%                        'surface_white_both.mat' contains a triangulation that corresponds with the
%                         SPM anatomical template in MNI coordinates
%   cfg.surfinflated   = string, file that contains the inflated surface (default = [])
%                        may require specifying a point-matching (uninflated) surffile
%   cfg.surfdownsample = number (default = 1, i.e. no downsampling)
%   cfg.projmethod     = projection method, how functional volume data is projected onto surface
%                        'nearest', 'project', 'sphere_avg', 'sphere_weighteddistance'
%   cfg.projvec        = vector (in mm) to allow different projections that
%                        are combined with the method specified in cfg.projcomb
%   cfg.projcomb       = 'mean', 'max', method to combine the different projections
%   cfg.projweight     = vector of weights for the different projections (default = 1)
%   cfg.projthresh     = implements thresholding on the surface level
%                        for example, 0.7 means 70% of maximum
%   cfg.sphereradius   = maximum distance from each voxel to the surface to be
%                        included in the sphere projection methods, expressed in mm
%   cfg.distmat        = precomputed distance matrix (default = [])
%
% The following parameters apply to surface-plotting independent of whether
% an interpolation is required
%   cfg.camlight       = 'yes' or 'no' (default = 'yes')
%   cfg.renderer       = 'painters', 'zbuffer',' opengl' or 'none' (default = 'opengl')
%                        note that when using opacity the OpenGL renderer is required.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure.
%
% See also FT_SOURCEANALYSIS, FT_SOURCEGRANDAVERAGE, FT_SOURCESTATISTICS,
% FT_VOLUMELOOKUP, FT_READ_ATLAS, FT_READ_MRI

% TODO have to be built in:
%   cfg.marker        = [Nx3] array defining N marker positions to display (orig: from sliceinterp)
%   cfg.markersize    = radius of markers (default = 5)
%   cfg.markercolor   = [1x3] marker color in RGB (default = [1 1 1], i.e. white) (orig: from sliceinterp)
%   white background option

%
% Author: Sarang S. Dalal
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

%% SPM check
if(~exist('spm_image','file'))
    error('Nutmegtrip requires a full version of SPM12 or SPM8 in your MATLAB path (fieldtrip/external/spm8 does not suffice). It may be downloaded from http://www.fil.ion.ucl.ac.uk/spm')
end

%%
% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble provenance
ft_preamble trackconfig
ft_preamble debug
ft_preamble loadvar functional

% the abort variable is set to true or false in ft_preamble_init
if ft_abort
    return
end

% this is not supported any more as of 26/10/2011
if ischar(functional)
    error('please use cfg.inputfile instead of specifying the input variable as a string');
end

    % set the defaults for all methods
    cfg.funparameter  = ft_getopt(cfg, 'funparameter',  []);
    cfg.vecparameter  = ft_getopt(cfg, 'vecparameter',  []);
    cfg.maskparameter = ft_getopt(cfg, 'maskparameter', []);
    
        
    if isfield(cfg, 'atlas') && ~isempty(cfg.atlas)
        % the atlas lookup requires the specification of the coordsys
        functional     = ft_checkdata(functional, 'datatype', {'volume', 'source'}, 'feedback', 'yes', 'hasunit', 'yes', 'hascoordsys', 'yes');
    else
        % check if the input functional is valid for this function, a coordsys is not directly needed
        functional     = ft_checkdata(functional, 'datatype', {'volume', 'source'}, 'feedback', 'yes', 'hasunit', 'yes');
    end
    
    if(isfield(functional,'dim'))
        functional = rmfield(functional,'dim');
    end
    
    % determine the type of functional
    issource = ft_datatype(functional, 'source');
    isvolume = ft_datatype(functional, 'volume');
    
    % set the defaults for all methods
    cfg.method        = ft_getopt(cfg, 'method',        'ortho');
    cfg.funparameter  = ft_getopt(cfg, 'funparameter',  []);
    cfg.maskparameter = ft_getopt(cfg, 'maskparameter', []);
    cfg.downsample    = ft_getopt(cfg, 'downsample',    1);
    cfg.title         = ft_getopt(cfg, 'title',         '');
    cfg.atlas         = ft_getopt(cfg, 'atlas',         []);
    cfg.marker        = ft_getopt(cfg, 'marker',        []);
    cfg.markersize    = ft_getopt(cfg, 'markersize',    5);
    cfg.markercolor   = ft_getopt(cfg, 'markercolor',   [1 1 1]);
    
    if ~isfield(cfg, 'anaparameter')
        if isfield(functional, 'anatomy')
            cfg.anaparameter = 'anatomy';
        else
            cfg.anaparameter = [];
        end
    end
    
    % set the common defaults for the functional data
    cfg.funcolormap   = ft_getopt(cfg, 'funcolormap',   'auto');
    cfg.funcolorlim   = ft_getopt(cfg, 'funcolorlim',   'auto');
    
    % set the common defaults for the statistical data
    cfg.opacitymap    = ft_getopt(cfg, 'opacitymap',    'auto');
    cfg.opacitylim    = ft_getopt(cfg, 'opacitylim',    'auto');
    cfg.roi           = ft_getopt(cfg, 'roi',           []);
    
    % set the defaults per method
    
    % ortho
    cfg.location            = ft_getopt(cfg, 'location',            'auto');
    cfg.locationcoordinates = ft_getopt(cfg, 'locationcoordinates', 'head');
    cfg.crosshair           = ft_getopt(cfg, 'crosshair',           'yes');
    cfg.colorbar            = ft_getopt(cfg, 'colorbar',            'yes');
    cfg.axis                = ft_getopt(cfg, 'axis',                'on');
    cfg.queryrange          = ft_getopt(cfg, 'queryrange',          3);
    
    if isfield(cfg, 'TTlookup'),
        error('TTlookup is old; now specify cfg.atlas, see help!');
    end
    
    % slice
    cfg.nslices    = ft_getopt(cfg, 'nslices',    20);
    cfg.slicedim   = ft_getopt(cfg, 'slicedim',   3);
    cfg.slicerange = ft_getopt(cfg, 'slicerange', 'auto');
    
    % surface
    cfg.downsample     = ft_getopt(cfg, 'downsample',     1);
    cfg.surfdownsample = ft_getopt(cfg, 'surfdownsample', 1);
    cfg.surffile       = ft_getopt(cfg, 'surffile', 'surface_white_both.mat'); % use a triangulation that corresponds with the collin27 anatomical template in MNI coordinates
    cfg.surfinflated   = ft_getopt(cfg, 'surfinflated',  []);
    cfg.sphereradius   = ft_getopt(cfg, 'sphereradius',  []);
    cfg.projvec        = ft_getopt(cfg, 'projvec',       1);
    cfg.projweight     = ft_getopt(cfg, 'projweight',    ones(size(cfg.projvec)));
    cfg.projcomb       = ft_getopt(cfg, 'projcomb',      'mean'); % or max
    cfg.projthresh     = ft_getopt(cfg, 'projthresh',    []);
    cfg.projmethod     = ft_getopt(cfg, 'projmethod',    'nearest');
    cfg.distmat        = ft_getopt(cfg, 'distmat',       []);
    cfg.camlight       = ft_getopt(cfg, 'camlight',      'yes');
    cfg.renderer       = ft_getopt(cfg, 'renderer',      'opengl');
    % if isequal(cfg.method,'surface')
    % if ~isfield(cfg, 'projmethod'),
    % error('specify cfg.projmethod');
    % end
    % end
    
    % for backward compatibility
    if strcmp(cfg.location, 'interactive')
        cfg.location = 'auto';
    end
    
    % ensure that old and unsupported options are not being relied on by the end-user's script
    % instead of specifying cfg.coordsys, the user should specify the coordsys in the functional data
    cfg = ft_checkconfig(cfg, 'forbidden', {'units', 'inputcoordsys', 'coordinates'});
    cfg = ft_checkconfig(cfg, 'deprecated', 'coordsys');

    cfg.title         = ft_getopt(cfg, 'title',         '');
    cfg.atlas         = ft_getopt(cfg, 'atlas',         []);
    cfg.topoplot      = ft_getopt(cfg, 'topoplot',  '');


    %% start building the figure
    global st
    
    % things get messy if there's an SPM window already open
    spmfigh = spm_figure('FindWin');
    if(~isempty(spmfigh))
        close(spmfigh)
    end
    
    spm_image('init',cfg.mripath); % load/reload structural MRI
    nmt_spmfig_setup(cfg);
    

% if funparameter is a simple string, convert to 1-element cell
if iscell(cfg.funparameter)
    funparameters = cfg.funparameter;
else
    funparameters{1} = cfg.funparameter;
end

for funidx = 1:length(funparameters)
    cfg.funparameter = funparameters{funidx}
    % ensure that old and unsupported options are not being relied on by the end-user's script
    % instead of specifying cfg.coordsys, the user should specify the coordsys in the functional data
    cfg = ft_checkconfig(cfg, 'renamedval', {'funparameter', 'avg.pow', 'pow'});
    cfg = ft_checkconfig(cfg, 'renamedval', {'funparameter', 'avg.coh', 'coh'});
    cfg = ft_checkconfig(cfg, 'renamedval', {'funparameter', 'avg.mom', 'mom'});
    cfg = ft_checkconfig(cfg, 'renamedval', {'funparameter', 'avg.aa', 'aa'});
    cfg = ft_checkconfig(cfg, 'renamedval', {'funparameter', 'avg.itc', 'itc'});
    cfg = ft_checkconfig(cfg, 'renamedval', {'funparameter', 'avg.tf', 'tf'});
    cfg = ft_checkconfig(cfg, 'renamedval', {'maskparameter', 'avg.pow', 'pow'});
    cfg = ft_checkconfig(cfg, 'renamedval', {'maskparameter', 'avg.coh', 'coh'});
    cfg = ft_checkconfig(cfg, 'renamedval', {'maskparameter', 'avg.mom', 'mom'});
    cfg = ft_checkconfig(cfg, 'renamedval', {'maskparameter', 'avg.aa', 'aa'});
    cfg = ft_checkconfig(cfg, 'renamedval', {'maskparameter', 'avg.itc', 'itc'});
    cfg = ft_checkconfig(cfg, 'renamedval', {'maskparameter', 'avg.tf', 'tf'});
    
    
    
    % select the functional and the mask parameter
    cfg.funparameter  = parameterselection(cfg.funparameter, functional);
    cfg.maskparameter = parameterselection(cfg.maskparameter, functional);
    % only a single parameter should be selected
    try, cfg.funparameter  = cfg.funparameter{1};  end
    try, cfg.maskparameter = cfg.maskparameter{1}; end
    
    switch(getdimord(functional,cfg.funparameter))
        case '{pos}_freq_time'
            cfg.plottype      = ft_getopt(cfg, 'plottype', 'tf'); % default plot type is "time-freq"
        otherwise
            cfg.plottype      = ft_getopt(cfg, 'plottype', 'ts'); % default plot type is "time series"
    end
    
    if isvolume && cfg.downsample~=1
        % optionally downsample the anatomical and/or functional volumes
        tmpcfg = keepfields(cfg, {'downsample'});
        tmpcfg.parameter = {cfg.funparameter, cfg.maskparameter, cfg.anaparameter};
        functional = ft_volumedownsample(tmpcfg, functional);
        [cfg, functional] = rollback_provenance(cfg, functional);
    end
    
    %%% make the local variables:
    if isfield(functional, 'dim')
        dim = functional.dim;
    else
        dim = [size(functional.pos,1) 1];
    end
    
    hasatlas = ~isempty(cfg.atlas);
    if hasatlas
        if ischar(cfg.atlas)
            % initialize the atlas
            [p, f, x] = fileparts(cfg.atlas);
            fprintf(['reading ', f,' atlas coordinates and labels\n']);
            atlas = ft_read_atlas(cfg.atlas);
        else
            atlas = cfg.atlas;
        end
    end
    
    hasroi = ~isempty(cfg.roi);
    if hasroi
        if ~hasatlas
            error('specify cfg.atlas which belongs to cfg.roi')
        else
            % get the mask
            tmpcfg          = [];
            tmpcfg.roi      = cfg.roi;
            tmpcfg.atlas    = cfg.atlas;
            tmpcfg.inputcoord = functional.coordsys;
            roi = ft_volumelookup(tmpcfg,functional);
        end
    end
    
    % %%% anaparameter
    hasana = 1; % by definition, you's got ana if you're using this function :-)
    % if isempty(cfg.anaparameter);
    %   hasana = 0;
    %   fprintf('not plotting anatomy\n');
    % elseif isfield(functional, cfg.anaparameter)
    %   hasana = 1;
    %   ana = getsubfield(functional, cfg.anaparameter);
    %   % convert integers to single precision float if neccessary
    %   if isa(ana, 'uint8') || isa(ana, 'uint16') || isa(ana, 'int8') || isa(ana, 'int16')
    %     fprintf('converting anatomy to double\n');
    %     ana = double(ana);
    %   end
    %   fprintf('scaling anatomy to [0 1]\n');
    %   dmin = min(ana(:));
    %   dmax = max(ana(:));
    %   ana  = (ana-dmin)./(dmax-dmin);
    % else
    %   warning('do not understand cfg.anaparameter, not plotting anatomy\n')
    %   hasana = 0;
    % end
    
    %%% funparameter
    % has fun?
    if ~isempty(cfg.funparameter)
        if issubfield(functional, cfg.funparameter)
            hasfun = 1;
            tmpfun = getsubfield(functional, cfg.funparameter);
        else
            error('cfg.funparameter not found in functional');
        end
    else
        hasfun = 0;
        fprintf('no functional parameter\n');
    end
    
    % handle the dimensions of functional data
    
    if hasfun
        dimord = getdimord(functional, cfg.funparameter);
        dimtok = tokenize(dimord, '_');
        
        if strcmp(dimtok{1}, '{pos}')
            tmpdim = getdimsiz(functional, cfg.funparameter);
            fun = nan(tmpdim);
            insideindx = find(functional.inside);
            
            tmpfun = cell2mat(tmpfun);
            switch(dimord)
                case {'{pos}_freq_time','{pos}_ori_time','{pos}_unknown_time'} % ori is hack... FT thinks that 3 freq bands is an ori!
                    tmpfun = reshape(tmpfun,[size(tmpfun,1)/length(insideindx) length(insideindx) size(tmpfun,2)]);
                    tmpfun = permute(tmpfun, [2 1 3]);
            end
            
            %    tmpfun = real(tmpfun);
            %    tmpfun = (tmpfun ./ repmat(mean(abs(tmpfun(:,:,1:200)),3),1,1,476));
            
            %    tmpfun = sum(tmpfun,2);
            %    fun = fun(insideindx,1,:);
            
            fun(insideindx,:,:) = tmpfun; % replace the cell-array functional with a normal array
            clear tmpfun;
            dimtok{1} = 'pos';  % update the description of the dimensions
            dimord([1 5]) = []; % remove the { and }
        else
            fun = tmpfun;
            clear tmpfun;
        end
        
        
        if strcmp(dimord, 'pos_rgb')
            % treat functional data as rgb values
            if any(fun(:)>1 | fun(:)<0)
                % scale
                tmpdim = size(fun);
                nvox   = prod(tmpdim(1:end-1));
                tmpfun = reshape(fun,[nvox tmpdim(end)]);
                m1     = max(tmpfun,[],1);
                m2     = min(tmpfun,[],1);
                tmpfun = (tmpfun-m2(ones(nvox,1),:))./(m1(ones(nvox,1),:)-m2(ones(nvox,1),:));
                fun    = reshape(tmpfun, tmpdim);
                clear tmpfun
            end
            qi      = 1;
            hasfreq = 0;
            hastime = 0;
            
            doimage = 1;
            fcolmin = 0;
            fcolmax = 1;
            
        else
            % determine scaling min and max (fcolmin fcolmax) and funcolormap
            if ~isa(fun, 'logical')
                funmin = min(fun(:));
                funmax = max(fun(:));
            else
                funmin = 0;
                funmax = 1;
            end
            % smart automatic limits
            if isequal(cfg.funcolorlim,'auto')
                if sign(funmin)>-1 && sign(funmax)>-1
                    cfg.funcolorlim = 'zeromax';
                elseif sign(funmin)<1 && sign(funmax)<1
                    cfg.funcolorlim = 'minzero';
                else
                    cfg.funcolorlim = 'maxabs';
                end
            end
            if ischar(cfg.funcolorlim)
                % limits are given as string
                if isequal(cfg.funcolorlim,'maxabs')
                    fcolmin = -max(abs([funmin,funmax]));
                    fcolmax =  max(abs([funmin,funmax]));
                    if isequal(cfg.funcolormap,'auto'); cfg.funcolormap = 'default'; end;
                elseif isequal(cfg.funcolorlim,'zeromax')
                    fcolmin = 0;
                    fcolmax = funmax;
                    if isequal(cfg.funcolormap,'auto'); cfg.funcolormap = 'hot'; end;
                elseif isequal(cfg.funcolorlim,'minzero')
                    fcolmin = funmin;
                    fcolmax = 0;
                    if isequal(cfg.funcolormap,'auto'); cfg.funcolormap = 'cool'; end;
                else
                    error('do not understand cfg.funcolorlim');
                end
            else
                % limits are numeric
                fcolmin = cfg.funcolorlim(1);
                fcolmax = cfg.funcolorlim(2);
                % smart colormap
                if isequal(cfg.funcolormap,'auto')
                    if sign(fcolmin) == -1 && sign(fcolmax) == 1
                        cfg.funcolormap = 'default';
                    else
                        if fcolmin < 0
                            cfg.funcolormap = 'cool';
                        else
                            cfg.funcolormap = 'hot';
                        end
                    end
                end
            end % if ischar
            clear funmin funmax
            
            % FIXME should this not be done earlier in the code?
            % ensure that the functional data is real
            if ~isreal(fun)
                warning('functional data is complex, taking absolute value');
                fun = abs(fun);
            end
            
            if ndims(fun)>3 || prod(dim)==size(fun,1)
                if strcmp(dimord, 'pos_freq_time') || strcmp(dimord, 'pos_ori_time')
                    % functional contains time-frequency representation
                    qi      = [1 1];
                    hasfreq = numel(functional.freq)>1;
                    hastime = numel(functional.time)>1;
                    %fun     = reshape(fun, [dim numel(functional.freq) numel(functional.time)]);
                    
                    fun = permute(fun,[1 3 2]); % reorder to pos_time_freq
                    dimord = 'pos_time_freq';
                elseif strcmp(dimord, 'pos_time')
                    % functional contains evoked field
                    qi      = 1;
                    hasfreq = 0;
                    hastime = numel(functional.time)>1;
                    fun     = squeeze(fun);
                    %        fun     = reshape(fun, [dim numel(functional.time)]);
                elseif strcmp(dimord, 'pos_freq')
                    % functional contains frequency spectra
                    qi      = 1;
                    hasfreq = numel(functional.freq)>1;
                    hastime = 0;
                    fun     = reshape(fun, [dim numel(functional.freq)]);
                else
                    qi      = 1;
                    hasfreq = 0;
                    hastime = 0;
                    fun     = reshape(fun, dim);
                end
            else
                % do nothing
                qi      = 1;
                hasfreq = 0;
                hastime = 0;
            end
            
            doimage = 0;
        end % if dimord has rgb or something else
    else
        % there is no functional data
        qi      = 1;
        hasfreq = 0;
        hastime = 0;
        
        doimage = 0;
        fcolmin = 0; % needs to be defined for callback
        fcolmax = 1;
    end % handle fun
    
    %%% maskparameter
    % has mask?
    if ~isempty(cfg.maskparameter)
        if issubfield(functional, cfg.maskparameter)
            if ~hasfun
                error('you can not have a mask without functional data')
            else
                hasmsk = 1;
                msk = getsubfield(functional, cfg.maskparameter);
                if islogical(msk) % otherwise sign() not posible
                    msk = double(msk);
                end
            end
        else
            error('cfg.maskparameter not found in functional');
        end
    else
        hasmsk = 0;
        fprintf('no masking parameter\n');
    end
    
    % handle mask
    if hasmsk
        % reshape to match fun
        if strcmp(dimord, 'pos_time_freq')
            % functional contains timefrequency representation
            msk = permute(msk,[1 3 2]); % reorder to pos_time_freq
            msk     = reshape(msk, [dim numel(functional.time) numel(functional.freq)]);
        elseif strcmp(dimord, 'pos_time')
            % functional contains evoked field
            msk     = reshape(msk, [dim numel(functional.time)]);
        elseif strcmp(dimord, 'pos_freq')
            % functional contains frequency spectra
            msk     = reshape(msk, [dim numel(functional.freq)]);
        else
            msk     = reshape(msk, dim);
        end
        
        % determine scaling and opacitymap
        mskmin = min(msk(:));
        mskmax = max(msk(:));
        % determine the opacity limits and the opacity map
        % smart limits: make from auto other string, or equal to funcolorlim if funparameter == maskparameter
        if isequal(cfg.opacitylim,'auto')
            if isequal(cfg.funparameter,cfg.maskparameter)
                cfg.opacitylim = cfg.funcolorlim;
            else
                if sign(mskmin)>-1 && sign(mskmax)>-1
                    cfg.opacitylim = 'zeromax';
                elseif sign(mskmin)<1 && sign(mskmax)<1
                    cfg.opacitylim = 'minzero';
                else
                    cfg.opacitylim = 'maxabs';
                end
            end
        end
        if ischar(cfg.opacitylim)
            % limits are given as string
            switch cfg.opacitylim
                case 'zeromax'
                    opacmin = 0;
                    opacmax = mskmax;
                    if isequal(cfg.opacitymap,'auto'), cfg.opacitymap = 'rampup'; end;
                case 'minzero'
                    opacmin = mskmin;
                    opacmax = 0;
                    if isequal(cfg.opacitymap,'auto'), cfg.opacitymap = 'rampdown'; end;
                case 'maxabs'
                    opacmin = -max(abs([mskmin, mskmax]));
                    opacmax =  max(abs([mskmin, mskmax]));
                    if isequal(cfg.opacitymap,'auto'), cfg.opacitymap = 'vdown'; end;
                otherwise
                    error('incorrect specification of cfg.opacitylim');
            end % switch opacitylim
        else
            % limits are numeric
            opacmin = cfg.opacitylim(1);
            opacmax = cfg.opacitylim(2);
            if isequal(cfg.opacitymap,'auto')
                if sign(opacmin)>-1 && sign(opacmax)>-1
                    cfg.opacitymap = 'rampup';
                elseif sign(opacmin)<1 && sign(opacmax)<1
                    cfg.opacitymap = 'rampdown';
                else
                    cfg.opacitymap = 'vdown';
                end
            end
        end % handling opacitylim and opacitymap
        clear mskmin mskmax
    else
        opacmin = [];
        opacmax = [];
    end
    
    % prevent outside fun from being plotted
    if hasfun && isfield(functional,'inside') && ~hasmsk
        hasmsk = 1;
        msk = ones(size(fun));
        cfg.opacitymap = 'rampup';
        opacmin = 0;
        opacmax = 1;
        % make intelligent mask
        if isequal(cfg.method,'surface')
            msk(functional.inside) = 1;
        else
            if hasana
%                msk(functional.inside) = 0.5; % so anatomy is visible
% FIXME: this lets nutmegtrip display at proper colorscale, but is the 0.5 functionality desired?
                msk(functional.inside) = 1;
            else
                msk(functional.inside) = 1;
            end
        end
    end
    
    % if region of interest is specified, mask everything besides roi
    if hasfun && hasroi && ~hasmsk
        hasmsk = 1;
        msk = roi;
        cfg.opacitymap = 'rampup';
        opacmin = 0;
        opacmax = 1;
    elseif hasfun && hasroi && hasmsk
        msk = roi .* msk;
        opacmin = [];
        opacmax = []; % has to be defined
    elseif hasroi
        error('you can not have a roi without functional data')
    end
    
    
    
    %%% set color and opacity mapping for this figure
    if hasfun
        colormap(cfg.funcolormap);
        cfg.funcolormap = colormap;
    end
    if hasmsk
        cfg.opacitymap = alphamap(cfg.opacitymap);
        alphamap(cfg.opacitymap);
        if ndims(fun)>3 && ndims(msk)==3
            siz = size(fun);
            msk = repmat(msk, [1 1 1 siz(4:end)]);
        end
    end
    
    if(~isfield(cfg,'colormap'))
        cfg.colormap = jet(64);
    elseif(ischar(cfg.colormap))
        cfg.colormap = feval(cfg.colormap,64);
    end

    %% *********************************************************************************
    % SPM8 expects everything in mm
    functional = ft_convert_units(functional,'mm');
    
    
    if ~isempty(cfg.funparameter)
        if issubfield(functional, cfg.funparameter)
            hasfun = 1;
 
            cfg.inside_idx = find(functional.inside);
            
            st.nmt.pos = functional.pos;
            
            st.nmt.fun{funidx} = fun;
            clear fun;
            
            if(hastime && ~hasfreq)
                % voxels x time
                st.nmt.time = functional.time;
                if(isfield(functional,'freqbands'))
                    st.nmt.freq = functional.freqbands;
                else
                    st.nmt.freq = [0 inf];
                end
                
                if(~isfield(cfg,'time') && ~isfield(cfg,'vox'))
                    [~,peakind] = max(abs(st.nmt.fun{funidx}(:)));
                    [peakvox_idx,peaktime_idx] = ind2sub(size(st.nmt.fun{funidx}),peakind);
                    cfg.time_idx(1) = peaktime_idx;
                    cfg.vox_idx = peakvox_idx;
                end
                
                if(~isfield(cfg,'time') && isfield(cfg,'vox'))
                    [~,peaktime_idx] = max(abs(st.nmt.fun{funidx}(cfg.vox_idx,:)));
                    cfg.time_idx(1) = peaktime_idx;
                end
                
                if(isfield(cfg,'time') && ~isfield(cfg,'vox'))
                    [~,peakvox_idx] = max(abs(st.nmt.fun{funidx}(cfg.time_idx,:)));
                    cfg.vox_idx = peakvox_idx;
                end
                
                % move MRI crosshairs to desired/peak voxel
                spm_orthviews('Reposition',st.nmt.pos(cfg.vox_idx,:))
                
                if(length(cfg.time_idx(1)) == 1)
                    cfg.time_idx(2) = cfg.time_idx(1);
                end
                
                cfg.freq_idx = [1 1];
            elseif(hastime && hasfreq)
                % voxels x frequency x time
                st.nmt.time = functional.time;
                st.nmt.freq = functional.freqbands;
                
                if(~isfield(cfg,'time') && ~isfield(cfg,'vox'))
                    [~,peakind] = max(abs(st.nmt.fun{funidx}(:)));
                    [peakvox_idx,peaktime_idx,peakfreq_idx] = ind2sub(size(st.nmt.fun{funidx}),peakind);
                    cfg.time_idx(1) = peaktime_idx;
                    cfg.freq_idx(1) = peakfreq_idx;
                    cfg.vox_idx = peakvox_idx;
                end
                
                if(~isfield(cfg,'time') && isfield(cfg,'vox'))
                    [~,peaktime_idx] = max(abs(st.nmt.fun{funidx}(cfg.vox_idx,:)));
                    cfg.time_idx(1) = peaktime_idx;
                end
                
                if(isfield(cfg,'time') && ~isfield(cfg,'vox'))
                    [~,peakvox_idx] = max(abs(st.nmt.fun{funidx}(cfg.time_idx,:)));
                    cfg.vox_idx = peakvox_idx;
                end
                
                % move MRI crosshairs to desired/peak voxel
                spm_orthviews('Reposition',st.nmt.pos(cfg.vox_idx,:))
                
                if(length(cfg.time_idx(1)) == 1)
                    cfg.time_idx(2) = cfg.time_idx(1);
                end
                if(length(cfg.freq_idx(1)) == 1)
                    cfg.freq_idx(2) = cfg.freq_idx(1);
                end
            else
                cfg.time_idx = [1 1]; % no time dimension in this case, e.g., 'pow'
                cfg.freq_idx = [1 1]; % frequency dimension is singleton in this case
                
                st.nmt.freq = [0 inf]; % dummy frequencies to make later functions happy
            end
            
            set(st.nmt.gui.f1,'String',num2str(st.nmt.freq(:,1)));
            set(st.nmt.gui.f2,'String',num2str(st.nmt.freq(:,2)));
            
            st.nmt.cfg = cfg;
            st.nmt.msk = reshape(msk,size(st.nmt.fun{funidx}));
            
        else
            error('cfg.funparameter not found in functional');
        end
    else
        hasfun = 0;
        fprintf('no functional parameter\n');
    end
    
    switch(cfg.topoplot)
        case 'timelock'
            st.nmt.timelock = functional.timelock;
        case {'spatialfilter'}
            st.nmt.spatialfilter = functional.filter;
        case {'leadfield','leadfieldX','leadfieldY','leadfieldZ','leadfieldori'}
            st.nmt.grid = functional.grid;
    end
    
    nmt_spm_plot(cfg);
    nmt_update_panel(funidx);
    nmt_image;
    
    if ~isempty(cfg.vecparameter)
        if issubfield(functional, cfg.vecparameter)
            oritmp = getsubfield(functional, cfg.vecparameter);
            Nvoxels = length(insideindx)
            Nsamples = size(oritmp{insideindx(1)},2)
            ori = zeros(Nvoxels,3,Nsamples);
            for ii=1:Nvoxels
                ori(ii,:,:) = oritmp{insideindx(ii)};
            end
            
            st.nmt.ori = ori;
            nmt_sourceoriplot;
        else
            error('cfg.vecparameter not found in functional');
        end
    end
end