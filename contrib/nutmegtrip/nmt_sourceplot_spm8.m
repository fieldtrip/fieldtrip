function nmt_sourceplot_spm8(cfg,functional)
% FT_SOURCEPLOT_SPM8
% plots functional source reconstruction data on slices or on
% a surface, optionally as an overlay on anatomical MRI data, where
% statistical data can be used to determine the opacity of the mask. Input
% data comes from FT_SOURCEANALYSIS, FT_SOURCEGRANDAVERAGE or statistical
% values from FT_SOURCESTATISTICS.
%
% Use as
%   ft_sourceplot_spm8(cfg, data)
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
%   **TODO** cfg.maskparameter = string, field in the data to be used for opacity masking of fun data (default = [])
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
%   **TODO** cfg.atlas         = string, filename of atlas to use (default = []) see FT_READ_ATLAS
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
%   **TODO** cfg.colorbar      = 'yes' or 'no' (default = 'yes')
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble provenance
ft_preamble trackconfig
ft_preamble debug
ft_preamble loadvar functional

% the abort variable is set to true or false in ft_preamble_init
if abort
    return
end

% ensure that old and unsupported options are not being relied on by the end-user's script
% instead of specifying cfg.coordsys, the user should specify the coordsys in the functional data
cfg = ft_checkconfig(cfg, 'forbidden', {'units', 'inputcoordsys', 'coordinates'});
cfg = ft_checkconfig(cfg, 'deprecated', 'coordsys');
cfg = ft_checkconfig(cfg, 'renamedval', {'funparameter', 'avg.pow', 'pow'});
cfg = ft_checkconfig(cfg, 'renamedval', {'funparameter', 'avg.coh', 'coh'});
cfg = ft_checkconfig(cfg, 'renamedval', {'funparameter', 'avg.mom', 'mom'});
cfg = ft_checkconfig(cfg, 'renamedval', {'maskparameter', 'avg.pow', 'pow'});
cfg = ft_checkconfig(cfg, 'renamedval', {'maskparameter', 'avg.coh', 'coh'});
cfg = ft_checkconfig(cfg, 'renamedval', {'maskparameter', 'avg.mom', 'mom'});


functional = ft_checkdata(functional, 'datatype', 'source', 'hasunit', 'yes');

%% SPM8 expects everything in mm
functional = ft_convert_units(functional,'mm');

%%
global st

%if(isempty(spm_figure('FindWin')))
spm_image('init',cfg.mripath); % load/reload structural MRI
nmt_spmfig_setup(cfg)
%end

%cfg.maskparameter = 'mask';
%cfg.funparameter = 'pow'

if ~isempty(cfg.funparameter)
    if issubfield(functional, cfg.funparameter)
        hasfun = 1;
        
        st.nmt.pos = functional.pos;
        
        switch(cfg.funparameter)
            case {'mom'}
                 fun = cell2mat(getsubfield(functional,cfg.funparameter));
                
                inside = find(functional.inside);
                st.nmt.fun = nan(length(functional.inside),size(fun,2));
                st.nmt.fun(inside,:) = fun;
                
                st.nmt.time = functional.time;
                
                
                
                if(~isfield(cfg,'time') & ~isfield(cfg,'vox'))
                    [~,peakind] = max(abs(st.nmt.fun(:)));
                    [peakvoxind,peaktimeind] = ind2sub(size(st.nmt.fun),peakind);
                    cfg.time(1) = peaktimeind;
                    cfg.vox = peakvoxind;
                end

                if(~isfield(cfg,'time') & isfield(cfg,'vox'))
                    [~,peaktimeind] = max(abs(st.nmt.fun(cfg.vox,:)));
                    cfg.time(1) = peaktimeind;
                end

                if(isfield(cfg,'time') & ~isfield(cfg,'vox'))
                    [~,peakvoxind] = max(abs(st.nmt.fun(cfg.time,:)));
                    cfg.vox = peakvoxind;
                end
                
                % move MRI crosshairs to desired/peak voxel
                spm_orthviews('Reposition',st.nmt.pos(cfg.vox,:))
                
                if(length(cfg.time(1)) == 1)
                    cfg.time(2) = cfg.time(1);
                end
            otherwise
                st.nmt.fun = getsubfield(functional, cfg.funparameter);
                if isfield(functional, 'time')
                  st.nmt.time = functional.time;
                end
                if isfield(functional, 'freq')
                  st.nmt.freq = functional.freq;
                end
        end
        
        st.nmt.cfg = cfg;
    else
        error('cfg.funparameter not found in functional');
    end
else
    hasfun = 0;
    fprintf('no functional parameter\n');
end

nmt_spm8_plot(cfg)
nmt_image