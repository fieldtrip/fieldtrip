function fiducial=plot_ortho(data)
%
% PLOT_ORTHO visualizes an MRI 3D volume in 3 orthogonal projections
%
% Use as
%   hs = plot_ortho(data, varargin)
%
% Some defaults for the additional arguments:
%
%   'location'            = location of cut, (default = 'auto')
%                              'auto', 'center' if only anatomy, 'max' if functional data
%                              'min' and 'max' position of min/max funparameter
%                              'center' of the brain
%                               [x y z], coordinates in voxels or head, see locationcoordinates
%   'locationcoordinates' = coordinate system used in location, 'head' or 'voxel' (default = 'head')
%                              'head', headcoordinates from anatomical MRI
%                              'voxel', voxelcoordinates
%   'crosshair'           = 'yes' or 'no' (default = 'yes')
%   'axis'                = 'on' or 'off' (default = 'on')
%   'interactive'         = 'yes' or 'no' (default = 'no')
%                               in interactive mode cursor click determines location of cut
%   'queryrange'          = number, in atlas voxels (default 3)
%   'funcolorlim'         = color range of the functional data (default = 'auto')
%                               [min max]
%                               'maxabs', from -max(abs(funparameter)) to +max(abs(funparameter))
%                               'zeromax', from 0 to max(abs(funparameter))
%                               'minzero', from min(abs(funparameter)) to 0
%                               'auto', if funparameter values are all positive: 'zeromax',
%                               all negative: 'minzero', both possitive and negative: 'maxabs'
%   'funparameter'        = string, field in data with the functional parameter of interest (default = [])
%   'inputcoord'          = 'mni' or 'tal', coordinate system of data used to lookup the label from the atlas
%   'colorbar'            = 'yes' or 'no' (default = 'yes')
%
% Example
%   figure, plot_ortho(data,'colorbar','no','interactive','yes','axis','off')
%
% Copyright (C) 2009, Cristiano Micheli 
%
% $Log: plot_ortho.m,v $
% Revision 1.3  2009/04/21 13:24:13  crimic
% modified help
%
% Revision 1.2  2009/04/20 11:22:51  crimic
% first implementation
%

% get the optional input arguments
location            = keyval('location',  varargin); if isempty(location),location='auto';end
locationcoordinates = keyval('locationcoordinates',  varargin);  if isempty(locationcoordinates),locationcoordinates='head';end
crosshair1          = keyval('crosshair',  varargin);  if isempty(crosshair1),crosshair1='yes';end
axis1               = keyval('axis',  varargin);  if isempty(axis1),axis1='yes';end
interactive         = keyval('interactive',  varargin);  if isempty(interactive),interactive='no';end
queryrange          = keyval('queryrange',  varargin);  if isempty(queryrange),queryrange=3;end
funcolorlim         = keyval('funcolorlim',  varargin);  if isempty(funcolorlim),funcolorlim='auto';end
funparameter        = keyval('funparameter',  varargin);  if isempty(funparameter),funparameter=[];end
colorbar1           = keyval('colorbar',  varargin);  if isempty(colorbar1),colorbar1='yes';end

% initialize empty output structure
fiducial.nas = [];
fiducial.lpa = [];
fiducial.rpa = [];

% check if it is a suitable volumetric dataset
if isfield(data,'transform')
  if ~isstr(location)
    if strcmp(locationcoordinates, 'head')
      % convert the headcoordinates location into voxel coordinates
      loc = inv(data.transform) * [location(:); 1];
      loc = round(loc(1:3));
    elseif strcmp(locationcoordinates, 'voxel')
      % the location is already in voxel coordinates
      loc = round(location(1:3));
    else
      error('you should specify locationcoordinates');
    end
  else
    if isequal(location,'auto')
      if hasfun
        if isequal(funcolorlim,'maxabs');
          loc = 'max';
        elseif isequal(funcolorlim, 'zeromax');
          loc = 'max';
        elseif isequal(funcolorlim, 'minzero');
          loc = 'min';
        else %if numerical
          loc = 'max';
        end
      else
        loc = 'center';
      end;
    else
      loc = location;
    end
  end

  % determine the initial intersection of the cursor (xi yi zi)
  if isstr(loc) && strcmp(loc, 'min')
    if isempty(funparameter)
      error('location is min, but no functional parameter specified');
    end
    [minval, minindx] = min(fun(:));
    [xi, yi, zi] = ind2sub(dim, minindx);
  elseif isstr(loc) && strcmp(loc, 'max')
    if isempty(funparameter)
      error('location is max, but no functional parameter specified');
    end
    [maxval, maxindx] = max(fun(:));
    [xi, yi, zi] = ind2sub(dim, maxindx);
  elseif isstr(loc) && strcmp(loc, 'center')
    xi = round(dim(1)/2);
    yi = round(dim(2)/2);
    zi = round(dim(3)/2);
  elseif ~isstr(loc)
    % using nearest instead of round ensures that the position remains within the volume
    xi = nearest(1:dim(1), loc(1));
    yi = nearest(1:dim(2), loc(2));
    zi = nearest(1:dim(3), loc(3));
  end

  % % do the actual plotting %%
  nas = [];
  lpa = [];
  rpa = [];
  interactive_flag = 1; % it happens at least once
  while(interactive_flag)
    interactive_flag = strcmp(interactive, 'yes');

    xi = round(xi); xi = max(xi, 1); xi = min(xi, dim(1));
    yi = round(yi); yi = max(yi, 1); yi = min(yi, dim(2));
    zi = round(zi); zi = max(zi, 1); zi = min(zi, dim(3));

    if interactive_flag
      fprintf('\n');
      fprintf('click with mouse button to reposition the cursor\n');
      fprintf('press n/l/r on keyboard to record a fiducial position\n');
      fprintf('press q on keyboard to quit interactive mode\n');
    end

    ijk = [xi yi zi 1]';
    xyz = data.transform * ijk;
    if hasfun && ~hasatlas
      val = fun(xi, yi, zi);
      fprintf('voxel %d, indices [%d %d %d], location [%.1f %.1f %.1f], value %f\n', sub2ind(dim, xi, yi, zi), ijk(1:3), xyz(1:3), val);
    elseif hasfun && hasatlas
      val = fun(xi, yi, zi);
      fprintf('voxel %d, indices [%d %d %d], %s coordinates [%.1f %.1f %.1f], value %f\n', sub2ind(dim, xi, yi, zi), ijk(1:3), inputcoord, xyz(1:3), val);
    elseif ~hasfun && ~hasatlas
      fprintf('voxel %d, indices [%d %d %d], location [%.1f %.1f %.1f]\n', sub2ind(dim, xi, yi, zi), ijk(1:3), xyz(1:3));
    elseif ~hasfun && hasatlas
      fprintf('voxel %d, indices [%d %d %d], %s coordinates [%.1f %.1f %.1f]\n', sub2ind(dim, xi, yi, zi), ijk(1:3), inputcoord, xyz(1:3));
    end

    if hasatlas
      % determine the anatomical label of the current position 
      lab = atlas_lookup(atlas, (xyz(1:3)), 'inputcoord', inputcoord, 'queryrange', queryrange);
      if isempty(lab)
        fprintf([f,' labels: not found\n']);
      else
        fprintf([f,' labels: '])
        fprintf('%s', lab{1});
        for i=2:length(lab)
          fprintf(', %s', lab{i});
        end
        fprintf('\n');
      end
    end

    % make vols and scales, containes volumes to be plotted (fun, ana, msk)
    vols = {};
    if hasana; vols{1} = ana; scales{1} = []; end; % needed when only plotting ana
    if hasfun; vols{2} = fun; scales{2} = [fcolmin fcolmax]; end;
    if hasmsk; vols{3} = msk; scales{3} = [opacmin opacmax]; end;

    if isempty(vols)
      % this seems to be a problem that people often have
      error('no anatomy is present and no functional data is selected, please check your funparameter');
    end

    h1 = subplot(2,2,1);
    [vols2D] = handle_ortho(vols, [xi yi zi], 2, dim);
    plot2D(vols2D, scales);
    xlabel('i'); ylabel('k'); axis(axis1);
    if strcmp(crosshair1, 'yes'), crosshair([xi zi]); end
    
    h2 = subplot(2,2,2);
    [vols2D] = handle_ortho(vols, [xi yi zi], 1, dim);
    plot2D(vols2D, scales);
    xlabel('j'); ylabel('k'); axis(axis1);
    if strcmp(crosshair1, 'yes'), crosshair([yi zi]); end

    h3 = subplot(2,2,3);
    [vols2D] = handle_ortho(vols, [xi yi zi], 3, dim);
    plot2D(vols2D, scales);
    xlabel('i'); ylabel('j'); axis(axis1);
    if strcmp(crosshair1, 'yes'), crosshair([xi yi]); end

    if strcmp(colorbar1,  'yes'),
      if hasfun
        % vectorcolorbar = linspace(fcolmin, fcolmax,length(funcolormap));
        % imagesc(vectorcolorbar,1,vectorcolorbar);colormap(funcolormap);
        subplot(2,2,4);
        % use a normal Matlab coorbar, attach it to the invisible 4th subplot
        caxis([fcolmin fcolmax]);
        hc = colorbar;
        set(hc, 'YLim', [fcolmin fcolmax]);
        set(gca, 'Visible', 'off');
      else
        warning('no colorbar possible without functional data')
      end
    end

    drawnow;

    if interactive_flag
      try, [d1, d2, key] = ginput(1); catch, key='q'; end
      if isempty(key)
        % this happens if you press the apple key
        % do nothing
      elseif key=='q'
        break;
      elseif key=='l'
        lpa = [xi yi zi];
      elseif key=='r'
        rpa = [xi yi zi];
      elseif key=='n'
        nas = [xi yi zi];			
      elseif key=='i' || key=='j' || key=='k' || key=='m'
        % update the view to a new position
        if     l1=='i' && l2=='k' && key=='i', zi = zi+1; 
        elseif l1=='i' && l2=='k' && key=='j', xi = xi-1;
        elseif l1=='i' && l2=='k' && key=='k', xi = xi+1;
        elseif l1=='i' && l2=='k' && key=='m', zi = zi-1;
        elseif l1=='i' && l2=='j' && key=='i', yi = yi+1; 
        elseif l1=='i' && l2=='j' && key=='j', xi = xi-1;
        elseif l1=='i' && l2=='j' && key=='k', xi = xi+1;
        elseif l1=='i' && l2=='j' && key=='m', yi = yi-1;
        elseif l1=='j' && l2=='k' && key=='i', zi = zi+1; 
        elseif l1=='j' && l2=='k' && key=='j', yi = yi-1;
        elseif l1=='j' && l2=='k' && key=='k', yi = yi+1;
        elseif l1=='j' && l2=='k' && key=='m', zi = zi-1;
	end;
      else
        % update the view to a new position
        l1 = get(get(gca, 'xlabel'), 'string');
        l2 = get(get(gca, 'ylabel'), 'string');
        switch l1,
          case 'i'
            xi = d1;
          case 'j'
            yi = d1;
          case 'k'
            zi = d1;
        end
        switch l2,
          case 'i'
            xi = d2;
          case 'j'
            yi = d2;
          case 'k'
            zi = d2;
        end
      end
    end % if interactive_flag
    if ~isempty(nas), fprintf('nas = [%f %f %f]\n', nas); fiducial.nas = nas; else fprintf('nas = undefined\n'); end
    if ~isempty(lpa), fprintf('lpa = [%f %f %f]\n', lpa); fiducial.lpa = lpa; else fprintf('lpa = undefined\n'); end
    if ~isempty(rpa), fprintf('rpa = [%f %f %f]\n', rpa); fiducial.rpa = rpa; else fprintf('rpa = undefined\n'); end	    
  end % while interactive_flag
end

function plot2D(vols2D, scales);
cla;
% put 2D volumes in fun, ana and msk
hasana = length(vols2D)>0 && ~isempty(vols2D{1});
hasfun = length(vols2D)>1 && ~isempty(vols2D{2});
hasmsk = length(vols2D)>2 && ~isempty(vols2D{3});

% the transpose is needed for displaying the matrix using the Matlab image() function
if hasana; ana = vols2D{1}'; end;
if hasfun; fun = vols2D{2}'; end;
if hasmsk; msk = vols2D{3}'; end;


if hasana
  % scale anatomy between 0 and 1
  fprintf('scaling anatomy\n');
  amin = min(ana(:));
  amax = max(ana(:));
  ana = (ana-amin)./(amax-amin);
  clear amin amax;
  % convert anatomy into RGB values
  ana = cat(3, ana, ana, ana);
  ha = imagesc(ana);
end
hold on

if hasfun
  hf = imagesc(fun);
  caxis(scales{2});
  % apply the opacity mask to the functional data
  if hasmsk
    % set the opacity
    set(hf, 'AlphaData', msk)
    set(hf, 'AlphaDataMapping', 'scaled')
    alim(scales{3});
  elseif hasana
    set(hf, 'AlphaData', 0.5)
  end
end

axis equal
axis tight
axis xy

function [vols2D] = handle_ortho(vols, indx, slicedir, dim);

% put 2Dvolumes in fun, ana and msk
if length(vols)>=1 && isempty(vols{1}); hasana=0; else ana=vols{1}; hasana=1; end;
if length(vols)>=2
  if isempty(vols{2}); hasfun=0; else fun=vols{2}; hasfun=1; end;
else hasfun=0; end
if length(vols)>=3
  if isempty(vols{3}); hasmsk=0; else msk=vols{3}; hasmsk=1; end;
else hasmsk=0; end

% select the indices of the intersection
xi = indx(1);
yi = indx(2);
zi = indx(3);

% select the slice to plot
if slicedir==1
  yi = 1:dim(2);
  zi = 1:dim(3);
elseif slicedir==2
  xi = 1:dim(1);
  zi = 1:dim(3);
elseif slicedir==3
  xi = 1:dim(1);
  yi = 1:dim(2);
end

% cut out the slice of interest
if hasana; ana = squeeze(ana(xi,yi,zi)); end;
if hasfun; fun = squeeze(fun(xi,yi,zi)); end;
if hasmsk; msk = squeeze(msk(xi,yi,zi)); end;

%put fun, ana and msk in vols2D
if hasana; vols2D{1} = ana; end;
if hasfun; vols2D{2} = fun; end;
if hasmsk; vols2D{3} = msk; end;
