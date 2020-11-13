classdef ft_plot_mesh_interactive<handle
  % FT_PLOT_MESH_INTERACTIVE is the lower-level implementation class
  % corresponding to FT_SOURCEPLOT_INTERACTIVE. See that higher-level
  % function for details on configuration options.
  %
  % Note that this class deliberately does not use plotting helper functions
  % like ft_plot_mesh etc., for performance reasons.
  %
  % Copyright (C) 2019 Eelke Spaak, Donders Institute, e.spaak@donders.ru.nl
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  properties
    % data properties
    tri; pos; time; data; unit; 
    
    % configuration options
    time_label; pow_label; data_labels; has_diff; clim;
    
    % not yet configurable
    timeplot_colourmap; colourmap;
    
    % graphics handles
    fig_surface; axes_surface; surfs_surface; camlight_surface;
    figs_time; axes_time; vlines_time;
    virt_elec_surfs;
    
    % currently displayed time point
    cur_time;
    
    % number of conditions/positions/time points
    ncond; npos; ntim;
    
    % atlas for looking up locations of virtual electrodes
    atlas;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function self = ft_plot_mesh_interactive(varargin)
      % Constructor
      
      self.tri = ft_getopt(varargin, 'tri');
      self.pos = ft_getopt(varargin, 'pos');
      self.time = ft_getopt(varargin, 'time');
      self.data = ft_getopt(varargin, 'data');
      self.unit = ft_getopt(varargin, 'unit');
      
      % configuration options
      self.time_label = ft_getopt(varargin, 'time_label', 'Time (s)');
      self.pow_label = ft_getopt(varargin, 'pow_label', 'Current density (a.u.)');
      self.data_labels = ft_getopt(varargin, 'data_labels',...
        arrayfun(@(x) sprintf('Input %d', x), 1:numel(self.data), 'uniformoutput', false));
      self.clim = ft_getopt(varargin, 'clim', []);
      % has_diff: treat the last input argument as special, and assign
      % different colour limits for the corresponding surface plot
      self.has_diff = ft_getopt(varargin, 'has_diff', false);
      self.atlas    = ft_getopt(varargin, 'atlas');
      self.colourmap = ft_getopt(varargin, 'colormap');
      
      if isempty(self.clim)
        if self.has_diff
          self.clim = [0 max(cellfun(@(x) max(x(:)), self.data(1:end-1)))*0.75];
        else
          self.clim = [0 max(cellfun(@(x) max(x(:)), self.data))*0.75];
        end
      end
      
      % we need brewermap
      ft_hastoolbox('brewermap', 1);
      if isempty(self.colourmap)
        self.colourmap = brewermap(64, 'YlOrRd');
      end
      self.axes_surface = [];
      self.surfs_surface = [];
      self.figs_time = [];
      
      self.virt_elec_surfs = [];
      
      self.ncond = numel(self.data);
      self.npos = size(self.pos, 1);
      self.ntim = numel(self.time);
      
      % an index into the "current" time point (in seconds)
      self.cur_time = 0;
      
      
      
      self.timeplot_colourmap = brewermap(8, 'Set2');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function show(self)
      % Show the interactive source viewer.
      self.init_surface_plots();
      self.add_time_plot(true(self.npos,1));
      title(self.axes_time(1), 'Average over cortex');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function init_surface_plots(self)
      % INIT_SURFACE_PLOTS initializes a single figure that displays
      % surface plots for all the functional data at a single time point.
      self.fig_surface = figure('Color', 'w');
      ft_colormap(self.fig_surface, self.colourmap);
      for k = 1:self.ncond
        % initialize axes and surface
        self.axes_surface(k) = subtightplot(1, numel(self.data), k);        
        
        hold(self.axes_surface(k), 'on');
        self.surfs_surface(k) = patch('Vertices', self.pos, 'Faces', self.tri, 'EdgeColor', 'None',...
          'FaceColor', 'interp', 'CDataMapping', 'scaled', 'FaceVertexCData',...
          self.data{k}(:,self.tim_to_ind(self.cur_time)));
        
        % set display properties
        lighting gouraud;
        shading interp;
        material shiny;
        axis vis3d equal off;
        set(self.axes_surface(k),'CameraViewAngleMode','Manual');
        set(self.axes_surface(k), 'CLim', self.clim);
        self.camlight_surface(k) = camlight('infinite');
        
        % treat an optional condition difference separately
        if self.has_diff && k == self.ncond
          tmp = self.data{k};
          lims = max(tmp(:)) * 0.75;
          set(self.axes_surface(k), 'CLim', [-lims lims]);
          ft_colormap(self.axes_surface(k), brewermap(64, '*RdBu'));
        end
        
        title(self.axes_surface(k), self.data_labels{k});
        
        self.init_surface_click_event(self.axes_surface(k), self.surfs_surface(k));
      end
      
      % link the axes
      Link = linkprop(self.axes_surface, {'CameraUpVector', 'CameraPosition',...
        'CameraTarget', 'XLim', 'YLim', 'ZLim'});
      setappdata(self.fig_surface, 'StoreTheLink', Link);
      
      % automatically reposition camlights when rotation has ended
      function rot_ended(~,~)
        for k = 1:numel(self.axes_surface)
          camlight(self.camlight_surface(k), 'infinite');
        end
      end
      h = rotate3d(self.fig_surface);
      h.ActionPostCallback = @rot_ended;
      h.Enable = 'on';
      % prevent the rotation manager from capturing events when shift is
      % down, because we want Shift+Click to add a virtual electrode (see
      % init_surface_click_event)
      function prevent = rotfilter(~,~)
        mod = get(self.fig_surface, 'CurrentModifier');
        prevent = ismember('shift', mod);
      end
      h.ButtonDownFilter = @rotfilter;
      
      % add button to set colour limits
      limit_button = uicontrol('Parent', self.fig_surface, 'Style', 'pushbutton',...
        'String', 'Rescale colours', 'Units', 'normalized', 'Position', [0.01, 0.01, 0.1, 0.05],...
        'Visible', 'on');
      function set_clim(~,~)
        timind = self.tim_to_ind(self.cur_time);
        if self.has_diff
          clim_conds = [0 max(cellfun(@(x) max(x(:,timind)), self.data(1:end-1)))*0.75];
          clim_diff = max(abs(self.data{end}(:,timind))*0.7);
          clim_diff = [-clim_diff clim_diff];
          for k = 1:self.ncond-1
            set(self.axes_surface(k), 'CLim', clim_conds);
          end
          set(self.axes_surface(end), 'Clim', clim_diff);
        else
          clim_conds = [0 max(cellfun(@(x) max(x(:,timind)), self.data(1:end-1)))*0.75];
          for k = 1:self.ncond
            set(self.axes_surface(k), 'CLim', clim_conds);
          end
        end
      end
      limit_button.Callback = @set_clim;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function init_surface_click_event(self, ax, surface)
      % INIT_SURFACE_CLICK_EVENT adds an on-click handler that adds a
      % "virtual electrode" to all surface plots, and opens a new figure
      % in which the time course of the functional data at that location is
      % plotted.
      function cb_plot_timecourse(~,~)
        point = get(ax, 'CurrentPoint');
        % get the intersection with the mesh (largely taken from
        % ft_sourcemovie)
        [ipos, d] = intersect_line(self.pos, self.tri, point(1,:), point(2,:));
        [~, ix]  = min(abs(d));
        if ~isempty(ix)
          dpos = self.pos - ipos(ix*ones(size(self.pos,1),1),:);
          index = nearest(sum(dpos.^2,2),0);
          newind = numel(self.figs_time) + 1;
          colour = self.timeplot_colourmap(mod(newind, 8)+1,:);
          
          thisfig_ind = self.add_time_plot(index);
          
          % render a little sphere at the location of the virtual electrode
          [x,y,z] = sphere();
          s.pos  = [x(:) y(:) z(:)].*4;
          s.unit = 'mm';
          s      = ft_convert_units(s, self.unit);
          x      = reshape(s.pos(:,1),sqrt(size(s.pos,1)).*[1 1]);
          y      = reshape(s.pos(:,2),sqrt(size(s.pos,1)).*[1 1]);
          z      = reshape(s.pos(:,3),sqrt(size(s.pos,1)).*[1 1]);
          for k = 1:self.ncond
            self.virt_elec_surfs(thisfig_ind,k) = surf(...
             self.axes_surface(k), x+self.pos(index,1), y+self.pos(index,2), z+self.pos(index,3),...
             'facecolor', colour, 'edgecolor', 'none', 'facelighting', 'gouraud');
            material shiny;
          end
          
          % set the new axis background color to reflect the color of the
          % virtual electrode
          set(self.axes_time(end), 'Color', [colour 0.1]);
          
          if ~isempty(self.atlas)
            if isfield(self.atlas, 'dim') && isfield(self.atlas, 'transform')
              % use atlas_lookup, which operates on volumetrically defined
              % atlases
              atlas_labels = atlas_lookup(self.atlas, self.pos(index,:), 'coordsys', 'mni');
            elseif isfield(self.atlas, 'pos') && size(self.atlas.pos,1)==size(self.pos,1)
              % assume that the topologies are the same, and use the index
              % of the selected vertex to look up the label.
              fn = fieldnames(self.atlas);
              atlas_labels = {};
              for i=1:numel(fn)
                if any(strcmp(fn, [fn{i} 'label']))
                  atlas_labels{end+1} = self.atlas.([fn{i} 'label']){self.atlas.(fn{i})(index)};
                end
              end
            else
              atlas_labels = [];
            end
            if numel(atlas_labels) > 1
              atlas_labels = unique(atlas_labels);
              tmp = sprintf('%s', strrep(atlas_labels{1}, '_', ' '));
              for i=2:length(atlas_labels)
                tmp = [tmp sprintf(', %s', strrep(atlas_labels{i}, '_', ' '))];
              end
              atlas_labels = tmp;
            else
              atlas_labels = '(no atlas label)';
            end
          else
            atlas_labels = '(no atlas available)';
          end
          
          title(self.axes_time(thisfig_ind), sprintf('At position (%.0f, %.0f, %.0f) mm\n%s',...
            self.pos(index,1), self.pos(index,2), self.pos(index,3), atlas_labels));
        end
      end
      set(surface, 'ButtonDownFcn', @cb_plot_timecourse);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function thisfig_ind = add_time_plot(self, posinds)
      % ADD_TIME_PLOT creates a new figure that shows the functional data
      % over time at the specified positional indices. If multiple indices
      % are provided, the average value over those positions is plotted.
      
      % set up figure and store handles
      fig = figure('Color', 'w');
      thisfig_ind = numel(self.figs_time) + 1;
      self.figs_time(thisfig_ind) = fig;
      ax = axes();
      self.axes_time(thisfig_ind) = ax;

      % get data to plot, ignoring a possible last special input argument,
      % and averaging over grid points
      plotdata = zeros(self.ntim, self.ncond);
      for k = 1:self.ncond-self.has_diff
        plotdata(:,k) = nanmean(self.data{k}(posinds,:), 1)';
      end
      
      pl = plot(ax, self.time, plotdata);
      
      % set up some appearance properties
      legend(self.data_labels(1:end-self.has_diff));
      xlabel(self.time_label);
      ylabel(self.pow_label);
      grid on;
      box off;
      
      % add a vertical line indicating the current time point we're
      % plotting in the surface plots
      hold(ax, 'on');
      self.vlines_time(thisfig_ind) = plot([self.cur_time self.cur_time], get(ax, 'ylim'),...
        'k:', 'Color', [0 0 0 0.5], 'LineWidth', 2, 'HandleVisibility', 'off');
      
      % add event handler: clicking somewhere should move the indicator and
      % update the surface plots
      is_down = false;
      function button_down(~,~)
        is_down = true;
        point = get(ax, 'CurrentPoint');
        self.fire_timepoint_change(point(1));
      end
      function button_up(~,~), is_down = false; end
      function mouse_dragged(~,~)
        if is_down
          point = get(ax, 'CurrentPoint');
          self.fire_timepoint_change(point(1));
        end
      end
      set(ax, 'ButtonDownFcn', @button_down);
      set(fig, 'WindowButtonUpFcn', @button_up);
      set(fig, 'WindowButtonMotionFcn', @mouse_dragged);
      % prevent the lines from capturing mouse events
      set(pl, 'PickableParts', 'none');
      
      % delete the handles properly when the user closes the figure
      function window_closed(~,~)
        try
          self.figs_time(thisfig_ind) = [];
          self.axes_time(thisfig_ind) = [];
          self.vlines_time(thisfig_ind) = [];
          delete(self.virt_elec_surfs(thisfig_ind,:));
          self.virt_elec_surfs(thisfig_ind,:) = [];
          delete(fig);
        catch
          delete(fig);
        end
      end
      set(fig, 'CloseRequestFcn', @window_closed);
      
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function fire_timepoint_change(self, t)
      % FIRE_TIMEPOINT_CHANGE handles the event triggered when the active
      % time point (should) change(s). Surface plots colors are updated,
      % and the vertical time point reference line in any time axis plots.
      self.cur_time = t;
      for k = 1:self.ncond
        set(self.surfs_surface(k), 'FaceVertexCData',...
          self.data{k}(:,self.tim_to_ind(t)));
      end
      for k = 1:numel(self.figs_time)
          set(self.vlines_time(k), 'XData', [t t]);
      end
      set(self.fig_surface, 'Name', sprintf('t = %.2f', t));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function ind = tim_to_ind(self, tim)
      % TIM_TO_IND takes a time point and returns the corresponding index
      % into the time axis.
      [~,ind] = min(abs(self.time-tim));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  end
  
end
