classdef ft_plot_mesh_interactive<handle
  % FT_PLOT_MESH_INTERACTIVE is the lower-level implementation class
  % corresponding to FT_SOURCEPLOT_INTERACTIVE. See that higher-level
  % function for details on configuration options.
  %
  % Note that this class deliberately does not use plotting helper functions
  % like ft_plot_mesh etc., for performance reasons.
  %
  % Copyright (C) 2019 Eelke Spaak, Donders Institute, e.spaak@donders.ru.nl
  % Copyright (C) 2024 Jan-Mathijs Schoffelen, Donders Institute
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  properties
    % data properties
    tri; pos; time; freq; data; unit; 
    
    % configuration options
    time_label; freq_label; pow_label; data_labels; has_diff; clim; clims
    
    % not yet configurable
    timeplot_colourmap; colourmap;
    
    % graphics handles
    fig_surface; axes_surface; surfs_surface; camlight_surface;
    figs_time; axes_time; vlines_time; rect_time;
    virt_elec_surfs;
    
    % currently displayed time point
    cur_time; t1; t2;
    cur_freq; f1; f2;

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
      self.freq = ft_getopt(varargin, 'freq');
      self.data = ft_getopt(varargin, 'data');
      self.unit = ft_getopt(varargin, 'unit');
      
      % configuration options
      self.time_label = ft_getopt(varargin, 'time_label', 'Time (s)');
      if ~isempty(self.freq) && numel(self.freq)>1
        self.freq_label = ft_getopt(varargin, 'freq_label', 'Frequency (Hz)');
      end
      self.pow_label = ft_getopt(varargin, 'pow_label', 'Current density (a.u.)');
      self.data_labels = ft_getopt(varargin, 'data_labels',...
        arrayfun(@(x) sprintf('Input %d', x), 1:numel(self.data), 'uniformoutput', false));
      self.clim = ft_getopt(varargin, 'clim', 'zeromax');
      % has_diff: treat the last input argument as special, and assign
      % different colour limits for the corresponding surface plot
      self.has_diff = ft_getopt(varargin, 'has_diff', false);
      self.atlas    = ft_getopt(varargin, 'atlas');
      self.colourmap = ft_getopt(varargin, 'colormap');
      
      n = numel(self.data);
      if self.has_diff
        n = numel(self.data)-1;
      end
      sat = 1;
      if ischar(self.clim) && strcmp(self.clim, 'maxmin')
        self.clims = [min(cellfun(@(x) min(x(:)), self.data(1:n)))*sat max(cellfun(@(x) max(x(:)), self.data(1:n)))*sat];
      elseif ischar(self.clim) && strcmp(self.clim, 'maxabs')
        tmpclim = [min(cellfun(@(x) min(x(:)), self.data(1:n)))*sat max(cellfun(@(x) max(x(:)), self.data(1:n)))*sat];
        self.clims = [-1 1].*max(abs(tmpclim));
      elseif ischar(self.clim) && strcmp(self.clim, 'zeromax')
        self.clims = [0 max(cellfun(@(x) max(x(:)), self.data(1:n)))*sat];
      elseif ischar(self.clim) && strcmp(self.clim, 'minzero')
        self.clims = [min(cellfun(@(x) min(x(:)), self.data(1:n)))*sat 0];
      else
        self.clims = self.clim;
      end
      
      % we need brewermap
      ft_hastoolbox('brewermap', 1);
      if isempty(self.colourmap)
        self.colourmap = brewermap(64, 'YlOrRd');
      end
      self.axes_surface  = [];
      self.surfs_surface = [];
      self.figs_time     = [];
      
      self.virt_elec_surfs = [];
      
      self.ncond = numel(self.data);
      self.npos  = size(self.pos, 1);
      self.ntim  = numel(self.time);
      
      % an index into the "current" time point (in seconds)
      self.cur_time = 0;
      self.t1 = 0;
      self.t2 = 0;
      
      if ~isempty(self.freq)
        self.cur_freq = mean(self.freq);
        self.f1 = self.freq(1);
        self.f2 = self.freq(end);
      end
 
      self.timeplot_colourmap = brewermap(8, 'Set2');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function show(self)
      % Show the interactive source viewer.
      self.init_surface_plots();
      if ~isempty(self.freq)
        self.add_tf_plot(true(self.npos,1));
      else
        self.add_time_plot(true(self.npos,1));
      end
      title(self.axes_time(1), 'Average over cortex');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function init_surface_plots(self)
      % INIT_SURFACE_PLOTS initializes a single figure that displays
      % surface plots for all the functional data at a single time point.
      self.fig_surface = figure('Color', 'w');
      p = get(self.fig_surface, 'position');
      p(1) = 10; % move the figure to the left
      set(self.fig_surface, 'position', p);

      ft_colormap(self.fig_surface, self.colourmap);
      for k = 1:self.ncond
        % initialize axes and surface
        self.axes_surface(k) = subtightplot(1, numel(self.data), k);        
        
        hold(self.axes_surface(k), 'on');
        if isempty(self.freq)
          datavector = self.data{k}(:,self.tim_to_ind(self.cur_time));
        else
          datavector = self.data{k}(:,self.frq_to_ind(self.cur_freq), self.tim_to_ind(self.cur_time));
        end
        self.surfs_surface(k) = patch('Vertices', self.pos, 'Faces', self.tri, 'EdgeColor', 'None',...
          'FaceColor', 'interp', 'CDataMapping', 'scaled', 'FaceVertexCData',...
          datavector);
        
        % set display properties
        lighting gouraud;
        shading interp;
        material dull;
        axis vis3d equal off;
        set(self.axes_surface(k),'CameraViewAngleMode','Manual');
        set(self.axes_surface(k), 'CLim', self.clims);
        self.camlight_surface(k) = camlight('infinite');
        
        % treat an optional condition difference separately
        if self.has_diff && k == self.ncond
          sat = 1;
          tmp = self.data{k};
          lims = max(tmp(:)) * sat;
          set(self.axes_surface(k), 'CLim', [-lims lims]);
          ft_colormap(self.axes_surface(k), brewermap(64, '*RdBu')); % FIXME: not sure why the colormap is hard-coded here
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
        
        if ~isempty(self.cur_freq)
          frqind = sort([self.frq_to_ind(self.f1) self.frq_to_ind(self.f2)]);
          timind = sort([self.tim_to_ind(self.t1) self.tim_to_ind(self.t2)]);
          ix = {':' frqind(1):frqind(2) timind(1):timind(2)};
        else
          timind = [self.tim_to_ind(self.t1) self.tim_to_ind(self.t2)];
          ix = {':' timind(1):timind(2)};
        end
        if self.has_diff
          n = numel(self.data)-1;
        else
          n = numel(self.data);
        end

        sat = 1;
        if (ischar(self.clim) && strcmp(self.clim, 'maxmin')) || isnumeric(self.clim)
          % FIXME the second condition is old default behavior I think (but
          % probably not optimal in its alignment with possible user expectations)
          clim_conds = [min(cellfun(@(x) min(nanmean(nanmean(x(ix{:}),3),2)), self.data(1:n)))*sat max(cellfun(@(x) max(nanmean(nanmean(x(ix{:}),3),2)), self.data(1:n)))*sat];
        elseif ischar(self.clim) && strcmp(self.clim, 'maxabs')
          tmpclim = [min(cellfun(@(x) min(nanmean(nanmean(x(ix{:}),3),2)), self.data(1:n)))*sat max(cellfun(@(x) max(nanmean(nanmean(x(ix{:}),3),2)), self.data(1:n)))*sat];
          clim_conds = [-1 1].*max(abs(tmpclim));
        elseif ischar(self.clim) && strcmp(self.clim, 'zeromax')
          clim_conds = [0 max(cellfun(@(x) max(nanmean(nanmean(x(ix{:}),3),2)), self.data(1:n)))*sat];
        elseif ischar(self.clim) && strcmp(self.clim, 'minzero')
          clim_conds = [min(cellfun(@(x) min(nanmean(nanmean(x(ix{:}),3),2)), self.data(1:n)))*sat 0];
        end
        for k = 1:n
          set(self.axes_surface(k), 'CLim', clim_conds);
        end

        if self.has_diff
          clim_diff = max(abs(self.data{end}(:,timind))*0.7);
          clim_diff = [-clim_diff clim_diff];
          set(self.axes_surface(end), 'Clim', clim_diff);
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
        [dum, ix]  = min(abs(d));
        if ~isempty(ix)
          dpos = self.pos - ipos(ix*ones(size(self.pos,1),1),:);
          index = nearest(sum(dpos.^2,2),0);
          newind = numel(self.figs_time) + 1;
          colour = self.timeplot_colourmap(mod(newind, 8)+1,:);
          
          if isempty(self.freq)
            thisfig_ind = self.add_time_plot(index);
          else
            thisfig_ind = self.add_tf_plot(index);
          end

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
                if any(strcmp(fn, [fn{i} 'label'])) && isfinite(self.atlas.(fn{i})(index))
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
      self.rect_time(thisfig_ind) = patch([self.t1 self.t1 self.t2 self.t2], [get(ax, 'ylim') flip(get(ax, 'ylim'),2)], 0.8.*[1 1 1],...
        'FaceAlpha', 0.3, 'Edgecolor', 'none');

      % add event handler: clicking somewhere should move the indicator and
      % update the surface plots
      is_down = false;
      is_up   = false;
      function button_down(~,~)
        is_down = true;
        is_up   = false;
        point   = get(ax, 'CurrentPoint');
        self.t1 = point(1);
        set(self.rect_time(thisfig_ind), 'XData', self.t1.*[1 1 1 1]);
      end
      function button_up(~,~)
        is_down = false; 
        is_up   = true;
        point   = get(ax, 'CurrentPoint');
        self.t2 = point(1);
        self.fire_timepoint_change([self.t1 self.t2]);
        set(self.vlines_time(thisfig_ind), 'Visible', 'on');
      end
      function mouse_dragged(~,~)
        if is_down
          point = get(ax, 'CurrentPoint');
          set(self.rect_time(thisfig_ind), 'XData', [self.t1 self.t1 point(1) point(1)]);
          set(self.vlines_time(thisfig_ind), 'Visible', 'off');
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
    
    function thisfig_ind = add_tf_plot(self, posinds)
      % ADD_TF_PLOT creates a new figure that shows the functional data
      % over time and frequencies at the specified positional indices. If multiple indices
      % are provided, the average value over those positions is plotted.
      
      % set up figure and store handles
      p = get(self.fig_surface, 'position');
      p(1) = p(1)+p(3)+10;
      fig = figure('Color', 'w', 'position', p);
      thisfig_ind = numel(self.figs_time) + 1;
      self.figs_time(thisfig_ind) = fig;
      ax = axes();
      self.axes_time(thisfig_ind) = ax;

      % get data to plot, FIXME this needs to gracefully deal with multiple
      % conditions
      if numel(self.data)>1
        error('only a single data input is supported for time frequency data');
      end
      plotdata = shiftdim(nanmean(self.data{1}(posinds,:,:), 1), 1);
      
      pl = imagesc(ax, self.time, self.freq, plotdata); axis xy
      if ischar(self.clim)
        switch self.clim
          case 'maxmin'
            cclim = [min(plotdata(:)) max(plotdata(:))];
          case 'maxabs'
            cclim = [-1 1].*max(abs(plotdata(:)));
          case 'zeromax'
            cclim = [0 max(plotdata(:))];
          case 'minzero'
            cclim = [min(plotdata(:)) 0];
        end
      else
        cclim = self.clim;
      end

      if all(isfinite(cclim))
        set(ax, 'CLim', cclim);
      end
      
      colorbar(ax);

      % set up some appearance properties
      xlabel(self.time_label);
      ylabel(self.freq_label);
      ft_colormap(self.colourmap);
      
      box off;
      
      self.rect_time(thisfig_ind) = patch([self.t1 self.t1 self.t2 self.t2], [self.f1 self.f2 self.f2 self.f1], 0.8.*[1 1 1],...
        'FaceAlpha', 0.3, 'Edgecolor', 'k');

      % add event handler: clicking somewhere should move the indicator and
      % update the surface plots
      is_down = false;
      is_up   = false;
      function button_down(~,~)
        is_down = true;
        is_up   = false;
        point   = get(ax, 'CurrentPoint');
        self.t1 = point(1);
        self.f1 = point(3);
        set(self.rect_time(thisfig_ind), 'XData', self.t1.*[1 1 1 1]);
        set(self.rect_time(thisfig_ind), 'YData', self.f1.*[1 1 1 1]);
      end
      function button_up(~,~)
        is_down = false; 
        is_up   = true;
        point   = get(ax, 'CurrentPoint');
        self.t2 = point(1);
        self.f2 = point(3);
        self.fire_tfpoint_change([self.f1 self.f2],[self.t1 self.t2]);
        %set(self.vlines_time(thisfig_ind), 'Visible', 'on');
      end
      function mouse_dragged(~,~)
        if is_down
          point = get(ax, 'CurrentPoint');
          set(self.rect_time(thisfig_ind), 'XData', [self.t1 self.t1 point(1) point(1)]);
          set(self.rect_time(thisfig_ind), 'YData', [self.f1 point(3) point(3) self.f1]);
          
          %set(self.vlines_time(thisfig_ind), 'Visible', 'off');
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
          %self.vlines_time(thisfig_ind) = [];
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
      self.cur_time = mean(t);
      for k = 1:self.ncond
        set(self.surfs_surface(k), 'FaceVertexCData',...
          mean(self.data{k}(:,self.tim_to_ind(t(1)):self.tim_to_ind(t(2))),2));
      end
      for k = 1:numel(self.figs_time)
          set(self.vlines_time(k), 'XData', [1 1].*mean(t));
      end
      set(self.fig_surface, 'Name', sprintf('t = %.2f - %.2f', t(1), t(2)));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function fire_tfpoint_change(self, f, t)
      % FIRE_TIMEPOINT_CHANGE handles the event triggered when the active
      % time point (should) change(s). Surface plots colors are updated,
      % and the vertical time point reference line in any time axis plots.
      t = sort(t);
      f = sort(f);
      self.cur_time = mean(t);
      self.cur_freq = mean(f);
      for k = 1:self.ncond
        set(self.surfs_surface(k), 'FaceVertexCData',...
          nanmean(nanmean(self.data{k}(:,self.frq_to_ind(f(1)):self.frq_to_ind(f(2)),self.tim_to_ind(t(1)):self.tim_to_ind(t(2))),2),3));
      end
      % for k = 1:numel(self.figs_time)
      %     set(self.vlines_time(k), 'XData', [1 1].*mean(t));
      % end
      set(self.fig_surface, 'Name', sprintf('t = %.2f - %.2f, f = %.2f - %.2f', t(1), t(2), f(1), f(2)));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function ind = tim_to_ind(self, tim)
      % TIM_TO_IND takes a time point and returns the corresponding index
      % into the time axis, or frequency axis, for spectra.
      [dum,ind] = min(abs(self.time-tim));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    function ind = frq_to_ind(self, frq)
      % FRQ_TO_IND takes a time point and returns the corresponding index
      % into the freq axis, for TF-data.
      [dum,ind] = min(abs(self.freq-frq));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  end
  
end
