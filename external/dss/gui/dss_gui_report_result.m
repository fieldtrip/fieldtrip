function dss_gui_report_result(S, state)
% Plots gathered data
% - components
% For each iteration (if data was gathered):
% - angles between consecutive projection vectors
% - angles between consecutive projection vector differences
% - value of the estimated objective function
%

% Copyright (C) 2004 Helsinki University of Technology
%   Laboratory of Computer and Information Science
%     Kosti Rytkönen (kosti@cis.hut.fi)
% $Id: dss_gui_report_result.m,v 1.1 2005/12/07 11:57:06 jaakkos Exp $

plot_count=sum([isfield(state.report, 'change') ...
		isfield(state.report, 'ddir') ...
		isfield(state.report, 'objective')]) + 1;
components = size(S, 1);
clf;

for component = 1:components
    
  plotindex=(component-1)*plot_count+1;

  subplot(components, plot_count, plotindex);
  plot(S(component, :));
  if component==1, title('Estimated signals'), end
  yscale = max(-min(S(component,:)), max(S(component,:)))*1.1;
  axis([0 size(S(component,:),2) -yscale yscale]);

  if ~isfield(state.report, 'iterations'); continue; end
  
  if (size(state.report.iterations) == [1 1])
    iterations = state.report.iterations;
  else
    iterations = state.report.iterations(component);
  end

  if iterations < 2; continue; end


  if isfield(state.report, 'change')
    plotindex = plotindex + 1;
    subplot(components, plot_count, plotindex);
    logplot(state.report.change(component, 1:iterations));
    %if component==1, title('Angle between consecutive projection vectors'), end
    if component==1, title('dw angle'), end
    xlabel('iteration');
    ylabel('angle');
    yscale = max(state.report.change(component,:));
    axis([1 iterations+1 -yscale*0.1 yscale*1.1]);
  end
  
    
  if isfield(state.report, 'ddir')
    plotindex = plotindex + 1;
    subplot(components, plot_count, plotindex);
    plot(state.report.ddir(component,2:iterations),'.');
    %if component==1, title('Angle between consecutive projection vector differences'), end
    if component==1, title('ddw angle'), end
    xlabel('iteration');
    ylabel('angle');
    axis([0 iterations+1 -10 190]);
  end

  if isfield(state.report, 'objective')
    plotindex = plotindex + 1;
    subplot(components, plot_count, plotindex);
    plot(state.report.objective(component,2:iterations));
    if component==1, title('Objective function'), end
    xlabel('iteration');
    ylabel('objective');
    %axis([0 iterations+1 -10 190]);
  end

end

