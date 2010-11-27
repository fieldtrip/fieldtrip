function icaplot(mode, varargin);
%ICAPLOT - plot signals in various ways
%
% ICAPLOT is mainly for plottinf and comparing the mixed signals and
% separated ica-signals.
%
% ICAPLOT has many different modes. The first parameter of the function
% defines the mode. Other parameters and their order depends on the
% mode. The explanation for the more common parameters is in the end.
%
% Classic
%     icaplot('classic', s1, n1, range, xrange, titlestr)
%
%     Plots the signals in the same manner as the FASTICA and FASTICAG
%     programs do. All the signals are plotted in their own axis.
%
% Complot
%     icaplot('complot', s1, n1, range, xrange, titlestr)
%
%     The signals are plotted on the same axis. This is good for
%     visualization of the shape of the signals. The scale of the signals 
%     has been altered so that they all fit nicely.
%
% Histogram
%     icaplot('histogram', s1, n1, range, bins, style)
%     
%     The histogram of the signals is plotted. The number of bins can be
%     specified with 'bins'-parameter. The style for the histograms can
%     be either 'bar' (default) of 'line'.
%
% Scatter
%     icaplot('scatter', s1, n1, s2, n2, range, titlestr, s1label,
%     s2label, markerstr)
%
%     A scatterplot is plotted so that the signal 1 is the 'X'-variable
%     and the signal 2 is the 'Y'-variable. The 'markerstr' can be used
%     to specify the maker used in the plot. The format for 'markerstr'
%     is the same as for Matlab's PLOT. 
%
% Compare
%     icaplot('compare', s1, n1, s2, n2, range, xrange, titlestr,
%     s1label, s2label)
%
%     This for for comparing two signals. The main used in this context
%     would probably be to see how well the separated ICA-signals explain 
%     the observed mixed signals. The s2 signals are first scaled with
%     REGRESS function.
%
% Compare - Sum
%     icaplot('sum', s1, n1, s2, n2, range, xrange, titlestr, s1label,
%     s2label)
%
%     The same as Compare, but this time the signals in s2 (specified by
%     n2) are summed together.
%
% Compare - Sumerror
%     icaplot('sumerror', s1, n1, s2, n2, range, xrange, titlestr,
%     s1label, s2label)
%     
%     The same as Compare - Sum, but also the 'error' between the signal
%     1 and the summed IC's is plotted.
%
%
% More common parameters
%     The signals to be plotted are in matrices s1 and s2. The n1 and n2
%     are used to tell the index of the signal or signals to be plotted
%     from s1 or s2. If n1 or n2 has a value of 0, then all the signals
%     from corresponding matrix will be plotted. The values for n1 and n2 
%     can also be vectors (like: [1 3 4]) In some casee if there are more
%     than 1 signal to be plotted from s1 or s2 then the plot will
%     contain as many subplots as are needed. 
%
%     The range of the signals to be plotted can be limited with
%     'range'-parameter. It's value is a vector ( 10000:15000 ). If range 
%     is 0, then the whole range will be plotted.
%
%     The 'xrange' is used to specify only the labels used on the
%     x-axis. The value of 'xrange' is a vector containing the x-values
%     for the plots or [start end] for begin and end of the range
%     ( 10000:15000 or [10 15] ). If xrange is 0, then value of range
%     will be used for x-labels.
%
%     You can give a title for the plot with 'titlestr'. Also the
%     's1label' and 's2label' are used to give more meaningfull label for 
%     the signals.
%
%     Lastly, you can omit some of the arguments from the and. You will
%     have to give values for the signal matrices (s1, s2) and the
%     indexes (n1, n2)

% @(#)$Id$

switch mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 'dispsig' is to replace the old DISPSIG
  % '' & 'classic' are just another names - '' quite short one :-)
 case {'', 'classic', 'dispsig'} 
  % icaplot(mode, s1, n1, range, xrange, titlestr)
  if length(varargin) < 1, error('Not enough arguments.'); end
  if length(varargin) < 5, titlestr = '';else titlestr = varargin{5}; end
  if length(varargin) < 4, xrange = 0;else xrange = varargin{4}; end
  if length(varargin) < 3, range = 0;else range = varargin{3}; end
  if length(varargin) < 2, n1 = 0;else n1 = varargin{2}; end
  s1 = varargin{1};
  range=chkrange(range, s1);
  xrange=chkxrange(xrange, range);
  n1=chkn(n1, s1);

  clf;
  
  numSignals = size(n1, 2);
  for i = 1:numSignals,
    subplot(numSignals, 1, i);
    plot(xrange, s1(n1(i), range));
  end
  subplot(numSignals,1, 1);
  if (~isempty(titlestr))
    title(titlestr);
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 case 'complot'
  % icaplot(mode, s1, n1, range, xrange, titlestr)
  if length(varargin) < 1, error('Not enough arguments.'); end
  if length(varargin) < 5, titlestr = '';else titlestr = varargin{5}; end
  if length(varargin) < 4, xrange = 0;else xrange = varargin{4}; end
  if length(varargin) < 3, range = 0;else range = varargin{3}; end
  if length(varargin) < 2, n1 = 0;else n1 = varargin{2}; end
  s1 = remmean(varargin{1});
  range=chkrange(range, s1);
  xrange=chkxrange(xrange, range);
  n1=chkn(n1, s1);
  
  for i = 1:size(n1, 2)
    S1(i, :) = s1(n1(i), range);
  end
  
  alpha = mean(max(S1')-min(S1'));
  for i = 1:size(n1,2)
    S2(i,:) = S1(i,:) - alpha*(i-1)*ones(size(S1(1,:)));
  end
  
  plot(xrange, S2');
  axis([min(xrange) max(xrange) min(min(S2)) max(max(S2)) ]);
  
  set(gca,'YTick',(-size(S1,1)+1)*alpha:alpha:0);
  set(gca,'YTicklabel',fliplr(n1));
  
  if (~isempty(titlestr))
    title(titlestr);
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 case 'histogram'
  % icaplot(mode, s1, n1, range, bins, style)
  if length(varargin) < 1, error('Not enough arguments.'); end
  if length(varargin) < 5, style = 'bar';else style = varargin{5}; end
  if length(varargin) < 4, bins = 10;else bins = varargin{4}; end
  if length(varargin) < 3, range = 0;else range = varargin{3}; end
  if length(varargin) < 2, n1 = 0;else n1 = varargin{2}; end
  s1 = varargin{1};
  range = chkrange(range, s1);
  n1 = chkn(n1, s1);
  
  numSignals = size(n1, 2);
  rows = floor(sqrt(numSignals));
  columns = ceil(sqrt(numSignals));
  while (rows * columns < numSignals)
    columns = columns + 1;
  end
  
  switch style
   case {'', 'bar'}
    for i = 1:numSignals,
      subplot(rows, columns, i);
      hist(s1(n1(i), range), bins);
      title(int2str(n1(i)));
      drawnow;
    end
    
   case 'line'
    for i = 1:numSignals,
      subplot(rows, columns, i);
      [Y, X]=hist(s1(n1(i), range), bins);
      plot(X, Y);
      title(int2str(n1(i)));
      drawnow;
    end
   otherwise
    fprintf('Unknown style.\n')
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 case 'scatter'
  % icaplot(mode, s1, n1, s2, n2, range, titlestr, xlabelstr, ylabelstr, markerstr)
  if length(varargin) < 4, error('Not enough arguments.'); end
  if length(varargin) < 9, markerstr = '.';else markerstr = varargin{9}; end
  if length(varargin) < 8, ylabelstr = 'Signal 2';else ylabelstr = varargin{8}; end
  if length(varargin) < 7, xlabelstr = 'Signal 1';else xlabelstr = varargin{7}; end
  if length(varargin) < 6, titlestr = '';else titlestr = varargin{6}; end
  if length(varargin) < 5, range = 0;else range = varargin{5}; end
  n2 = varargin{4};
  s2 = varargin{3};
  n1 = varargin{2};
  s1 = varargin{1};
  range = chkrange(range, s1);
  n1 = chkn(n1, s1);
  n2 = chkn(n2, s2);
  
  rows = size(n1, 2);
  columns = size(n2, 2);
  for r = 1:rows
    for c = 1:columns
      subplot(rows, columns, (r-1)*columns + c);
      plot(s1(n1(r), range),s2(n2(c), range),markerstr);
      if (~isempty(titlestr))
	title(titlestr);
      end
      if (rows*columns == 1)
	xlabel(xlabelstr);
	ylabel(ylabelstr);
      else 
	xlabel([xlabelstr ' (' int2str(n1(r)) ')']);
	ylabel([ylabelstr ' (' int2str(n2(c)) ')']);
      end
      drawnow;
    end
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 case {'compare', 'sum', 'sumerror'}
  % icaplot(mode, s1, n1, s2, n2, range, xrange, titlestr, s1label, s2label)
  if length(varargin) < 4, error('Not enough arguments.'); end
  if length(varargin) < 9, s2label = 'IC';else s2label = varargin{9}; end
  if length(varargin) < 8, s1label = 'Mix';else s1label = varargin{8}; end
  if length(varargin) < 7, titlestr = '';else titlestr = varargin{7}; end
  if length(varargin) < 6, xrange = 0;else xrange = varargin{6}; end
  if length(varargin) < 5, range = 0;else range = varargin{5}; end
  s1 = varargin{1};
  n1 = varargin{2};
  s2 = varargin{3};
  n2 = varargin{4};
  range = chkrange(range, s1);
  xrange = chkxrange(xrange, range);
  n1 = chkn(n1, s1);
  n2 = chkn(n2, s2);

  numSignals = size(n1, 2);
  if (numSignals > 1)
    externalLegend = 1;
  else
    externalLegend = 0;
  end
  
  rows = floor(sqrt(numSignals+externalLegend));
  columns = ceil(sqrt(numSignals+externalLegend));
  while (rows * columns < (numSignals+externalLegend))
    columns = columns + 1;
  end
  
  clf;
  
  for j = 1:numSignals
    subplot(rows, columns, j);
    switch mode
     case 'compare'
      plotcompare(s1, n1(j), s2,n2, range, xrange);
      [legendtext,legendstyle]=legendcompare(n1(j),n2,s1label,s2label,externalLegend);
     case 'sum'
      plotsum(s1, n1(j), s2,n2, range, xrange);
      [legendtext,legendstyle]=legendsum(n1(j),n2,s1label,s2label,externalLegend);
     case 'sumerror'
      plotsumerror(s1, n1(j), s2,n2, range, xrange);
      [legendtext,legendstyle]=legendsumerror(n1(j),n2,s1label,s2label,externalLegend);
    end
    
    if externalLegend
      title([titlestr ' (' s1label  ' ' int2str(n1(j)) ')']);
    else
      legend(char(legendtext));
      if (~isempty(titlestr))
	title(titlestr);
      end
    end
  end
  
  if (externalLegend)
    subplot(rows, columns, numSignals+1);
    legendsize = size(legendtext, 2);
    hold on;
    for i=1:legendsize
      plot([0 1],[legendsize-i legendsize-i], char(legendstyle(i)));
      text(1.5, legendsize-i, char(legendtext(i)));
    end
    hold off;
    axis([0 6 -1 legendsize]);
    axis off;
  end
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotcompare(s1, n1, s2, n2, range, xrange);
  style=getStyles;
  K = regress(s1(n1,:)',s2');
  plot(xrange, s1(n1,range), char(style(1)));
  hold on
  for i=1:size(n2,2)
    plotstyle=char(style(i+1));
    plot(xrange, K(n2(i))*s2(n2(i),range), plotstyle);
  end
  hold off

function [legendText, legendStyle]=legendcompare(n1, n2, s1l, s2l, externalLegend);
  style=getStyles;
  if (externalLegend)
    legendText(1)={[s1l ' (see the titles)']};
  else
    legendText(1)={[s1l ' ', int2str(n1)]};
  end
  legendStyle(1)=style(1);
  for i=1:size(n2, 2)
    legendText(i+1) = {[s2l ' ' int2str(n2(i))]};
    legendStyle(i+1) = style(i+1);
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotsum(s1, n1, s2, n2, range, xrange);
  K = diag(regress(s1(n1,:)',s2'));
  sigsum = sum(K(:,n2)*s2(n2,:));
  plot(xrange, s1(n1, range),'k-', ...
       xrange, sigsum(range), 'b-');

function [legendText, legendStyle]=legendsum(n1, n2, s1l, s2l, externalLegend);
  if (externalLegend)
    legendText(1)={[s1l ' (see the titles)']};
  else
    legendText(1)={[s1l ' ', int2str(n1)]};
  end
  legendText(2)={['Sum of ' s2l ': ', int2str(n2)]};
  legendStyle={'k-';'b-'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotsumerror(s1, n1, s2, n2, range, xrange);
  K = diag(regress(s1(n1,:)',s2'));
  sigsum = sum(K(:,n2)*s2(n2,:));
  plot(xrange, s1(n1, range),'k-', ...
       xrange, sigsum(range), 'b-', ...
       xrange, s1(n1, range)-sigsum(range), 'r-');

function [legendText, legendStyle]=legendsumerror(n1, n2, s1l, s2l, externalLegend);
  if (externalLegend)
    legendText(1)={[s1l ' (see the titles)']};
  else
    legendText(1)={[s1l ' ', int2str(n1)]};
  end
  legendText(2)={['Sum of ' s2l ': ', int2str(n2)]};
  legendText(3)={'"Error"'};
  legendStyle={'k-';'b-';'r-'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function style=getStyles;
  color = {'k','r','g','b','m','c','y'};
  line = {'-',':','-.','--'};
  for i = 0:size(line,2)-1
    for j = 1:size(color, 2)
      style(j + i*size(color, 2)) = strcat(color(j), line(i+1));
    end
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function range=chkrange(r, s)
  if r == 0
    range = 1:size(s, 2);
  else
    range = r;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xrange=chkxrange(xr,r);
  if xr == 0
    xrange = r;
  elseif size(xr, 2) == 2
    xrange = xr(1):(xr(2)-xr(1))/(size(r,2)-1):xr(2);
  elseif size(xr, 2)~=size(r, 2)
    error('Xrange and range have different sizes.');
  else
    xrange = xr;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n=chkn(n,s)
  if n == 0
    n = 1:size(s, 1);
  end