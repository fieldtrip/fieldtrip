function hh = violinplot(x,varargin)
%VIOLINPLOT Plot a vertical violinplot
%   
%  Description
%    VIOLINPLOT(X,Y,OPTIONS) Plot a vertical violin illustrating
%    the density of the data Y at location given by X. X is 1xM
%    vector and Y is NxM matrix. If X is NxM vector and Y is 1xM
%    matrix, the horizontal violin is plotted. 
%  
%    For each column of X and Y, one violinplot and 5%, 50% and 95%
%    lines presenting the estimated distribution of Y is plotted. 
%    The density estimate is made using LGPDENS.
%  
%    H=VIOLINPLOT(X,Y,OPTIONS) returns graphics handles.
%
%    OPTIONS is optional parameter-value pair
%      color  - color used can be given as character or RGB value vector
%               color is used for edgecolor and lighter version for facecolor 
%      range  - tells the estimation range, default is 
%               [min(min(x),mean(x)-3*std(x)), max(max(x),mean(x)+3*std(x))]
%      dots   - indicates whether extreme points outside violin are plotted
%               on or off (default)
%      cutoff - cutoff for not plotting thin tails. Default is
%               0.01, that is, 1% of the mass in both tails is not plotted
%
  
% Copyright (c) 2011 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.
  
  ip=inputParser;
  ip.FunctionName = 'VIOLINPLOT';
  ip.addRequired('x', @(x) isnumeric(x));
  ip.addOptional('y', [], @(x) isnumeric(x));
  ip.addParamValue('color','k', @(x) ischar(x) || ...
                   (isnumeric(x) & isequal(size(x),[1 3])));
  ip.addParamValue('range',[], @(x) isreal(x)&&(length(x)==2||length(x)==4));
  ip.addParamValue('cutoff',.01, @(x) isreal(x)&& isscalar(x) && x>0 && x<1);
  ip.addParamValue('dots', 'off', @(x) islogical(x) || ...
                   ismember(x,{'on' 'off'}))
  ip.parse(x,varargin{:});
  x=ip.Results.x;
  y=ip.Results.y;
  color=ip.Results.color;
  xrange=ip.Results.range;
  cutoff=ip.Results.cutoff;
  dots=ip.Results.dots;
  
  horiz=true;
  if isempty(y)
    y=x;
    x=1:size(y,2);
  elseif (isvector(y) && ~isvector(x)) || isscalar(y)
    horiz=false;
    tmp=x;
    x=y;
    y=tmp;
  end
  
  [n,m]=size(y);
  for i1=1:m
    [p,notused,yy]=lgpdens(y(:,i1),'gridn',200,'range',xrange);
    cp=cumsum(p./sum(p));
    qicutofflo=binsgeq(cp,cutoff/2);
    qicutoffhi=binsgeq(cp,1-cutoff/2);
    p([1:qicutofflo qicutoffhi:end])=[];
    yy([1:qicutofflo qicutoffhi:end])=[];
    cp=cumsum(p./sum(p));
    qi5=binsgeq(cp,0.05);
    qi50=binsgeq(cp,0.5);
    qi95=binsgeq(cp,0.95);
    p=p./max(p)/5;
    % if color was character, this will be rgb vector
    if horiz
      hp(1,i1)=patch(x(i1)+[p; -p(end:-1:1)],[yy; yy(end:-1:1)],color);
    else
      hp(1,i1)=patch([yy; yy(end:-1:1)],x(i1)+[p; -p(end:-1:1)],color);
    end
    if i1==1
      color=get(hp(1,i1),'facecolor');
    end
    set(hp(1,i1),'edgecolor',color);
    set(hp(1,i1),'facecolor',color.*0.2+[.8 .8 .8])
    if horiz
      % median line
      h(1,i1)=line([x(i1)-p(qi50) x(i1)+p(qi50)],[yy(qi50) yy(qi50)],'color',color,'linewidth',1);
      % 5% line
      h(1,i1)=line([x(i1)-p(qi5) x(i1)+p(qi5)],[yy(qi5) yy(qi5)],'color',color,'linewidth',1);
      % 95% line
      h(1,i1)=line([x(i1)-p(qi95) x(i1)+p(qi95)],[yy(qi95) yy(qi95)],'color',color,'linewidth',1);
    else
      % median line
      h(1,i1)=line([yy(qi50) yy(qi50)],[x(i1)-p(qi50) x(i1)+p(qi50)],'color',color,'linewidth',1);
      % 5% line
      h(1,i1)=line([yy(qi5) yy(qi5)],[x(i1)-p(qi5) x(i1)+p(qi5)],'color',color,'linewidth',1);
      % 95% line
      h(1,i1)=line([yy(qi95) yy(qi95)],[x(i1)-p(qi95) x(i1)+p(qi95)],'color',color,'linewidth',1);
    end
    if isequal(dots,'on')
      ydots=y(y(:,i1)>max(yy),i1);
      if ~isempty(ydots)
        if horiz
          line(x(i1),ydots,'marker','.','linestyle','none','color',color)
        else
          line(ydots,x(i1),'marker','.','linestyle','none','color',color)
        end
      end
      ydots=y(y(:,i1)<min(yy),i1);
      if ~isempty(ydots)
        if horiz
          line(x(i1),ydots,'marker','.','linestyle','none','color',color)
        else
          line(ydots,x(i1),'marker','.','linestyle','none','color',color)
        end
      end
    end
  end
  if nargout>0
    hh=[hp; h];
    hh=hh(:);
  end
  