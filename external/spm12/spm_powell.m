function [p,f] = spm_powell(p,xi,tolsc,func,varargin)
% Powell optimisation method
% FORMAT [p,f] = spm_powell(p,xi,tolsc,func,varargin)
%   p        - Starting parameter values
%   xi       - columns containing directions in which to begin searching
%   tolsc    - stopping criteria, optimisation stops when
%                sqrt(sum(((p-p_prev)./tolsc).^2))<1
%   func     - name of evaluated function
%   varargin - remaining arguments to func (after p)
%
%   p        - final parameter estimates
%   f        - function value at minimum
%__________________________________________________________________________
%
% Method is based on Powell's optimisation method described in
% Numerical Recipes (Press, Flannery, Teukolsky & Vetterling).
%__________________________________________________________________________
% Copyright (C) 2001-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_powell.m 4423 2011-08-04 16:28:51Z guillaume $


p = p(:);
f = feval(func,p,varargin{:});
for iter=1:512
    %if numel(p)>1, fprintf('iteration %d...\n', iter); end;            %-#
    ibig = numel(p); 
    pp   = p;
    fp   = f;
    del  = 0;
    for i=1:length(p)
        ft = f;
        [p,junk,f] = min1d(p,xi(:,i),func,f,tolsc,varargin{:});
        if abs(ft-f) > del,
            del  = abs(ft-f);
            ibig = i;
        end
    end
    if numel(p)==1 || sqrt(sum(((p(:)-pp(:))./tolsc(:)).^2))<1, return; end
    ft = feval(func,2.0*p-pp,varargin{:});
    if ft < f
        [p,xi(:,ibig),f] = min1d(p,p-pp,func,f,tolsc,varargin{:});
    end
end
warning('Too many optimisation iterations');


%==========================================================================
% function [p,pi,f] = min1d(p,pi,func,f,tolsc,varargin)
%==========================================================================
function [p,pi,f] = min1d(p,pi,func,f,tolsc,varargin)
% Line search for minimum.

global lnm % used in funeval
lnm      = struct('p',p,'pi',pi,'func',func,'args',[]);
lnm.args = varargin;

min1d_plot('Init', 'Line Minimisation','Function','Parameter Value');
min1d_plot('Set', 0, f);

tol      = 1/sqrt(sum((pi(:)./tolsc(:)).^2));
t        = bracket(f);
[f,pmin] = search(t,tol);
pi       = pi*pmin;
p        = p + pi;

%if length(p)<12,
%    for i=1:length(p), fprintf('%-8.4g ', p(i)); end;                  %-#
%    fprintf('| %.5g\n', f);                                            %-#
%else
%    fprintf('%.5g\n', f);                                              %-#
%end
min1d_plot('Clear');


%==========================================================================
% function f = funeval(p)
%==========================================================================
function f = funeval(p)
% Reconstruct parameters and evaluate.

global lnm % defined in min1d
pt = lnm.p+p.*lnm.pi;
f  = feval(lnm.func,pt,lnm.args{:});
min1d_plot('Set',p,f);


%==========================================================================
% function t = bracket(f)
%==========================================================================
function t = bracket(f)
% Bracket the minimum (t(2)) between t(1) and t(3)

gold   = (1+sqrt(5))/2; % Golden ratio

t(1)   = struct('p',0,'f',f);
t(2).p = 1;
t(2).f = funeval(t(2).p);

% if t(2) not better than t(1) then swap
if t(2).f > t(1).f
    t(3) = t(1);
    t(1) = t(2);
    t(2) = t(3);
end

t(3).p = t(2).p + gold*(t(2).p-t(1).p);
t(3).f = funeval(t(3).p);

while t(2).f > t(3).f

    % fit a polynomial to t
    tmp = cat(1,t.p)-t(2).p;
    pol = pinv([ones(3,1) tmp tmp.^2])*cat(1,t.f);

    % minimum is when gradient of polynomial is zero
    % sign of pol(3) (the 2nd deriv) should be +ve
    if pol(3)>0
        % minimum is when gradient of polynomial is zero
        d    = -pol(2)/(2*pol(3)+eps);

        % A very conservative constraint on the displacement
        if d > (1+gold)*(t(3).p-t(2).p),
            d = (1+gold)*(t(3).p-t(2).p);
        end
        u.p  = t(2).p+d;
    else
        % sign of pol(3) (the 2nd deriv) is not +ve
        % so extend out by golden ratio instead
        u.p  = t(3).p+gold*(t(3).p-t(2).p);
    end

    % FUNCTION EVALUATION
    u.f  = funeval(u.p);

    if (t(2).p < u.p) == (u.p < t(3).p)

        % u is between t(2) and t(3)
        if u.f < t(3).f
            % minimum between t(2) and t(3) - done
            t(1) = t(2);
            t(2) = u;
            return
        elseif u.f > t(2).f
            % minimum between t(1) and u - done
            t(3) = u;
            return;
        end
    end

    % Move all 3 points along
    t(1) = t(2);
    t(2) = t(3);
    t(3) = u;
end


%==========================================================================
% function [f,p] = search(t, tol)
%==========================================================================
function [f,p] = search(t, tol)
% Brent's method for line searching - given that minimum is bracketed

gold1 = 1-(sqrt(5)-1)/2;

% Current and previous displacements
d     = Inf;
pd    = Inf;

% sort t into best first order
[junk,ind] = sort(cat(1,t.f));
t   = t(ind);
brk = [min(cat(1,t.p)) max(cat(1,t.p))];

for iter=1:128
    % check stopping criterion
    if abs(t(1).p - 0.5*(brk(1)+brk(2)))+0.5*(brk(2)-brk(1)) <= 2*tol
        p = t(1).p;
        f = t(1).f;
        return;
    end

    % keep last two displacents
    ppd = pd;
    pd  = d;

    % fit a polynomial to t
    tmp = cat(1,t.p)-t(1).p;
    pol = pinv([ones(3,1) tmp tmp.^2])*cat(1,t.f);

    % minimum is when gradient of polynomial is zero
    d   = -pol(2)/(2*pol(3)+eps);
    u.p = t(1).p+d;

    % check so that displacement is less than the last but two,
    % that the displaced point is between the brackets
    % and that the solution is a minimum rather than a maximum
    eps2 = 2*eps*abs(t(1).p)+eps;
    if abs(d) > abs(ppd)/2 || u.p < brk(1)+eps2 || u.p > brk(2)-eps2 || pol(3)<=0
        % if criteria are not met, then golden search into the larger part
        if t(1).p >= 0.5*(brk(1)+brk(2)),
            d = gold1*(brk(1)-t(1).p);
        else
            d = gold1*(brk(2)-t(1).p);
        end
        u.p = t(1).p+d;
    end

    % FUNCTION EVALUATION
    u.f = funeval(u.p);

    % Insert the new point into the appropriate position and update
    % the brackets if necessary
    if u.f <= t(1).f
        if u.p >= t(1).p, brk(1)=t(1).p; else brk(2)=t(1).p; end
        t(3) = t(2);
        t(2) = t(1);
        t(1) = u;
    else
        if u.p < t(1).p, brk(1)=u.p; else brk(2)=u.p; end
        if u.f <= t(2).f
            t(3) = t(2);
            t(2) = u;
        elseif u.f <= t(3).f
            t(3) = u;
        end
    end
end


%==========================================================================
% function min1d_plot(action,arg1,arg2,arg3)
%==========================================================================
function min1d_plot(action,arg1,arg2,arg3)
% Visual output for line minimisation
persistent min1dplot

if ~nargin, action = 'Init'; end

% Find the Interactive window and exit if not
%--------------------------------------------------------------------------
fg = spm_figure('FindWin','Interactive');
if isempty(fg), return; end

%-Initialize
%--------------------------------------------------------------------------
if strcmpi(action,'init')
    if nargin<4, arg3 = 'Function';          end
    if nargin<3, arg2 = 'Value';             end
    if nargin<2, arg1 = 'Line minimisation'; end
    
    min1dplot = struct('pointer',get(fg,'Pointer'),...
                       'name',   get(fg,'Name'),...
                       'ax',     [],...
                       'buffer', get(fg,'DoubleBuffer'));
    min1d_plot('Clear');
    set(fg,'Pointer','Watch');
    set(fg,'DoubleBuffer','on');
    min1dplot.ax = axes('Position', [0.15 0.1 0.8 0.75],...
                        'Box',      'on',...
                        'Parent',   fg);
    lab = get(min1dplot.ax,'Xlabel');
    set(lab,'string',arg3,'FontSize',10);
    lab = get(min1dplot.ax,'Ylabel');
    set(lab,'string',arg2,'FontSize',10);
    lab = get(min1dplot.ax,'Title');
    set(lab,'string',arg1);
    line('Xdata',[], 'Ydata',[],...
        'LineWidth',2,'Tag','LinMinPlot',...
        'LineStyle','-','Marker','o',...
        'Parent',min1dplot.ax);
    drawnow;
    
%-Reset
%--------------------------------------------------------------------------
elseif strcmpi(action,'set')
    br = findobj(fg,'Tag','LinMinPlot');
    if ~isempty(br)
        [xd,indx] = sort([get(br,'Xdata') arg1]);
        yd = [get(br,'Ydata') arg2];
        yd = yd(indx);
        set(br,'Ydata',yd,'Xdata',xd);
        drawnow;
    end
    
%-Clear
%--------------------------------------------------------------------------
elseif strcmpi(action,'clear')
    fg = spm_figure('FindWin','Interactive');
    if isstruct(min1dplot)
        if ishandle(min1dplot.ax), delete(min1dplot.ax); end
        set(fg,'Pointer',min1dplot.pointer);
        set(fg,'Name',min1dplot.name);
        set(fg,'DoubleBuffer',min1dplot.buffer);
    end
    spm_figure('Clear',fg);
    drawnow;
end
