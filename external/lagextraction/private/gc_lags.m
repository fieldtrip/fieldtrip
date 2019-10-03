function [lags,use_ascend,dc,sc,E] = gc_lags(points,options)

% GC_LAGS       Run lag extraction using a binary graph cut
%
%   Run lag extraction on a reordered raster plot (2D Image) using a binary graph cut
%
%   SYNTAX
%       [LAGS,USE_ASCEND,DC,SC] = GC_LAGS(POINTS,OPTIONS)
%
%   DC is the data cost
%   SC is the smoothness cost
%

% $Id: gc_lags.m 4 2009-08-15 21:10:35Z gramfort $
% $LastChangedBy: gramfort $
% $LastChangedDate: 2009-08-15 17:10:35 -0400 (Sam, 15 aoû 2009) $
% $Revision: 4 $

if nargin<2
    options.null = 0;
end

if ~isfield(options, 'alpha')
    options.alpha = 0.01;
end
alpha = options.alpha;

if ~isfield(options, 'verbose')
    options.verbose = true;
end
verbose = options.verbose;

if ~isfield(options, 'use_lcurve')
    options.use_lcurve = false;
end
use_lcurve = options.use_lcurve;

if ~isfield(options, 'maxit') % maximum number of iterations for alpha search in l-curve
    options.maxit = 40;
end
maxit = options.maxit;

if ~isfield(options, 'dc_tol') % tolerance on data cost for alpha search in l-curve
    options.dc_tol = 1e-25;
end
dc_tol = options.dc_tol;

if isfield(options, 'null')
    options = rmfield(options,'null');
end

%%%%%%%% END OPTIONS %%%%%%%%%

points = - points; % min cut through maximum of points values
points = points - min(points(:)); % only positive capacities for the maxflow

npoints = size(points,1);

if options.disp_log
    disp(['Nb trials : ',num2str(npoints)]);
    disp(['Nb time samples : ',num2str(size(points,2))]);
    disp(['Lambda : ',num2str(alpha)]);
end

if use_lcurve
    disp('Computing best alpha with L-curve');
    verbose = false;
    sc = Inf;
    ht = java.util.Hashtable;

    % Add alpha = 0 to ht
    [lags0,use_ascend,dc0,sc0] = gc_lags_aux(0.0);

    % Find the biggest alpha
    while sc > 0
        fprintf('.');
        [lags,use_ascend,dc,sc] = gc_lags_aux(alpha);
        ht.put(alpha,[dc sc]);
        alpha = 3*alpha;
        % alpha = 2.5*alpha;
        % alpha = 5*alpha;
        % alpha = 2*alpha;
        % alpha = 1.5*alpha;
    end
    fprintf('\n Biggest alpha : %f\n',alpha);
    % Try to find best alpha iteratively
    goon = true;
    niter = 1;
    while goon
        [alphas,dcs,scs] = ht_to_array();

        % Make sure dcs and scs values are unique
        [tmp,dcs_idx] = unique(dcs,'first');
        [tmp,scs_idx] = unique(scs,'first');
        good_idx = sort(intersect(dcs_idx,scs_idx));
        dcs = dcs(good_idx);
        scs = scs(good_idx);
        alphas = alphas(good_idx);
        clear good_idx

        scs_intermediate = diff(dcs+alphas.*scs) ./ diff(alphas);
        dcs_intermediate = (dcs(1:end-1)+alphas(1:end-1).*scs(1:end-1)) - alphas(1:end-1) .* scs_intermediate;
        dcs_interp = linspace(min(dcs),max(dcs),500);
        scs_interp = interp1([dcs;dcs_intermediate],[scs;scs_intermediate],dcs_interp,'linear');

        l1 = 1./alphas(1:end-1);
        l2 = 1./alphas(2:end);
        angle_intermediate = acos((1+l1.*l2)./(sqrt(1+l1.^2).*sqrt(1+l2.^2))); % angle between neighboring segments

        angle_intermediate(1) = 0;
        angle_intermediate(end) = 0;
        midx = localmax(angle_intermediate);

        [tmp,best_midx] = max(angle_intermediate);

        best_dcs = dcs_intermediate(best_midx);
        best_scs = scs_intermediate(best_midx);

        best_dcs_all = best_dcs;
        best_scs_all = best_scs;

        % get alpha by computing harmonic mean of neighboring alphas (better than arithmetic mean)
        new_alphas = alphas(midx) .* alphas(midx+1) ./ (alphas(midx)+alphas(midx+1)) * 2;

        % Sort alphas using order of angle differences
        [tmp,I] = sort(angle_intermediate(midx));
        new_alphas = new_alphas(I);

        precision = 0;

        needs_update = false;
        for ii = 1:length(new_alphas)
            alpha = new_alphas(ii); % take in between alpha

            if isempty(ht.get(alpha)) % Check if solution with this alpha has already between computed
                disp(['    Computing solution for alpha : ',num2str(alpha)]);
                [lags,use_ascend,dc,sc] = gc_lags_aux(alpha);
                ht.put(alpha,[dc sc]);
                needs_update = true; % Say if at least one new alpha is added
            else
                sc_dc = ht.get(alpha);
                dc = sc_dc(1);
                sc = sc_dc(2);
            end
            precision = max(precision,min(abs(dcs-dc)));
        end

        if ~needs_update
            disp('Stopping search (all alphas already computed)');
            goon = false;
        end

        disp([' Error : ',num2str(precision)]);
        if precision < dc_tol | niter >= maxit
            goon = false;
        end
        niter = niter+1;

        if 1 % set best final alpha with interpolation
            if ~goon
                [alpha,best_dcs,best_scs] = get_best_alpha();
            end
        end

        if 1 % show l-curve
            smart_figure('Lag extraction L-curve'); clf
            title('Lag extraction L-curve')
            if goon
                color = {'g','b'};
                linewidth = 1;
            else
                color = {'r','r'};
                linewidth = 3;
            end
            USE_LOG = false;
            % USE_LOG = true;
            if USE_LOG % plot log-log or not
                loglog(dcs_interp(dcs_interp>0),scs_interp(dcs_interp>0),color{2});
                hold on
                loglog(dcs(dcs>0),scs(dcs>0),[color{2},'x']);
            else
                plot(dcs_interp,scs_interp,color{2},'LineWidth',linewidth);
                hold on
                plot(dcs,scs,[color{2},'x']);
            end

            for k=1:length(dcs)
                if alphas(k) > 0
                    flow = dcs(k)+alphas(k)*scs(k);
                    line([0 , flow],[flow/alphas(k), 0],'Color','k');
                end
            end

            axis tight
            line([min(dcs_interp)/100,best_dcs],[best_scs,best_scs], ...
                 'LineStyle','--','Color',color{1},'LineWidth',linewidth);
            line([best_dcs,best_dcs],[min(scs_interp)/100,best_scs], ...
                 'LineStyle','--','Color',color{1},'LineWidth',linewidth);
            hold off;
            xlabel('Data term')
            ylabel('Smoothness term')

            if 1 % show max curvature points
                hold on
                if USE_LOG
                    loglog(best_dcs_all,best_scs_all,'bo');
                else
                    plot(best_dcs_all,best_scs_all,'bo');
                end
                hold off
            end
        end
    end
    disp(['L-curve done with final alpha : ',num2str(alpha)]);
    disp(['         tol   : ',num2str(precision)]);
    disp(['         niter : ',num2str(niter)]);
else
    [lags,use_ascend,dc,sc] = gc_lags_aux(alpha);
end

E = dc+alpha*sc;

if options.disp_log
    disp('E = DC + alpha * SC');
    disp(sprintf('%f = %f + %f * %f',dc+alpha*sc,dc,alpha,sc));
end

% ****************************************************************************** %
function [best_alpha,best_dcs,best_scs] = get_best_alpha()

[local_alphas,local_dcs,local_scs] = ht_to_array();
offset = 1; % don't consider extremities of the curve
local_dcs = local_dcs(1+offset:end-offset);
local_scs = local_scs(1+offset:end-offset);
local_alphas = local_alphas(1+offset:end-offset);

local_dcs_interp = linspace(min(local_dcs),max(local_dcs),100);
local_scs_interp = interp_cubic_herm(local_dcs,local_scs,-1./local_alphas,local_dcs_interp);
local_curvature = curvature_2D_diff(local_dcs_interp,local_scs_interp);
local_curvature = abs(local_curvature); % use absolute value
local_curvature(1:5) = NaN;
local_curvature(end-5:end) = NaN;

% [tmp,idx] = max(local_curvature);
[tmp,idx] = min(local_curvature);
best_dcs = local_dcs_interp(idx);
best_scs = local_scs_interp(idx);

[tmp,local_alpha_idx] = min(abs(local_dcs-best_dcs));

best_alpha = local_alphas(local_alpha_idx);

if 0
    smart_figure('Lag extraction L-curve'); hold on
    plot(local_dcs_interp,local_scs_interp,'g')
end

if 0
    smart_figure('SC vs alpha');
    plot(local_alphas,local_scs)
    xlabel('alpha')
    ylabel('SC')
    hold on
    plot(best_alpha,best_scs,'ro');
    hold off
end

if verbose
    smart_figure('Curve interpolation');
    plot(local_dcs_interp,local_scs_interp);
    hold on
    plot(best_dcs,best_scs,'ro');
    hold off
    title('Curve interpolation');
end

end % end get_best_alpha function

% ****************************************************************************** %
function [alphas,dcs,scs] = ht_to_array()
    alphas = zeros(ht.size(),1);
    dcs = zeros(ht.size(),1);
    scs = zeros(ht.size(),1);
    keys = ht.keys();
    cnt = 1;
    while keys.hasNext()
       alpha = keys.next();
       dc_sc = ht.get(alpha);
       dc = dc_sc(1);
       sc = dc_sc(2);
       alphas(cnt) = alpha;
       dcs(cnt) = dc;
       scs(cnt) = sc;
       cnt = cnt+1;
    end;
    [tmp,I] = sort(alphas);
    alphas = alphas(I);
    dcs = dcs(I);
    scs = scs(I);

    dcs = dcs - dc0; % remove offset
end

% ****************************************************************************** % 
function idx = localmax(x)
% Find local maxima in an array

idx = find( diff( sign( diff([0; x(:); 0]) ) ) < 0 );

end % end function

% ****************************************************************************** %

function [lags,use_ascend,dc,sc] = gc_lags_aux(alpha)

if options.disp_log
    disp('Testing ascend')
end

order = [npoints:-1:1];
[lags_ascend,flow_ascend,labels_ascend] = gc_aux_mex(points(order,:),alpha);
lags_ascend = lags_ascend(order);
labels_ascend = labels_ascend(order,:);

% Compute manually resulting data cost and smoothness cost
%
% lag_ind = sub2ind(size(points(order,:)),[1:size(points,1)]',lags_ascend(:));
% dc = sum(points(lag_ind))
% sc = sum(abs(diff(lags_ascend)))
% flow_ascend
% fl = dc + alpha * sc

if options.disp_log
    disp('Testing descend')
end
[lags_descend,flow_descend,labels_descend] = gc_aux_mex(points,alpha);

if flow_ascend < flow_descend
    lags = lags_ascend;
    use_ascend = true;
    flow = flow_ascend;
else
    lags = lags_descend;
    use_ascend = false;
    flow = flow_descend;
end

sc = sum(abs(diff(lags))); % Smoothness cost
dc = flow - alpha * sc; % Data cost

end %  function

end
