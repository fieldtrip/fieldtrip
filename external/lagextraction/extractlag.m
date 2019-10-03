% extractlag() -
%
% Usage:
%   >>  [order,lags,coords,E_lags] = extractlag( points, options );
%
% Inputs:
%   points     - first input of the function
%   options    - Various options
%
% Outputs:
%   order     - permutation to apply to reorder points
%   lags      - value of computed lags
%   coords    - points coordinates in low dimensional space
%   E_lags    - energy of lag extraction procedure
%
% See also:
%   POP_EXTRACTLAG, EEGLAB

% (c) Copyright 2008-2009 Alexandre Gramfort. All Rights Reserved.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Id: extractlag.m 4 2009-08-15 21:10:35Z gramfort $
% $LastChangedBy: gramfort $
% $LastChangedDate: 2009-08-15 17:10:35 -0400 (Sam, 15 aoû 2009) $
% $Revision: 4 $

function [order,lags,coords,E_lags] = extractlag( points, options )

if nargin < 1
    help extractlag;
    return;
end;

if nargin<2
    options.null = 0;
end

sigmas = options.sigma;
best_E_lags = Inf;

disp(['- Running lag extraction']);
disp(['--- Nb trials : ',num2str(size(points,1))]);
disp(['--- Nb time samples : ',num2str(size(points,2))]);

E = [];
for ii=1:length(sigmas)
    options.sigma = sigmas(ii);
    %% Embedding
    [order,coords,points] = perform_embedding(points,options);

    if isfield(options, 'order') % Hack to force order
        order = options.order;
    end

    %% Extracting lags
    [lags,use_ascend,E_lags] = perform_extraction(points(order,:),options);

    if ii==1 || E_lags < best_E_lags
        best_lags       = lags;
        best_order      = order;
        best_use_ascend = use_ascend;
        best_E_lags     = E_lags;
        best_sigma      = options.sigma;
        best_ii         = ii;
    end
    E(ii) = E_lags;
end

lags       = best_lags;
order      = best_order;
use_ascend = best_use_ascend;
E_lags     = best_E_lags;

if ~use_ascend
    order = order(end:-1:1);
    lags = lags(end:-1:1);
end

if length(sigmas)>1 & options.verbose
    smart_figure('E lags vs sigmas');
    plot(sigmas,E);
    yl = ylim;
    hold on
    line([best_sigma best_sigma],[yl(1) E(best_ii)],'Color','k');
    hold off
    xlabel('Sigma');
    ylabel('E lags');
end

% keyboard

disp(['---------- Using sigma = ',num2str(best_sigma)]);

end