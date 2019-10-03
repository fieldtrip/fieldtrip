function [data_aligned,times_aligned] = perform_realign(data_reordered,times,lags)

% PERFORM_REALIGN Time series based on extracted lags
%
%   SYNTAX
%       [DATA_ALIGNED,TIMES_ALIGNED] = PERFORM_REALIGN(DATA_REORDERED,TIMES,LAGS)
%

% $Id: perform_realign.m 2 2009-06-16 19:24:10Z gramfort $
% $LastChangedBy: gramfort $
% $LastChangedDate: 2009-06-16 15:24:10 -0400 (Mar, 16 jui 2009) $
% $Revision: 2 $

mean_lag = fix(mean(lags));
data_aligned = inf*ones(size(data_reordered));

delta_max = 1;
for k=1:length(lags)
    lag = lags(k);
    delta = lag-mean_lag;
    delta_max = max(delta_max,abs(delta));
    if delta > 0
        point = data_reordered(k,delta:end);
        data_aligned(k,1:length(point)) = point;
    else
        point = data_reordered(k,1:(end+delta));
        data_aligned(k,(end-length(point)+1):end) = point;
    end
end

times_aligned = times - times(delta_max) + times(1);

