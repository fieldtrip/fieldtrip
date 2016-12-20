function fun = nmt_ts_intervalpower(fun,meantype)
% when time interval is selected, calculate desired representation of interval power
switch(meantype)
    case 'msq' % mean power
        fun = mean(fun.^2,3); % average over freq dimension first, if present
        fun = mean(fun.^2,2);
    case 'rms' % RMS power
        fun = sqrt(mean(fun.^2,3));
        fun = sqrt(mean(fun.^2,2));
    case 'mean'  % simple average
        fun = mean(fun,3);
        fun = mean(fun,2);
end
