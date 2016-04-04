function fun = nmt_ts_intervalpower(fun,meantype);
% when time interval is selected, calculate desired representation of interval power
switch(meantype)
    case 'msq' % mean power
        fun = mean(fun.^2,2);
    case 'rms' % RMS power
        fun = sqrt(mean(fun.^2,2));
    case 'mean'  % simple average
        fun = mean(fun,2);
end