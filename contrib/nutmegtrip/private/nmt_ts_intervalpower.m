function fun = nmt_ts_intervalpower(fun);
% when time interval is selected, calculate desired representation of interval power
switch(2)
    case 1 % mean power
        fun = mean(fun.^2,2);
    case 2 % RMS power
        fun = sqrt(mean(fun.^2,2));
    case 3  % simple average
        fun = mean(fun,2);
end