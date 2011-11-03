% FT_PREAMBLE_CALLINFO is a helper script that records the time and memory
% at the start of the function. This is to be used together with
% FT_POSTAMBLE_CALLINFO which will record the time and memory at the end of
% the function.

% record the start time and memory
ftFuncTimer = tic();
ftFuncClock = clock();
ftFuncMem   = memtic();
