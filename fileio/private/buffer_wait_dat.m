function available = buffer_wait_dat(selection, host, port)

% BUFFER_WAIT_DAT implementation that is also backwards compatibility with ft buffer version 1
%
% Use as
%   available = buffer_wait_dat(selection, host, port)
% where
%   selection(1) = nsamples, 0 indicates not to wait
%   selection(2) = nevents,  0 indicates not to wait
%   selection(3) = timeout in seconds
%
% It returns a structure with the available nsamples and nevents.

selection(1) = selection(1)-1;
selection(2) = selection(2)-1;
selection(3) = selection(3)*1000; % in miliseconds

% check if backwards compatibilty mode is required
try
  % the WAIT_DAT request waits until it has more samples or events
  % the following should work for buffer version 2
  available = buffer('WAIT_DAT', selection, host, port);
  
catch
  % the error means that the buffer is version 1, which does not support the WAIT_DAT request
  % the wait_dat can be implemented by polling the buffer
  
  timeout   = selection(3)*1000; % in seconds
  stopwatch = tic;
  
  % results are retrieved in the order written to the buffer
  orig = buffer('GET_HDR', [], host, port);
  
  if timeout > 0
    % wait maximal timeout seconds until more than nsamples samples or nevents events have been received
    while toc(stopwatch)<timeout
      if nsamples == -1 && nevents == -1,             break, end
      if nsamples ~= -1 && orig.nsamples >= nsamples, break, end
      if nevents  ~= -1 && orig.nevents  >= nevents,  break, end
      orig = buffer('GET_HDR', [], host, port);
      pause(0.001);
    end
  else
    % no reason to wait
  end
  
  available.nsamples = orig.nsamples;
  available.nevents  = orig.nevents;
end % try buffer v1 or v2
