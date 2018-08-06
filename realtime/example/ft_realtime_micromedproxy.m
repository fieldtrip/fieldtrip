% FT_REALTIME_MICROMEDPROXY reads continuous data from a Micromed acquisition system
% through the BCI interface and writes it to a FieldTrip buffer.
%
% The FieldTrip buffer is a network transparent server that allows the acquisition
% client to stream data to it. An analysis client can connect to read the data upon
% request. Multiple clients can connect simultaneously, each analyzing a specific
% aspect of the data concurrently.
%
% Use as
%   ft_realtime_micromedproxy(cfg)
%
% The configuration should contain
%   cfg.feedback             = 'yes' or 'no' (default = 'no')
%
% The target to write the data to is configured as
%   cfg.target.datafile      = string, target destination for the data (default = 'buffer://localhost:1972')
%   cfg.target.dataformat    = string, default is determined automatic
%   cfg.target.eventfile     = string, target destination for the events (default = 'buffer://localhost:1972')
%   cfg.target.eventformat   = string, default is determined automatic
%
% To stop this realtime function, you have to press Ctrl-C
%
% This function is released together with the FieldTrip toolbox as
% precompiled *.p file. The original code to this function is not released
% as open source under the GLP, but is to Micromed users available upon
% request. You can contact Robert Oostenveld (FCDC) or Cristiano Rizzo
% (Micromed) to get a copy of the original source code.
%
% Copyright (C) 2009, Erik Aarnoutse, UMC Utrecht, The Netherlands
% Copyright (C) 2009, Cristiano Rizzo, Micromed A.S., Italy
% Copyright (C) 2009, Robert Oostenveld, Donders Institute Nijmegen, The Netherlands
%
% See also FT_REALTIME_SIGNALPROXY, FT_REALTIME_SIGNALVIEWER
