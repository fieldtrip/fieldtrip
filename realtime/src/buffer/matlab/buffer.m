function [varargout] = buffer(varargin)

% BUFFER manages and accesses the realtime data acquisition buffer
% This function is implented as a mex file.
%
% Use as
%   retval = buffer(cmd, detail, host, port)
%
% To read data from a buffer server over the network
%   hdr = buffer('get_hdr', [],     host, port)
%   dat = buffer('get_dat', datsel, host, port)
%   evt = buffer('get_evt', evtsel, host, port)
%
% The selection for data and events should be zero-offset and contain
%   datsel = [begsample endsample]
%   evtsel = [begevent  endevent ]
%
% To write data to a buffer server over the network
%   buffer('put_hdr', hdr, host, port)
%   buffer('put_dat', dat, host, port)
%   buffer('put_evt', evt, host, port)
%
% To implement a local buffer server and have other clients
% connect to it at a specified network port
%   buffer('tcpserver', 'init', [], port)
%   buffer('tcpserver', 'exit', [], port)

% Copyright (C) 2008-2010, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

error('Could not locate mex file');
