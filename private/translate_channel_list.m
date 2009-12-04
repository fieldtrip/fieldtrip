function [select] = translate_channel_list(varargin)

% TRANSLATE_CHANNEL_LIST is an obsolete function
% This function is replaced by CHANNELSELECTION

warning('please use the new framework function channelselection.m');
select = channelselection(varargin{:});
