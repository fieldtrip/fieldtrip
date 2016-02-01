function [info] = mne_ex_rt(mne_rt_server_ip, mne_rt_server_commandPort, mne_rt_server_fiffDataPort, p_nBuffers)
%
%   An example of a mne_rt_server real-time connection
%
%   function mne_ex_rt(mne_rt_server_ip, mne_rt_server_commandPort ,mne_rt_server_fiffDataPort)
%
%   mne_rt_server_ip           - IP of the running mne_rt_server
%   mne_rt_server_commandPort  - Command port of the mne_rt_server
%   mne_rt_server_fiffDataPort - Fiff data port of the mne_rt_server
%
%   Returns the measurement info
%   
%

%
%   Author : Christoph Dinh, Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%

if nargin == 3
    p_nBuffers = 100;
elseif nargin ~= 4
    error(me,'Incorrect number of arguments');
end

%% add dynamic java path
javaaddpath(fileparts(mfilename('fullpath')));

%% create command client
t_cmdClient = mne_rt_cmd_client(mne_rt_server_ip, mne_rt_server_commandPort);

%% create data client
t_dataClient = mne_rt_data_client(mne_rt_server_ip, mne_rt_server_fiffDataPort);

%% set data client alias -> for convinience (optional)
t_dataClient.setClientAlias('mne_ex_matlab'); % used in option 2 later on

%% example commands
t_helpInfo = t_cmdClient.sendCommand('help');
fprintf('### Help ###\n%s',t_helpInfo);
t_clistInfo = t_cmdClient.sendCommand('clist');
fprintf('### Client List ###\n%s',t_clistInfo);
t_conInfo = t_cmdClient.sendCommand('conlist');
fprintf('### Connector List ###\n%s',t_conInfo);

%% read meas info
% Option 1
t_aliasOrId = t_dataClient.getClientId();
% Option 2
%t_aliasOrId = 'mne_ex_matlab';
t_cmdClient.requestMeasInfo(t_aliasOrId);
t_measInfo = t_dataClient.readInfo();

%% start measurement
t_cmdClient.requestMeas(t_aliasOrId);

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

figure;
hold on;
h_old=plot(0,0);

t_bIsRunning = true;
t_iCount = 0;
while (t_bIsRunning)
    fprintf('read buffer...');
    [kind, t_matRawBuffer] = t_dataClient.readRawBuffer(t_measInfo.nchan);
    %
    %   You can add your own miracle here
    %
    if(kind == FIFF.FIFF_DATA_BUFFER)
        fprintf('(%d channels x %d samples) [done]\r\n', size(t_matRawBuffer,1), size(t_matRawBuffer,2));
        h = plot(t_matRawBuffer');
        delete(h_old);
        h_old = h;
        drawnow;
    elseif (kind == FIFF.FIFF_BLOCK_END && t_matRawBuffer == FIFF.FIFFB_RAW_DATA)
        t_bIsRunning = false;
    end
    
    t_iCount = t_iCount+1;
    if(t_iCount >= p_nBuffers)
        t_cmdClient.stopAll();
        t_bIsRunning = false;
    end 
end

%% close the sockets
t_cmdClient.close();
t_dataClient.close();

end