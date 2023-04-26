classdef mne_rt_cmd_client < mne_rt_client
    %MNE_RT_CMD_CLIENT is a class to send commands to mne_rt_server using
    %the command port
    %
    %   This class is inherited of mne_rt_client.
    %
    %   Author : Christoph Dinh, Matti Hamalainen, MGH Martinos Center
    %   License : BSD 3-clause
    %
    methods
        % =================================================================
        %% mne_rt_cmd_client
        function obj = mne_rt_cmd_client(host, port, numOfRetries)
            %
            % obj = mne_rt_cmd_client(host, port, numOfRetries)
            %
            % Constructor
            %
            % host          - ip adress of the mne_rt_Server
            % port          - comand port of the mne_rt_Server
            % numOfRetries  - number of connection retries, when a
            %                 connection fails

            %   Author : Christoph Dinh, Matti Hamalainen, MGH Martinos Center
            %   License : BSD 3-clause
            %
            if (nargin < 3)
                numOfRetries = 20; % set to -1 for infinite
            end
            obj = obj@mne_rt_client(host, port, numOfRetries);%Superclass call
        end % mne_rt_cmd_client
        
        % =================================================================
        %% sendCommand
        function [info] = sendCommand(obj,p_sCommand)
            %
            % [info] = sendCommand(obj,p_sCommand)
            %
            % Sends a command to the mne_rt_server command port
            %
            % p_sCommand    - the command which should be send to the
            %                 mne_rt_server over the command port

            %   Author : Christoph Dinh, Matti Hamalainen, MGH Martinos Center
            %   License : BSD 3-clause
            %
            
            import java.net.Socket
            import java.io.*
            
            t_sCommand = sprintf('%s\n',p_sCommand);
            info = [];
            
            if ~isempty(obj.m_DataOutputStream) && ~isempty(obj.m_DataInputStream)
                
                obj.m_DataOutputStream.writeBytes(t_sCommand);
                obj.m_DataOutputStream.flush;

                % read data from the socket - wait a short time first
                pause(0.5);
                bytes_available = obj.m_DataInputStream.available;

                info = zeros(1, bytes_available, 'uint8');
                for i = 1:bytes_available
                    info(i) = obj.m_DataInputStream.readByte;
                end

                info = char(info);              
            end
        end
        
        % =================================================================
        %% requestMeasInfo
        function requestMeasInfo(obj, AliasOrId)
            %
            % requestMeasInfo(obj, AliasOrId)
            %
            % For convinience, creates the measinfo command and requeststo
            % send the measinfo to the specified ID or Alias
            %
            % AliasOrId    - ID or Alias to which the info should be send
            %
            
            %   Author : Christoph Dinh, Matti Hamalainen, MGH Martinos Center
            %   License : BSD 3-clause
            %
            import java.net.Socket
            import java.io.*
            
            if(ischar(AliasOrId))
                command = sprintf('measinfo %s', AliasOrId);
            elseif(isnumeric(AliasOrId))
                command = sprintf('measinfo %d', AliasOrId);
            else
                error('unknown format for AliasOrId');
            end
            
            obj.sendCommand(command);
        end
        
        % =================================================================
        %% requestMeas
        function requestMeas(obj, AliasOrId)
            %
            % requestMeas(obj, AliasOrId)
            %
            % For convinience, creates the measurement command
            %
            % AliasOrId    - ID or Alias to which the measurement should be send
            %
            
            %   Author : Christoph Dinh, Matti Hamalainen, MGH Martinos Center
            %   License : BSD 3-clause
            %
            
            import java.net.Socket
            import java.io.*
            
            if(ischar(AliasOrId))
                command = sprintf('start %s', AliasOrId);
            elseif(isnumeric(AliasOrId))
                command = sprintf('start %d', AliasOrId);
            else
                error('unknown format for AliasOrId');
            end
            
            obj.sendCommand(command);
        end
        
        % =================================================================
        %% stopAll
        function stopAll(obj)
            %
            % stopMeas(obj, AliasOrId)
            %
            % For convinience, stops all measurements
            %
            
            %   Author : Christoph Dinh, Matti Hamalainen, MGH Martinos Center
            %   License : BSD 3-clause
            %
            
            import java.net.Socket
            import java.io.*

            command = sprintf('stop-all');
            
            obj.sendCommand(command);
        end
        
    end
end

