classdef mne_rt_client < handle
    %MNE_RT_CLIENT is a base class for a mne_rt_server real-time connection
    %
    %   This class provides methods for establishing a mne_rt_server
    %   connection. It stores the established connection until connection
    %   is closed.
    %
    %   Author : Christoph Dinh, Matti Hamalainen, MGH Martinos Center
    %   License : BSD 3-clause
    %
    
    properties (Access = public)
        m_DataInputStream = [];
        m_DataOutputStream = [];
    end
    
    properties (Access = private)
        m_TcpSocket = [];
        m_InputStream = [];
        m_OutputStream = [];
        m_BufferedInputStream = []; %to increase performance use a buffer in between
        m_iNumOfRetries = -1;
    end
    
    methods
        % =================================================================
        %% mne_rt_client
        function obj = mne_rt_client(host, port, number_of_retries)
            %
            % obj = mne_rt_client(host, port, numOfRetries)
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
            obj.init(host, port, number_of_retries);
        end % mne_rt_client
        
        % =================================================================
        %% init
        function result = init(obj, host, port, numOfRetries) 
            %
            % result = init(obj, host, port, numOfRetries) 
            %
            % Init the connection
            %
            % host          - ip adress of the mne_rt_Server
            % port          - comand port of the mne_rt_Server
            % numOfRetries  - number of connection retries, when a
            %                 connection fails

            %   Author : Christoph Dinh, Matti Hamalainen, MGH Martinos Center
            %   License : BSD 3-clause
            %
            import java.net.Socket
            import java.io.*

            if (nargin < 3)
                obj.m_iNumOfRetries = 20; % set to -1 for infinite
            end

            retry = 0;
            obj.close();
            
            result      = false;

            while true

                retry = retry + 1;
                if ((obj.m_iNumOfRetries > 0) && (retry > obj.m_iNumOfRetries))
                    fprintf(1, 'Too many retries\n');
                    break;
                end

                try
                    fprintf(1, 'Retry %d connecting to %s:%d\n', ...
                            retry, host, port);

                    % throws if unable to connect
                    obj.m_TcpSocket = Socket(host, port);

                    obj.m_InputStream = obj.m_TcpSocket.getInputStream;
                    obj.m_BufferedInputStream = BufferedInputStream(obj.m_InputStream);
                    obj.m_DataInputStream = DataInputStream(obj.m_BufferedInputStream);

                    obj.m_OutputStream = obj.m_TcpSocket.getOutputStream;
                    obj.m_DataOutputStream = DataOutputStream(obj.m_OutputStream);
                    
                    fprintf(1, 'Connected to server\n');
                    
                    result = true;
                    
                    break;

                catch
                    obj.close();

                    % pause before retrying
                    pause(1);
                end
            end
        end % init

        % =================================================================
        %% close
        function close(obj)
            %
            % Close the connection
            %

            %   Author : Christoph Dinh, Matti Hamalainen, MGH Martinos Center
            %   License : BSD 3-clause
            %
            if ~isempty(obj.m_TcpSocket)
                obj.m_TcpSocket.close();
            end
            obj.m_TcpSocket = [];
            obj.m_InputStream = [];
            obj.m_OutputStream = [];
            obj.m_BufferedInputStream = []; 
            obj.m_DataInputStream = [];
            obj.m_DataOutputStream = [];
        end % close
    end
end


