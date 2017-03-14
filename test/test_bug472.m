function test_bug472

% MEM 1500mb
% WALLTIME 00:13:01

% TEST buffer.mexa64 buffer.mexmaci buffer.mexw64 buffer.mexglx buffer.mexmaci64 buffer.mexmac buffer.mexw32

global ft_default;
ft_default.feedback = 'no';

% start without a buffer
ft_destroy_buffer

% number of attempts
cnt = 1; 

% use default url
url = 'buffer://localhost:1972';

% dummy header
hdr.Fs = 256;
hdr.nChans = 100;

stopwatch = tic;
% run the test for 10 minutes
% after a few attempts, MATLAB crashes
while (toc(stopwatch)<600)
    disp(['counter: ' num2str(cnt)]); cnt = cnt + 1;
    ft_create_buffer(1972);
    ft_write_data(url,rand(hdr.nChans,10000), 'header', hdr, 'dataformat', 'fcdc_buffer', 'append', false);
    ft_read_header(url)
    ft_destroy_buffer
    pause(0.1); % probably meaningless
end

disp('the test completed without detected problems');

