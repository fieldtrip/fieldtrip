function testft
% start without a buffer
ft_destroy_buffer
% number of attempts
cnt = 1; 
% use default url
url = 'buffer://localhost:1972';
% dummy header
hdr.Fs = 256;
hdr.nChans = 100;
% after a few attempts, Matlab crashes
while(1)    
    disp(['counter: ' num2str(cnt)]); cnt = cnt + 1;
    ft_create_buffer(1972);
    ft_write_data(url,rand(hdr.nChans,10000), 'header', hdr, 'dataformat', 'fcdc_buffer', 'append', false);
    ft_read_header(url)
    ft_destroy_buffer
    pause(0.1); % probably meaningless
end
