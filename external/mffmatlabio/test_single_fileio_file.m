% Script to test a single File-IO file

file =  '/Users/arno/GoogleDrive/EGI/Net Station Files from EGI/Unprocessed Continuous/32 channels/NIA_333ms_HCGSN32_test01.mff';
tmpfile = '/tmp/test.mff';
hdr = ft_read_header(file, 'headerformat', 'egi_mff_v3');
evt = ft_read_event( file, 'eventformat', 'egi_mff_v3', 'header', hrd);
data = ft_read_data( file, 'dataformat', 'egi_mff_v3', 'header', hrd);

ft_write_data(tmpfile, data, 'header', hdr, 'event', evt, 'dataformat', 'mff');
