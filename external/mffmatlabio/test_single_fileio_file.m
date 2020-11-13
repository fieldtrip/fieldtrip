% Script to test a single File-IO file

file =  '/Users/arno/GoogleDrive/EGI/Net Station Files from EGI/Unprocessed Continuous/32 channels/NIA_333ms_HCGSN32_test01.mff';
file =  '/Users/arno/GoogleDrive/EGI/Net Station Files from EGI/OtherFilesDavid/01_024 0531 1145_seg_fil_bcr_ave_WITH_AUTONOMOUS.mff';
tmpfile = '/Users/arno/temp/test.mff';
hdr = ft_read_header(file, 'headerformat', 'egi_mff_v3');
evt = ft_read_event( file, 'eventformat', 'egi_mff_v3', 'header', hdr);
data = ft_read_data( file, 'dataformat', 'egi_mff_v3', 'header', hdr);
%evt = ft_read_event( file, 'eventformat', 'egi_mff_v3');
%data = ft_read_data( file, 'dataformat', 'egi_mff_v3');

ft_write_data(tmpfile, data, 'header', hdr, 'event', evt, 'dataformat', 'mff');
