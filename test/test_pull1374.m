function test_pull1374

# WALLTIME 00:10:00
# MEM 2gb
# DEPENDENCY xsens_mvnx

filename=dccnpath('/home/common/matlab/fieldtrip/data/test/original/motion/xsens/pull1374.mvnx');

ft_read_header(filename);
ft_read_data(filename);
