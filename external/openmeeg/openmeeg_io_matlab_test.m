% =============================
% = Testing standard matrices =
% =============================

data = randn(4,9);

om_save_full(data,'test.txt','ascii');
om_save_full(data,'test.bin','binary');
om_save_full(data,'test.mat','matlab');

data_txt = om_load_full('test.txt','ascii');
data_bin = om_load_full('test.bin','binary');
data_mat = om_load_full('test.mat','matlab');

norm(data_txt - data)
norm(data_bin - data)
norm(data_mat - data)

delete 'test.txt'
delete 'test.bin'
delete 'test.mat'

% =============================
% = Testing symmetric matrices =
% =============================

randn('seed',0);
data = randn(5,5);
data = (data+data')/2;

om_save_sym(data,'test.txt','ascii');
om_save_sym(data,'test.bin','binary');
om_save_sym(data,'test.mat','matlab');

data_txt = om_load_sym('test.txt','ascii');
data_bin = om_load_sym('test.bin','binary');
data_mat = om_load_sym('test.mat','matlab');

norm(data_txt - data)
norm(data_bin - data)
norm(data_mat - data)

delete 'test.txt'
delete 'test.bin'
delete 'test.mat'

% =============================
% = Testing sparse matrices =
% =============================

data = sprand(5,5,0.5);

om_save_sparse(data,'test.txt','ascii');
om_save_sparse(data,'test.bin','binary');
om_save_sparse(data,'test.mat','matlab');

data_txt = om_load_sparse('test.txt','ascii');
data_bin = om_load_sparse('test.bin','binary');
data_mat = om_load_sparse('test.mat','matlab');

norm(full(data_txt - data))
norm(full(data_bin - data))
norm(full(data_mat - data))

delete 'test.txt'
delete 'test.bin'
delete 'test.mat'

