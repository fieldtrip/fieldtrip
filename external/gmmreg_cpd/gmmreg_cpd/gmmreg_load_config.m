function [config] = gmmreg_load_config(f_config)

model_file = ml_GetPrivateProfileString('Common','model', f_config);
scene_file = ml_GetPrivateProfileString('Common','scene', f_config);
ctrl_pts_file = ml_GetPrivateProfileString('Common','ctrl_pts', f_config);

config.model = load(model_file);
config.scene = load(scene_file);
config.ctrl_pts = load(ctrl_pts_file);

s_sigma = ml_GetPrivateProfileString('gmmreg_cpd_tps_grbf','sigma', f_config);
s_outliers = ml_GetPrivateProfileString('gmmreg_cpd_tps_grbf','outliers', f_config);
s_lambda = ml_GetPrivateProfileString('gmmreg_cpd_tps_grbf','lambda', f_config);
s_anneal_rate = ml_GetPrivateProfileString('gmmreg_cpd_tps_grbf','anneal', f_config);
s_tol = ml_GetPrivateProfileString('gmmreg_cpd_tps_grbf','tol', f_config);
s_emtol = ml_GetPrivateProfileString('gmmreg_cpd_tps_grbf','emtol', f_config);
s_max_iter = ml_GetPrivateProfileString('gmmreg_cpd_tps_grbf','max_iter', f_config);
s_max_em_iter = ml_GetPrivateProfileString('gmmreg_cpd_tps_grbf','max_em_iter', f_config);
s_beta = ml_GetPrivateProfileString('gmmreg_cpd_tps_grbf','beta', f_config);

config.init_sigma = str2num(s_sigma);
config.outliers = str2num(s_outliers);
config.lambda = str2num(s_lambda);
config.beta = str2num(s_beta);
config.anneal_rate = str2num(s_anneal_rate);
config.tol = str2num(s_tol);
config.emtol = str2num(s_emtol);
config.max_iter = str2num(s_max_iter);
config.max_em_iter = str2num(s_max_em_iter);
