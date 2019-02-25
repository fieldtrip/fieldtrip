config = gmmreg_load_config('./fish_partial.ini');
config.motion = 'grbf';
config.init_param = zeros(25,2);
[fp,fm] = gmmreg_cpd(config);
DisplayPoints(fm,config.scene,2);