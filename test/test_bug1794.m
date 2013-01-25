function test_bug1794

return

cd('/home/electromag/johzum/fieldtrip-dev/test')

load bug1794cfg cfg; 
% find attached file bugcfg containing a cfg structure with a
% precomputed grid,  ready for input in ft_prepare_leadfield

cfgmm=cfg;
cfgmm.grid=ft_convert_units(cfg.grid,'mm');
cfgmm.vol=ft_convert_units(cfg.vol,'mm');
cfgmm.grad=ft_convert_units(cfg.grad,'mm');

cfgcm=cfg;
cfgcm.grid=ft_convert_units(cfg.grid,'cm');
cfgcm.vol=ft_convert_units(cfg.vol,'cm');
cfgcm.grad=ft_convert_units(cfg.grad,'cm');

cfgmix1=cfg;
cfgmix1.grid=ft_convert_units(cfg.grid,'cm');
cfgmix1.vol=ft_convert_units(cfg.vol,'mm');
cfgmix1.grad=ft_convert_units(cfg.grad,'mm');

cfgmix2=cfg;
cfgmix2.grid=ft_convert_units(cfg.grid,'mm');
cfgmix2.vol=ft_convert_units(cfg.vol,'cm');
cfgmix2.grad=ft_convert_units(cfg.grad,'cm');

% % first using current (today 20130123)
% cd('/home/common/matlab/fieldtrip')
% restoredefaultpath;
% ft_defaults;
gridLFmmcur= ft_prepare_leadfield(cfgmm);
gridLFcmcur= ft_prepare_leadfield(cfgcm);
gridLFmix1= ft_prepare_leadfield(cfgmix1);
gridLFmix2= ft_prepare_leadfield(cfgmix2);

gridLFmmcur.leadfield{gridLFmmcur.inside(1)}(1,:)
gridLFcmcur.leadfield{gridLFcmcur.inside(1)}(1,:)
gridLFmix1.leadfield{gridLFmix1.inside(1)}(1,:)
gridLFmix2.leadfield{gridLFmix2.inside(1)}(1,:)


%% Testing different versions
cd('/home/common/matlab/fieldtrip-20120630')
restoredefaultpath;
ft_defaults;
gridLFmmjune= ft_prepare_leadfield(cfgmm);
gridLFcmjune= ft_prepare_leadfield(cfgcm);


% Here Add fieltrip-20120426 directory

cd('/home/common/matlab/fieldtrip')
restoredefaultpath;
clear ft_defaults
ft_defaults;
gridLFcmcur= ft_prepare_leadfield(cfg);


gridLFcmcur.leadfield{gridLFcmcur.inside(1)}(1,:)
gridLFmmcur.leadfield{gridLFmmcur.inside(1)}(1,:)
gridLFcmjune.leadfield{gridLFcmjune.inside(1)}(1,:)
gridLFmmjune.leadfield{gridLFmmjune.inside(1)}(1,:)



% % Here Add fieltrip-20120426 directory
% ft_defaults;
% gridLFold= ft_prepare_leadfield(cfg);
% restoredefaultpath;
% % Here Add fieltrip-20121025 directory
% ft_defaults;
% cfg.grid=ft_convert_units(cfg.grid,'cm');
% cfg.vol=ft_convert_units(cfg.vol,'cm');
% cfg.grad=ft_convert_units(cfg.grad,'cm');
% gridLFnew= ft_prepare_leadfield(cfg);

