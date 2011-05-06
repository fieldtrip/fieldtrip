function sb_write_par(cfg,file_name)

% sb_write_par writes a parameter file for SimBio on disk
%
% Use as
%   sb_write_dip(cfg,file_name)
% 
% file_name is the name of the parameter file (without extension)
% The resuired field are:
%   cfg.cond        conductivities file
%   cfg.labels      labels, an array of integers relative to the tissue
%                    types derived from the segmented MRI
%   cfg.tensorlabel [optional] labels, an array of integers to map
%                   conductivity anisotropy to the different tissues
%                   (not implemented yet)
%
% Defaults conductivities are expressed in S/m according to:
% Oostendorp T, Delbeke J, Stegeman DF,
% The conductivity of the Human Skull: Results of In Vivo and In Vitro
% Measurements. IEEE Biomed. Eng, Vol. 47, No.11, 2000
%   and
% Baumann SB, Wozny DR, Kelly SK, Meno FM,
% The Electrical Conductivity of Human Cerebrospinal Fluid at Body
% Temperature. IEEE Biomed.Eng, Vol. 44, No.3, 1997, pp.220-223
% 
% Default labels, conductivities are as following:
% skull 1 , 0.022
% extra 3
% white 6 , 0.14
% grey  7 , 0.33
% csf   8 , 1.79
% skin  ? , 0.33 
% 
% Copyright (C) 2011, Felix Lucka and Cristiano Micheli

ft_defaults

% # Labels in the head model corresponding to the different tissues
% # The first value corresponds to the electrode which will be added to the mesh
% Here, the labels from the MR image are used to assign the conductivities
% given above to the compartments
cfg.labels = ft_getopt(cfg,'labels');
labels = [1000 cfg.labels(:)'];

% # Conductivities of fem head model
% # The first value corresponds to the electrodes which will be
% # added to the mesh 
cfg.cond   = ft_getopt(cfg,'cond');
cond = [1.0 cfg.cond(:)'];

% # Tissue labels in the head model for which the tensor valued coductivity
% # should be used if available
tenslabels = zeros(1,length(cfg.cond)+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SimBio general parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% #residuum of the forward solution
TOL = 0.1e-7;
opt.toleranceforwardsolution  = TOL;
% # degree of Gaussian integration (2 is recommended)
opt.degreeofintegration       = 2;
% # different analytical solutions for the EEG problem
opt.analyticalsolution        = 1;
% # should METIS repartition for parallel run? (1 = yes?)
opt.metisrepartitioning       = 1;
% # SOLVER (1:Jakobi-CG, 2:IC(0)-CG, 3:AMG-CG, 4:PILUTS(ILDLT)-CG)
opt.solvermethod              = 2;
% # use or use not lead field basis approach (1 = use it?)
opt.associativity             = 1; 
% # NeuroFEM Solver
% # parameter file for Pebbles solver
opt.pebbles                   = 'pebbles.inp'; % only for 3:AMG-CG method
% ONLY for MultiRHS-Pebbles: number of right-hand-sides which are solved simultaneously
opt.pebblesnumrhs             = 1;
% # SOURCE SIMULATION'
% # threshold (percentage of the greatest dipole strength) of all dipoles to appear in the result files
opt.dipolethreshold           = 1;
% # blurring single loads to adjacent nodes by means of the Gaussian (normal) distribution
opt.sourcesimulationepsilondirac= TOL;
% # DIPOLE MODELING
% # dipole modeling; weighting of the source distribution with the power of the declared value
opt.dipolemodelingsmoothness  = 2;
% # power of the dipole moments to be considered
opt.dipolemodelingorder       = 2;
% # necessary internal scaling factor; should be larger than twice the element edge length
opt.dipolemodelingscale       = 20.000;
% # Lagrangian multiplier for the (inverse) dipole modeling
opt.dipolemodelinglambda      = TOL*100;
% # source-sink separation of the analytical dipole model
opt.dipolemodelingdistance    = 1;
% # use rango dipole model (0 = not use)?
opt.dipolemodelingrango       = 0;
% # dipole model (0 = Blurred, 1 = Zenger, 2 = Rango)
opt.dipolemodel               = 0;
% # Monopole/Dipole
% # calculation of the potential distribution for spherical, homogeneous structures 
% # by means of a analytical description (0 = not use?)
opt.analyticaldipole          = 0;
% # forward solution computation with monopoles or dipoles (1 = dipoles)
opt.dipolesources             = 1;
% # spread a nodal load to adjacent nodes (1 = spread yes?)
opt.sourcesimulation          = 1;
% # to compare analytical and numerical solutions, an integral average value
% # of the potential is subtracted from the analytical results because they
% # are not related to a mean potential (in contrast to the numerical solution)
% # (1 = correction yes?) 
opt.averagecorrection         = 0;
% # Number of materials (par #22)
opt.nummaterials              = 7; %numel(cond);
% # Conductivities (par #23)
opt.conductivities            = cond;  
% # Labels (par #24)
opt.labels                    = labels;  
% # Tensor Labels (par #25)
opt.condtensor                = tenslabels;
% #ReferenceData file: 1= Vista, 2= ASA, 3= ASCII
opt.ReferenceDataIn           = 2;
% #LeadfieldMatrix input file: 1= Vista, 2= ASA 
opt.LeadfieldIn               = 1;
% #SensorConfiguration file: 1= Vista, 2= ASA 
opt.SensorconfigurationIn     = 2;
% #Source Space Grid: 1= Vista, 2= ASA, 3= CAUCHY
opt.SourceSpaceIn             = 2;
% #FEM HeadGrid: 1= Vista, 2= ASA, 3= CAUCHY
opt.HeadGridIn                = 1;
% #Output File; Format: 1= Vista, 2 = ASA, 3 = Vista + ASA, 4 = Vista +
% # Ascii, 5 = Cauchy, 6 = SCIRun 
opt.ResultOut                 = 2;

% fill in the options cell
fnames = fieldnames(opt);
for ii=1:numel(fnames)
  options{ii} = getfield(opt,fnames{ii});
end

str = sprintf(['#Parameter file: FEM for source simulation\n\n', ...
  '[NeuroFEMSimulator]\n', ...
  'toleranceforwardsolution= %0.6E\n', ...
  'degreeofintegration= %d\n', ...
  'analyticalsolution= %d\n\n', ...
  'metisrepartitioning= %d\n\n', ...
  'solvermethod= %d\n', ...
  'associativity= %d\n', ...
  'pebbles= %s\n', ...
  'pebblesnumrhs= %d\n\n', ...
  'dipolethreshold= %0.6E\n', ...
  'sourcesimulationepsilondirac= %0.5E\n\n', ...
  'dipolemodelingsmoothness= %d\n', ...
  'dipolemodelingorder= %d\n', ...
  'dipolemodelingscale= %0.3f\n', ...
  'dipolemodelinglambda= %0.5E\n', ...
  'dipolemodelingdistance= %0.4f\n', ...
  'dipolemodelingrango= %d\n', ...
  'dipolemodel= %d\n\n', ...
  'analyticaldipole= %d\n', ...
  'dipolesources= %d\n', ...
  'sourcesimulation= %d\n', ...
  'averagecorrection= %d\n\n', ...
  '[NeuroFEMGridGenerator]\n', ...
  'nummaterials= %d\n', ...
  'conductivities= \n', ...
  '%s\n', ... 
  'labels=\n', ... 
  '%s\n\n', ... 
  'tensorlabels=\n', ... 
  '%s\n\n', ...
  '[FileFormats]   \n', ... 
  'ReferenceDataIn= %d\n', ... 
  'LeadfieldIn= %d\n', ... 
  'SensorconfigurationIn= %d \n', ... 
  'SourceSpaceIn= %d\n', ... 
  'HeadGridIn= %d\n', ... 
  'ResultOut= %d\n\n'], ...
  options{1},options{2},options{3},options{4},options{5},options{6},options{7},options{8}, ...
  options{9},options{10},options{11},options{12},options{13},options{14},options{15},options{16}, ...
  options{17},options{18},options{19},options{20},options{21},options{22},num2str(options{23}), ...
  num2str(options{24}),num2str(options{25}), ...  
  options{26},options{27},options{28},options{29},options{30},options{31}); 

  % write the par file on disk
  fid = fopen(file_name,'w');
  fprintf(fid,str);
  fclose(fid);

