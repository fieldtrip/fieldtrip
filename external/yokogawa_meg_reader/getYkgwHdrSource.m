% Get header of the source information
% [ Yokogawa MEG Reader toolbox for MATLAB ]
%
% brief:
%   This function retrieves information of the sources in the specified file.
%
% usage:
%   source = getYkgwHdrSource(filepath)
%
% arguments:
%   filepath                  : file path
%
% return values:
%  source                     : structure array of analyzed source information.
%                                Note : Sources are arranged in order of estimated time.
%    .type                    : double, Type of source  
%                                DipoleModel             = 1;
%                                DistributedSourceModel  = 2;
%    .time                    : double, Analyzed Time [second] from 1970.1.1    
%    .sample_no               : double, Time sample index of source     
%    .channel_list            : row vector(double), Channel number (0 origin) list which used to estimate   
%    .model                   : structure of conductor model       
%       .type                 : double, Conductor model type    
%                                UNKNOWN_MODEL        = -1;
%                                NO_MODEL             =  0;
%                                SPHERICAL_MODEL      =  1;
%                                LAYERED_MODEL        =  2;
%      If Conductor model type is SPHERICAL_MODEL,                     
%       .cx                   : double, x coordinate of spherical center position on MEG coordinate [meter]     
%       .cy                   : double, y coordinate of spherical center position on MEG coordinate [meter]     
%       .cz                   : double, z coordinate of spherical center position on MEG coordinate [meter]     
%       .radius               : double, radius of spherical conductor on MRI coordinate [meter] 
%      If Conductor model type is LAYERED_MODEL,                       
%       .ax                   : double, Coefficient 'ax' of planar equation 'ax * x + ay * y + az * z = c'      
%       .ay                   : double, Coefficient 'ay' of planar equation 'ax * x + ay * y + az * z = c'      
%       .az                   : double, Coefficient 'az' of planar equation 'ax * x + ay * y + az * z = c'      
%       .c                    : double, Coefficient 'c' of planar equation 'ax * x + ay * y + az * z = c'       
%    .algorithm               : structure of conductor algorithm
%       .magnetic_field_calc  : double, Algorithm of magnetic field calculation 
%                                BiotSavartLaw        = 1;
%                                SarvasLaw            = 2;
%                                MagneticDipoleLaw    = 3;
%       .variable_restraint   : double, Algorithm of  variable restraint        
%                                NoRestraint          = 0;
%                                PositionRestraint    = 1;
%                                DirectionRestraint   = 2;
%                                IntensityRestraint   = 3;
%       .optimization         : double, Algorithm of  optimization      
%                                GradientAlgorithm                = 1;
%                                LeadFieldReconstructionAlgorithm = 2;
%                                ManualSetAlgorithm               = 3;
%                                UserAlgorithm                    = 4;
%    .filter                  : structure of spectral filter setting       
%       .hpf , .lpf           : structure of high-pass / low-pass filter setting  
%          .enable            :   bool, Does this filter enable?        
%          .cutoff_frequency  : double, Cutoff frequency [Hz]   
%          .window_type       : double, Window type     
%                                NoWindow        = 0;
%                                HanningWindow   = 1;
%                                HammingWindow   = 2;
%          .width             : double, Filter width    
%       .bpf, .bef            : structure of band-pass / band-eliminate filter setting  
%          .enable            : bool, Does this filter enable?        
%          .low_frequency     : double, Low frequency [Hz]      
%          .high_frequency    : double, High frequency [Hz]     
%          .window_type       : double, Window type     
%          .width             : double, Filter width    
%       .moveave              : structure of moving average setting        
%          .enable            :   bool, Does this filter enable?        
%          .width             : double, Filter width    
%       .baseadj              : structure of baseline adjustment setting   
%          .enable            :   bool, Does this filter enable?        
%          .type              : double, Type of baseline adjustment     
%                                PretriggerBaselineAdjust     = 0;
%                                PosttriggerBaselineAdjust    = 1;
%                                AllRangeBaselineAdjust       = 2;
%                                ExplicitBaselineAdjust       = 3;
%          .start_time        : double, Start time [millisecond]        
%          .end_time          : double, End time [millisecond]  
%    .gof                     : double, Goodness-of-fit (GOF)   
%    .correlation             : double, Corrlation Coefficiency   
%    .label                   : double, Label       
%    .comment                 : string, Comment 
%    .total_intensity         : double, Total intensity of sources      
%    .dipole_count            : double, Number of dipole sources        
%    .dipole_list             : structure array of dipole sources   
%       .x                    : double, x coordinate of dipole position on MEG coordinate [meter]       
%       .y                    : double, y coordinate of dipole position on MEG coordinate [meter]       
%       .z                    : double, z coordinate of dipole position on MEG coordinate [meter]       
%       .zdir                 : double, Dipole orientation from z-axis [degree] 
%       .xdir                 : double, Dipole orientation from z-axis [degree] 
%       .intensity            : double, Dipole intensity (moment) [Ampere Meter]        
%
% rivision history
%   2 : 2011.03.03 : Structure fields (confidence volume, confidence ratio, reference_no)
%                     which are not used were removed.
%   1 : 2011.02.14 : 1st argument is modified from file ID to file path.
%   0 : 2010.06.24 : first release
% 
% Copyright (C) 2010-2011 Yokogawa Electric Corporation, All Rights Reserved.
