function [ signal, states, parameters, total_samples ] = load_bcidat( varargin )
%LOAD_BCIDAT Load BCI2000 data files into Matlab workspace variables.
%
%  [ signal, states, parameters, total_samples ] = load_bcidat( 'filename1', 'filename2', ... )
%
%  loads signal, state, and parameter data from the files whose names are given
%  as function arguments.
%
%  Examples for loading multiple files:
%    files = dir( '*.dat' );
%    [ signal, states, parameters ] = load_bcidat( files.name );
%
%    files = struct( 'name', uigetfile( 'MultiSelect', 'on' ) );
%    [ signal, states, parameters ] = load_bcidat( files.name );
%
%
%  For multiple files, number of channels, states, and signal type must be
%  consistent.
%
%  By default, signal data will be in raw A/D units, and will be represented by the 
%  smallest Matlab data type that accommodates them.
%  To obtain signal data calibrated into physical units (microvolts),
%  specify '-calibrated' as an option anywhere in the argument list.
%
%  The 'states' output variable will be a Matlab struct with BCI2000 state
%  names as struct member names, and the number of state value entries matching
%  the first dimension of the 'signal' output variable.
%
%  The 'parameters' output variable will be a Matlab struct with BCI2000
%  parameter names as struct member names.
%  Individual parameter values are represented as cell arrays of strings in a 
%  'Value' struct member, and additionally as numeric matrices in a 'NumericValue'
%  struct member. When there is no numeric interpretation possible, the 
%  corresponding matrix entry will be NaN. For nested matrices, no NumericValue
%  field is provided.
%  If multiple files are given, parameter values will match the ones contained 
%  in the first file.
%
%  Optionally, sample ranges may be specified for individual files:
%  [ signal, states, parameters ] = load_bcidat( 'filename', [first last] )
%  will load a subset of samples defined by first and last sample index.
%  Specifying [0 0] for an empty sample range allows to read states and 
%  parameters from a file without reading sample data:
%  [ ignored, states, parameters ] = load_bcidat( 'filename', [0 0] );
%
%  The 'total_samples' output variable reports the total number of samples
%  present in all files.
%
%
%  The load_bcidat function is part of the BCI2000 project. 
%  (C) 2000-2008, BCI2000 Project
%  http://www.bci2000.org

%  This is a help file documenting the functionality contained in
%  load_bcimat.mex.
%  $Id$
%
error( 'There is no load_bcidat mex file for your platform available.' );
