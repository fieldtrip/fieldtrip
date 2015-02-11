function menu_about(handle, eventdata, varargin)

% this is a call-back for the FieldTrip menu

msgbox( {
  'FieldTrip is the MATLAB toolbox for MEG and EEG analysis that is being developed at the Centre for Cognitive Neuroimaging of the Donders Institute for Brain, Cognition and Behaviour together with collaborating institutes. The FieldTrip software is released as open source under the GNU general public license.'
  ''
  sprintf('This is FieldTrip, version %s', ft_version)
  } );
