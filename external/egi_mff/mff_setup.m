function mff_setup

% MFF_SETUP ensures that the JVM is available and that the MFF java
% implementation is on the java classpath.

% this toolbox requires JAVA to be running and properly configured
if ~usejava('jvm')
  error('the EGI_MFF format requires that the java virtual machine is available, see http://fieldtrip.fcdonders.nl/getting_started/egi')
end

% check that the java archive (jar file) is on the path
% the user can add it to the static path using classpath.txt
if ~(exist('com.egi.services.mff.api.MFFFactory', 'class')==8)
  % try to add it to the dynamic classpath
  p = fullfile(fileparts(mfilename('fullpath')), 'java', 'MFF-1.0.jar');
  warning('adding %s to your Java path', p);
  javaclasspath(p);
end

% check once more, give an error if it still does not exist
if ~(exist('com.egi.services.mff.api.MFFFactory', 'class')==8)
  error('the EGI_MFF format requires "MFF-1.0.jar" on your classpath, see http://fieldtrip.fcdonders.nl/getting_started/egi')
end

