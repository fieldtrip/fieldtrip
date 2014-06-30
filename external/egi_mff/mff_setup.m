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
    jarpath = fullfile(fileparts(mfilename('fullpath')), 'java', 'M*.jar');
    jarlist = dir(jarpath);
    if size(jarlist, 1) == 1
        % try to add it to the dynamic classpath
        jarpath = fullfile(fileparts(jarpath), jarlist(1).name);
        warning('adding %s to your Java path', jarpath);
        javaclasspath(jarpath);
    else
        error('the EGI_MFF format requires one MFF .jar file on your classpath, see http://fieldtrip.fcdonders.nl/getting_started/egi')
    end
end

% check once more, give an error if it still does not exist
if ~(exist('com.egi.services.mff.api.MFFFactory', 'class')==8)
  error('the EGI_MFF format requires one MFF .jar file on your classpath, see http://fieldtrip.fcdonders.nl/getting_started/egi')
end
