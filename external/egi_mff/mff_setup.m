function mff_setup

% MFF_SETUP ensures that the JVM is available and that the MFF java
% implementation is on the java classpath.

% this toolbox requires JAVA to be running and properly configured
if ~usejava('jvm')
  error('the EGI_MFF format requires that the java virtual machine is available, see http://www.fieldtriptoolbox.org/getting_started/egi')
end

% check that the java archive (jar file) is on the path
% the user can add it to the static path using classpath.txt
if ~(exist('com.egi.services.mff.api.MFFFactory', 'class')==8)
    jarpath = fullfile(fileparts(mfilename('fullpath')), 'java', 'M*.jar');
    jarlist = dir(jarpath);
    if size(jarlist, 1) == 1
        
        % Some MATLAB versions have a bug in javaclasspath that causes global
        % variables to be erased. Hence we keep a local copy and restore it afterwards.
        % See http://bugzilla.fcdonders.nl/show_bug.cgi?id=2133

        var = whos('global');
        glob = [var.global];
        name = {var.name};
        name = name(glob);
        for i=1:length(name)
        eval(sprintf('global %s', name{i}));
        eval(sprintf('localcopy.%s = %s;', name{i}, name{i}));
        end
        
        % try to add it to the dynamic classpath
        jarpath = fullfile(fileparts(jarpath), jarlist(1).name);
        warning('adding %s to your Java path', jarpath);
        javaclasspath(jarpath);
        
        % The following pertains to http://bugzilla.fcdonders.nl/show_bug.cgi?id=2133
        % also available as http://bit.ly/1o7QlRk, for which we cannot provide a solution. 
        warning('In some versions of Matlab, javaclasspath clears all global variables due to a bug in MATLAB. Please execute %s.m immediately after starting MATLAB. See http://bit.ly/1o7QlRk', mfilename('fullpath'));

        % restore the global variables
        for i=1:length(name)
        eval(sprintf('global %s', name{i}));
        eval(sprintf('%s = localcopy.%s;', name{i}, name{i}));
        end
        clear localcopy

    else
        error('the EGI_MFF format requires one MFF .jar file on your classpath, see http://www.fieldtriptoolbox.org//getting_started/egi')
    end
end

% check once more, give an error if it still does not exist
if ~(exist('com.egi.services.mff.api.MFFFactory', 'class')==8)
  error('the EGI_MFF format requires one MFF .jar file on your classpath, see http://www.fieldtriptoolbox.org/getting_started/egi')
end
