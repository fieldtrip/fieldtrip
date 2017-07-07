if strcmp(computer, 'PCWIN')
  mex midiIn.c  -I../portmidi/pm_common libportmidi_win32.a C:\MinGW\lib\libwinmm.a
  mex midiOut.c -I../portmidi/pm_common libportmidi_win32.a C:\MinGW\lib\libwinmm.a
elseif strcmp(computer, 'PCWIN64')
  mex midiIn.c  -I../portmidi/pm_common libportmidi_win64.a C:\MinGW64\x86_64-w64-mingw32\lib\libwinmm.a
  mex midiOut.c -I../portmidi/pm_common libportmidi_win64.a C:\MinGW64\x86_64-w64-mingw32\lib\libwinmm.a
elseif strcmp(computer, 'GLNX86')
  mex midiIn.c  -I../portmidi/pm_common ../portmidi/libportmidi_linux32.a
  mex midiOut.c -I../portmidi/pm_common ../portmidi/libportmidi_linux32.a
elseif strcmp(computer, 'GLNXA64')
  mex midiIn.c  -I../portmidi/pm_common ../portmidi/libportmidi_linux64.a
  mex midiOut.c -I../portmidi/pm_common ../portmidi/libportmidi_linux64.a
elseif strcmp(computer, 'MACI64')
  % portmidi can be installed from MacPorts with 'sudo port install portmidi' or
  % downloaded and compiled by hand.
  %
  % compiling the mex file requires linking to a number of frameworks. In principle
  % LDFLAGS="$LDFLAGS -framework CoreFoundation -framework CoreAudio -framework CoreMidi -framework Carbon'
  % should copy the existing flags over, but it does not.
  %
  % using the -v option I figured out the difference on the command line with and without the LDFLAG option
  % '-arch x86_64 -Wl,-syslibroot,/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk -mmacosx-version-min=10.7 -bundle -Wl,-exported_symbols_list,/Applications/MATLAB_R2014b.app/extern/lib/maci64/mexFunction.map'
  mex midiIn.c  -I/opt/local/include -L/opt/local/lib -lportmidi_s LDFLAGS='-framework CoreFoundation -framework CoreAudio -framework CoreMidi -framework Carbon -arch x86_64 -Wl,-syslibroot,/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk -mmacosx-version-min=10.7 -bundle -Wl,-exported_symbols_list,/Applications/MATLAB_R2014b.app/extern/lib/maci64/mexFunction.map'
  mex midiOut.c -I/opt/local/include -L/opt/local/lib -lportmidi_s LDFLAGS='-framework CoreFoundation -framework CoreAudio -framework CoreMidi -framework Carbon -arch x86_64 -Wl,-syslibroot,/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk -mmacosx-version-min=10.7 -bundle -Wl,-exported_symbols_list,/Applications/MATLAB_R2014b.app/extern/lib/maci64/mexFunction.map'
else
  ft_error('Unsupported platform');
end
