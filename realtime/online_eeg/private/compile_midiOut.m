if strcmp(computer, 'PCWIN')
  mex midiOut.c -I../portmidi/pm_common libportmidi_win32.a C:\MinGW\lib\libwinmm.a
elseif strcmp(computer, 'PCWIN64')
  mex midiOut.c -I../portmidi/pm_common libportmidi_win64.a C:\MinGW64\x86_64-w64-mingw32\lib\libwinmm.a
elseif strcmp(computer, 'GLNX86')
  mex midiOut.c -I../portmidi/pm_common ../portmidi/libportmidi_linux32.a
elseif strcmp(computer, 'GLNXA64')
  mex midiOut.c -I../portmidi/pm_common ../portmidi/libportmidi_linux64.a
elseif strcmp(computer, 'MACI64')
  mex midiOut.c -I/opt/local/include -L/opt/local/lib -lportmidi_s 
else
  error('Unsupported platform');
end