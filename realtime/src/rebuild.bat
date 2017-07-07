:: Simple bat file to recompile on Windows
:: use "make" on Cygwin, "mingw32-make" on MinGW

@echo off
SET MAKE=mingw32-make --quiet %1 %2 %3 %4 %5 %6 %7 %8 %9

echo Building buffer and ODM...

PUSHD buffer\src
%MAKE% || goto FAILED
POPD

PUSHD buffer\cpp
%MAKE% || goto FAILED
POPD

ECHO Building test

PUSHD buffer\test
%MAKE% || goto FAILED
POPD

FOR /D %%G IN ("utilities\*") DO (
  ECHO Building %%G...
  PUSHD %%G
  %MAKE%
  POPD
  )

:: Blacklisting seems more trouble than it is worth... so keep an eye on the output of make.
FOR /D %%G IN ("acquisition\*") DO (
  ECHO Building %%G...
  PUSHD %%G
  %MAKE%
  POPD
  )

ECHO Compilation of acquisition software complete.

goto EOF
:failed
ECHO Building of libraries or utilities failed :/
POPD

:EOF
