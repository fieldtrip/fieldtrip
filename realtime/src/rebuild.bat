:: Simple bat file to recompile on Windows usign mingw32.
:: TODO: add checking of output
@echo off
SET MAKE=mingw32-make --quiet %1 %2 %3 %4 %5 %6 %7 %8 %9

echo Building buffer and ODM...
PUSHD buffer\src
%MAKE% || goto FAILED
POPD

PUSHD buffer\cpp
%MAKE% || goto FAILED
POPD

ECHO Building utilities
PUSHD utilities\buffer || goto FAILED
%MAKE%
POPD

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