@echo off
rem ==========================================================================================
rem  ConformalStereo - double-click this file to start the program.
rem
rem  The program needs Python 3.8 or newer and nothing else.
rem    1. If Python is already on this computer, it is used.
rem    2. If not, a private copy of Python is fetched from python.org into the folder
rem       "runtime" next to this file (about 11 MB, one time only).
rem  Nothing is installed on the computer and no administrator rights are needed.
rem  To remove everything, delete this folder.
rem
rem  For testing: "ConformalStereo.bat --private-python" skips step 1.
rem ==========================================================================================
setlocal EnableExtensions DisableDelayedExpansion
title ConformalStereo
cd /d "%~dp0"

set "CS_APP=%~dp0app\run.py"
set "CS_RUNTIME=%~dp0runtime"
set "CS_PYDIR=%~dp0runtime\python"
set "CS_PYVER=3.12.10"
set "CS_PY="
set "CS_ARGS="

if not exist "%CS_APP%" goto broken

if /i "%~1"=="--private-python" goto private
set CS_ARGS=%*

rem ---- a private copy fetched on an earlier start
if not exist "%CS_PYDIR%\python.exe" goto installed
call :try "%CS_PYDIR%\python.exe"
if defined CS_PY goto run

:installed
rem ---- a Python that is already on this computer (the test also rejects the Microsoft Store placeholder)
call :try py -3
if defined CS_PY goto run
call :try python
if defined CS_PY goto run
call :try python3
if defined CS_PY goto run

:private
if not exist "%CS_PYDIR%\python.exe" goto fetch
call :try "%CS_PYDIR%\python.exe"
if defined CS_PY goto run

:fetch
echo.
echo   ConformalStereo needs Python, and none was found on this computer.
echo   A private copy is fetched now from python.org, about 11 MB, one time only.
echo   It goes into the folder "runtime" next to this file. Nothing is installed.
echo.
set "CS_ARCH=amd64"
if /i "%PROCESSOR_ARCHITECTURE%"=="ARM64" set "CS_ARCH=arm64"
if /i "%PROCESSOR_ARCHITECTURE%"=="x86" set "CS_ARCH=win32"
if /i "%PROCESSOR_ARCHITEW6432%"=="AMD64" set "CS_ARCH=amd64"
if /i "%PROCESSOR_ARCHITEW6432%"=="ARM64" set "CS_ARCH=arm64"
set "CS_URL=https://www.python.org/ftp/python/%CS_PYVER%/python-%CS_PYVER%-embed-%CS_ARCH%.zip"
set "CS_ZIP=%~dp0runtime\python-embed.zip"
if not exist "%CS_RUNTIME%" mkdir "%CS_RUNTIME%"
if not exist "%CS_RUNTIME%" goto nowrite
if exist "%CS_ZIP%" del "%CS_ZIP%"
echo   Downloading %CS_URL%
powershell -NoProfile -ExecutionPolicy Bypass -Command "$ErrorActionPreference = 'Stop'; $ProgressPreference = 'SilentlyContinue'; try { [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; Invoke-WebRequest -Uri $env:CS_URL -OutFile $env:CS_ZIP -UseBasicParsing } catch { exit 1 }"
if not exist "%CS_ZIP%" curl.exe --location --fail --silent --show-error --output "%CS_ZIP%" "%CS_URL%"
if not exist "%CS_ZIP%" goto fetchfailed
echo   Unpacking ...
powershell -NoProfile -ExecutionPolicy Bypass -Command "$ErrorActionPreference = 'Stop'; try { if (Test-Path -LiteralPath $env:CS_PYDIR) { Remove-Item -LiteralPath $env:CS_PYDIR -Recurse -Force }; Expand-Archive -LiteralPath $env:CS_ZIP -DestinationPath $env:CS_PYDIR -Force } catch { exit 1 }"
if exist "%CS_PYDIR%\python.exe" goto unpacked
if not exist "%CS_PYDIR%" mkdir "%CS_PYDIR%"
tar.exe -xf "%CS_ZIP%" -C "%CS_PYDIR%"
:unpacked
del "%CS_ZIP%" >nul 2>nul
if not exist "%CS_PYDIR%\python.exe" goto fetchfailed
call :try "%CS_PYDIR%\python.exe"
if not defined CS_PY goto fetchfailed
echo   Python %CS_PYVER% is ready.
goto run

:run
%CS_PY% "%CS_APP%" %CS_ARGS%
if errorlevel 1 goto failed
goto end

:try
rem  %* is a command that should start Python. It is accepted if it is version 3.8 or newer
rem  and can load the parts of Python that the program uses.
%* -c "import sys, json, zipfile, http.server, webbrowser, xml.etree.ElementTree; sys.exit(0 if sys.version_info >= (3, 8) else 1)" >nul 2>nul
if errorlevel 1 goto :eof
set CS_PY=%*
goto :eof

:fetchfailed
echo.
echo   Python could not be fetched. Check the internet connection and start again.
echo   If this computer may not reach python.org, install Python from
echo   https://www.python.org/downloads/ and start this file again.
echo.
pause
exit /b 1

:nowrite
echo.
echo   This folder cannot be written to. Move the folder ConformalStereo_App to a place
echo   of your own, for example the Desktop or Documents, and start again.
echo.
pause
exit /b 1

:broken
echo.
echo   The folder "app" is missing next to this file.
echo   Unpack the whole folder ConformalStereo_App first, then start this file from there.
echo.
pause
exit /b 1

:failed
echo.
echo   ConformalStereo stopped with an error. The lines above say what happened.
echo.
pause
exit /b 1

:end
endlocal
exit /b 0
