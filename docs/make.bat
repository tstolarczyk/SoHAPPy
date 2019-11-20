@ECHO OFF

pushd %~dp0

REM Command file for Sphinx documentation
REM sphinx-build -b html -D graphviz_dot=C:\graphviz\bin\dot.exe . _build/html
if "%SPHINXBUILD%" == "" (
set SPHINXBUILD=sphinx-build

)
REM Source of the documentation
set SOURCEDIR=.
set BUILDDIR=build

if "%1" == "" goto help

%SPHINXBUILD% >NUL 2>NUL
if errorlevel 9009 (
	echo.
	echo.The 'sphinx-build' command was not found. Make sure you have Sphinx
	echo.installed, then set the SPHINXBUILD environment variable to point
	echo.to the full path of the 'sphinx-build' executable. Alternatively you
	echo.may add the Sphinx directory to PATH.
	echo.
	echo.If you don't have Sphinx installed, grab it from
	echo.http://sphinx-doc.org/
	exit /b 1
)

REM remove previous file otherwise cannot be overwritten

del %BUILDDIR% 
del api
%SPHINXBUILD% -M %1 %SOURCEDIR% %BUILDDIR% %SPHINXOPTS%
goto end

:help
%SPHINXBUILD% -M help %SOURCEDIR% %BUILDDIR% %SPHINXOPTS%
REM - for inheritance diagrams, which are quite useless in this code
REM %SPHINXBUILD% -b html -D graphviz_dot=C:\Users\stolar\AppData\Local\Continuum\Anaconda3\Library\bin\dot.bat . build/html
:end
popd
