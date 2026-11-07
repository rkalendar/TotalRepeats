@echo off
REM Compiles TotalRepeats together with its test suite and runs the suite.
REM Requires only a JDK (25+); there are no external test dependencies.
REM
REM   run-tests.bat
REM
REM Exits non-zero if any check fails, so it can gate a build or a CI job.

setlocal
set HERE=%~dp0
set OUT=%HERE%build\test-classes

if exist "%OUT%" rmdir /s /q "%OUT%"
mkdir "%OUT%"

javac -d "%OUT%" "%HERE%src\*.java" "%HERE%tests\*.java"
if errorlevel 1 exit /b 1

java -cp "%OUT%" TestRunner
exit /b %errorlevel%
