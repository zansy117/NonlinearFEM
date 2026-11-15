@echo off
cd /d "%~dp0"
python gui\dynamics.py
if errorlevel 1 pause
