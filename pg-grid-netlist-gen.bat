@echo off
REM Batch script wrapper for pg-grid-netlist-gen on Windows
REM This allows running the command without activating the virtual environment

set SCRIPT_DIR=%~dp0
set VENV_EXE=%SCRIPT_DIR%.venv\Scripts\pg-grid-netlist-gen.exe

if not exist "%VENV_EXE%" (
    echo Error: Entry point not found at %VENV_EXE%
    echo Please run: uv pip install -e .
    exit /b 1
)

"%VENV_EXE%" %*
