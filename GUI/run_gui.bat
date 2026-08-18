@echo off
rem Double-click to start the MovingBoundaryMinerals.jl GUI (Windows).
julia -t auto --project=%~dp0 %~dp0app.jl
pause
