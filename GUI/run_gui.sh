#!/usr/bin/env bash
# Double-click (or run: ./run_gui.sh) to start the MovingBoundaryMinerals.jl GUI.
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"
julia -t auto --project="$DIR" "$DIR/app.jl"
