#!/usr/bin/env bash
# Double-click in Finder to start the MovingBoundaryMinerals.jl GUI (macOS).
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"
julia -t auto --project="$DIR" "$DIR/app.jl"
