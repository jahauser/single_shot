#!/usr/bin/env bash
set -euo pipefail
mkdir -p output logs
exec julia --project=. jobs/run.jl "$@"