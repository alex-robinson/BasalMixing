#!/usr/bin/env bash
# runme executable for the BasalMixing ensemble.
#
# runme stages this script together with the Julia sources, data/, and
# Project.toml into the run directory, then runs it there with the staged TOML
# parameter file as the first argument. `Pkg.instantiate()` below resolves the
# environment from Project.toml (writing a platform-correct Manifest.toml in the
# run directory), then it samples the posterior and writes chain.jld2 there.
#
# Plotting is deliberately NOT done here: figures often overlay two runs (e.g.
# combined vs 81Kr-only), so they live outside any single run. Produce them
# separately from the repo with plot_basalmixing_ensemble.jl (see its header);
# figures go to the repo's plots/ folder.
set -euo pipefail

parfile="${1:-basalmixing.toml}"

# One Julia thread per MCMC chain. On SLURM this comes from --cpus-per-task
# (set via `runme --omp N`, or the `omp` value in .runme/config.toml); off the
# cluster it falls back to OMP_NUM_THREADS, then to 4.
nthreads="${SLURM_CPUS_PER_TASK:-${OMP_NUM_THREADS:-4}}"

echo "=== BasalMixing ensemble: sampling (parfile=${parfile}, threads=${nthreads}) ==="
# --startup-file=no is required: a personal ~/.julia/config/startup.jl runs
# `Pkg.activate(myanalysis)` unconditionally, which otherwise hijacks --project=.
# and makes instantiate resolve the wrong environment (no Manifest in the run dir).
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no -t "${nthreads}" --project=. run_basalmixing_ensemble.jl "${parfile}"

echo "=== Done: chain.jld2 written. Plot separately from the repo. ==="
