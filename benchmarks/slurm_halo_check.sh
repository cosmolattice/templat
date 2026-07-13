#!/bin/bash
#SBATCH --job-name=halo_check
#SBATCH --output=halo_check_%j.log
#SBATCH --time=00:30:00
#SBATCH --nodes=2                 # >=2 nodes: only a real network shows the true message cost
#SBATCH --ntasks-per-node=32      # adjust to the cores/node on your partition
#SBATCH --cpus-per-task=1
# #SBATCH --partition=<your_partition>
# #SBATCH --account=<your_account>

# ---------------------------------------------------------------------------
# Halo-exchange scaling check on the cluster.
# Builds bench-halo_scaling once at a chosen grid size, then runs a rank sweep.
# The benchmark itself prints, per run:
#   - kick vs halo_update vs multi-component (cplx_halo, su2_halo) timings
#   - the overlap CEILING (kick/halo geometry, analytic bound only)
# Run me with:  sbatch benchmarks/slurm_halo_check.sh
# ---------------------------------------------------------------------------
set -euo pipefail

# --- environment: load whatever gives you cmake + an MPI + FFTW (edit for your site) ---
# module load cmake gcc openmpi fftw
export OMP_NUM_THREADS=1          # one thread per rank; scaling is over MPI ranks here

REPO="${SLURM_SUBMIT_DIR:-$PWD}"
BUILD="$REPO/build-halo"
NGRID="${NGRID:-256}"            # override at submit: NGRID=512 sbatch ...
NSTEPS="${NSTEPS:-100}"

# --- configure + build (kokkos is cached after the first build; size changes only recompile the bench) ---
cmake -S "$REPO" -B "$BUILD" \
  -DCMAKE_BUILD_TYPE=Release \
  -DTEMPLAT_BENCH=ON -DMPI=ON -DPARAFAFT=ON -DNOTHREADING=ON \
  -DHALO_NGRID="$NGRID" -DHALO_NSTEPS="$NSTEPS"
cmake --build "$BUILD" --target bench-halo_scaling -j 8

BIN="$BUILD/benchmarks/bench-halo_scaling"

# --- rank sweep. IMPORTANT: force the ranks to SPREAD across all allocated nodes, otherwise srun
#     packs them onto one node and you only ever measure the (break-even) single-node case. ---
NNODES="$SLURM_NNODES"
TOTAL=$(( SLURM_NNODES * SLURM_NTASKS_PER_NODE ))
echo "### halo_check: NGRID=$NGRID NSTEPS=$NSTEPS  nodes=$NNODES  total_ranks=$TOTAL"

run_np () {  # run_np <total_ranks> <nodes>
  local NP="$1" NODES="$2"
  (( NP % NODES != 0 )) && NODES=1                 # keep it even; fall back to 1 node
  local PERNODE=$(( NP / NODES ))
  echo ""
  echo "==================== ranks = $NP  (${NODES} node(s), ${PERNODE}/node) ===================="
  srun --nodes="$NODES" --ntasks="$NP" --ntasks-per-node="$PERNODE" --cpu-bind=cores "$BIN"
}

# Single-node baseline vs the same total spread over all nodes — the difference is the network effect.
run_np 1      1
run_np 16     1
run_np "$TOTAL" 1          # all ranks packed on one node (your production config)
run_np "$TOTAL" "$NNODES"  # same ranks spread across nodes -> this is where overlap can pay off

echo ""
echo "### done. Compare the last two runs (1 node vs ${NNODES} nodes): the difference is the network cost."
echo "### cplx_halo/su2_halo vs halo_update shows the per-component cost of multi-component fields."
