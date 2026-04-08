// check_gpu_mpi.cu — Diagnose GPU-aware MPI and GPUDirect RDMA availability.
//
// Compile:
//   nvcc -o check_gpu_mpi check_gpu_mpi.cu -lmpi -I$(mpi_include_path)
//   or: mpicxx -x cu -o check_gpu_mpi check_gpu_mpi.cu
//
// Run (needs exactly 2 ranks):
//   mpirun -np 2 ./check_gpu_mpi
//
// What it does:
//   1. Tries MPI_Send/Recv with a device pointer. If MPI segfaults or errors,
//      GPU-aware MPI is not available.
//   2. If step 1 works, benchmarks device-pointer transfers vs explicit
//      device→host copy + host transfer to detect GPUDirect RDMA.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <vector>
#include <numeric>
#include <algorithm>

#include <mpi.h>
#include <cuda_runtime.h>

#define CUDA_CHECK(call)                                                       \
  do {                                                                         \
    cudaError_t err = (call);                                                  \
    if (err != cudaSuccess) {                                                  \
      fprintf(stderr, "[rank %d] CUDA error at %s:%d: %s\n", rank, __FILE__,  \
              __LINE__, cudaGetErrorString(err));                              \
      MPI_Abort(MPI_COMM_WORLD, 1);                                           \
    }                                                                          \
  } while (0)

static int rank, nranks;

// ---------------------------------------------------------------------------
// Ping-pong: returns median one-way time in microseconds.
// If use_device_ptr is true, MPI touches the device pointer directly.
// Otherwise we manually stage through a host buffer.
// ---------------------------------------------------------------------------
double ping_pong(size_t nbytes, int niters, bool use_device_ptr) {
  void *d_buf = nullptr;
  void *h_buf = nullptr;

  CUDA_CHECK(cudaMalloc(&d_buf, nbytes));
  CUDA_CHECK(cudaMemset(d_buf, rank + 1, nbytes));
  h_buf = malloc(nbytes);

  // Warm up
  for (int i = 0; i < 5; ++i) {
    if (rank == 0) {
      if (use_device_ptr) {
        MPI_Send(d_buf, nbytes, MPI_BYTE, 1, 0, MPI_COMM_WORLD);
        MPI_Recv(d_buf, nbytes, MPI_BYTE, 1, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
      } else {
        CUDA_CHECK(cudaMemcpy(h_buf, d_buf, nbytes, cudaMemcpyDeviceToHost));
        MPI_Send(h_buf, nbytes, MPI_BYTE, 1, 0, MPI_COMM_WORLD);
        MPI_Recv(h_buf, nbytes, MPI_BYTE, 1, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        CUDA_CHECK(cudaMemcpy(d_buf, h_buf, nbytes, cudaMemcpyHostToDevice));
      }
    } else {
      if (use_device_ptr) {
        MPI_Recv(d_buf, nbytes, MPI_BYTE, 0, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        MPI_Send(d_buf, nbytes, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
      } else {
        MPI_Recv(h_buf, nbytes, MPI_BYTE, 0, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        CUDA_CHECK(cudaMemcpy(d_buf, h_buf, nbytes, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(h_buf, d_buf, nbytes, cudaMemcpyDeviceToHost));
        MPI_Send(h_buf, nbytes, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
      }
    }
  }

  CUDA_CHECK(cudaDeviceSynchronize());
  MPI_Barrier(MPI_COMM_WORLD);

  // Timed iterations
  std::vector<double> times(niters);
  for (int i = 0; i < niters; ++i) {
    MPI_Barrier(MPI_COMM_WORLD);
    auto t0 = std::chrono::high_resolution_clock::now();

    if (rank == 0) {
      if (use_device_ptr) {
        MPI_Send(d_buf, nbytes, MPI_BYTE, 1, 0, MPI_COMM_WORLD);
        MPI_Recv(d_buf, nbytes, MPI_BYTE, 1, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
      } else {
        CUDA_CHECK(cudaMemcpy(h_buf, d_buf, nbytes, cudaMemcpyDeviceToHost));
        MPI_Send(h_buf, nbytes, MPI_BYTE, 1, 0, MPI_COMM_WORLD);
        MPI_Recv(h_buf, nbytes, MPI_BYTE, 1, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        CUDA_CHECK(cudaMemcpy(d_buf, h_buf, nbytes, cudaMemcpyHostToDevice));
      }
    } else {
      if (use_device_ptr) {
        MPI_Recv(d_buf, nbytes, MPI_BYTE, 0, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        MPI_Send(d_buf, nbytes, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
      } else {
        MPI_Recv(h_buf, nbytes, MPI_BYTE, 0, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        CUDA_CHECK(cudaMemcpy(d_buf, h_buf, nbytes, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(h_buf, d_buf, nbytes, cudaMemcpyDeviceToHost));
        MPI_Send(h_buf, nbytes, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
      }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    // Half round-trip = one-way latency
    times[i] =
        std::chrono::duration<double, std::micro>(t1 - t0).count() / 2.0;
  }

  std::sort(times.begin(), times.end());
  double median = times[niters / 2];

  CUDA_CHECK(cudaFree(d_buf));
  free(h_buf);
  return median;
}

// ---------------------------------------------------------------------------
// Step 1: Can MPI handle a device pointer at all?
// We fork a child so that a segfault doesn't kill the main process.
// ---------------------------------------------------------------------------
bool test_gpu_aware_mpi() {
  // Simple correctness test: send a known pattern via device pointer.
  const int N = 256;
  int *d_buf = nullptr;
  std::vector<int> h_buf(N);

  CUDA_CHECK(cudaMalloc(&d_buf, N * sizeof(int)));

  if (rank == 0) {
    // Fill device buffer with pattern
    std::vector<int> pattern(N);
    for (int i = 0; i < N; ++i) pattern[i] = i * 42 + 7;
    CUDA_CHECK(
        cudaMemcpy(d_buf, pattern.data(), N * sizeof(int), cudaMemcpyHostToDevice));

    MPI_Send(d_buf, N * sizeof(int), MPI_BYTE, 1, 99, MPI_COMM_WORLD);
  } else {
    CUDA_CHECK(cudaMemset(d_buf, 0, N * sizeof(int)));
    MPI_Recv(d_buf, N * sizeof(int), MPI_BYTE, 0, 99, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);

    // Copy back and verify
    CUDA_CHECK(
        cudaMemcpy(h_buf.data(), d_buf, N * sizeof(int), cudaMemcpyDeviceToHost));
    for (int i = 0; i < N; ++i) {
      if (h_buf[i] != i * 42 + 7) {
        CUDA_CHECK(cudaFree(d_buf));
        return false;
      }
    }
  }

  CUDA_CHECK(cudaFree(d_buf));
  return true;
}

// ---------------------------------------------------------------------------
// Query CUDA-aware MPI support at compile/runtime if the MPI implementation
// exposes the standard MPIX query (OpenMPI 2.0+, MVAPICH2, Cray MPICH).
// ---------------------------------------------------------------------------
void print_compile_time_hints() {
#if defined(MPIX_CUDA_AWARE_SUPPORT)
  if (rank == 0)
    printf("  MPIX_CUDA_AWARE_SUPPORT compile-time flag: %s\n",
           MPIX_CUDA_AWARE_SUPPORT ? "YES" : "NO");
#else
  if (rank == 0)
    printf("  MPIX_CUDA_AWARE_SUPPORT: not defined by this MPI implementation\n");
#endif

  // Runtime query (OpenMPI)
#if defined(MPIX_Query_cuda_support)
  // function form
  int supported = MPIX_Query_cuda_support();
  if (rank == 0)
    printf("  MPIX_Query_cuda_support(): %s\n", supported ? "YES" : "NO");
#else
  // Try weak symbol — may resolve at link time with some MPIs
  if (rank == 0)
    printf("  MPIX_Query_cuda_support(): not available\n");
#endif
}

// ---------------------------------------------------------------------------
int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);

  if (nranks != 2) {
    if (rank == 0)
      fprintf(stderr, "Error: this tool requires exactly 2 MPI ranks.\n"
                       "  mpirun -np 2 %s\n", argv[0]);
    MPI_Finalize();
    return 1;
  }

  // Print GPU info
  int dev;
  cudaDeviceProp prop;
  CUDA_CHECK(cudaGetDevice(&dev));
  CUDA_CHECK(cudaGetDeviceProperties(&prop, dev));
  printf("[rank %d] GPU %d: %s (CC %d.%d)\n", rank, dev, prop.name,
         prop.major, prop.minor);
  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 0) {
    printf("\n===== GPU-Aware MPI Diagnostic =====\n\n");
    printf("--- Step 1: Compile-time hints ---\n");
  }
  print_compile_time_hints();

  if (rank == 0)
    printf("\n--- Step 2: Functional test (MPI_Send/Recv with device pointer) ---\n");

  MPI_Barrier(MPI_COMM_WORLD);

  bool gpu_aware = false;
  // We can't easily catch a segfault inside MPI, so we just try it.
  // If MPI is not GPU-aware, this will likely crash — which is itself
  // a clear diagnostic. To be safe, set CUDA_MANAGED_FORCE_DEVICE_ALLOC=1
  // or use cudaMallocManaged as a fallback test first.
  gpu_aware = test_gpu_aware_mpi();

  int all_ok = gpu_aware ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &all_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  gpu_aware = (all_ok == 1);

  if (rank == 0) {
    if (gpu_aware)
      printf("  Result: GPU-aware MPI is FUNCTIONAL (data verified correctly)\n");
    else
      printf("  Result: GPU-aware MPI FAILED (data corruption detected)\n");
  }

  if (!gpu_aware) {
    if (rank == 0) {
      printf("\n===== VERDICT =====\n");
      printf("  GPU-aware MPI: NO\n");
      printf("  GPUDirect RDMA: N/A\n");
      printf("\nGPU-aware MPI is not working. Possible causes:\n");
      printf("  - MPI not built with CUDA support (e.g. OpenMPI needs --with-cuda)\n");
      printf("  - UCX not built with CUDA (--with-cuda for UCX)\n");
      printf("  - Missing: export UCX_MEMTYPE_CACHE=n (sometimes needed)\n");
      printf("  - MVAPICH2: ensure MV2_USE_CUDA=1 is set\n");
    }
    MPI_Finalize();
    return 1;
  }

  // -----------------------------------------------------------------------
  // Step 3: Benchmark to detect GPUDirect RDMA
  // -----------------------------------------------------------------------
  if (rank == 0)
    printf("\n--- Step 3: Performance benchmark (device-ptr vs staged host copy) ---\n");

  MPI_Barrier(MPI_COMM_WORLD);

  struct BenchResult {
    size_t bytes;
    double gpu_direct_us;
    double staged_us;
    double ratio;
  };

  // Test multiple sizes to get a clearer picture
  std::vector<size_t> sizes = {
      4 * 1024,         // 4 KB
      64 * 1024,        // 64 KB
      1024 * 1024,      // 1 MB
      16 * 1024 * 1024, // 16 MB
  };
  constexpr int NITERS = 50;

  std::vector<BenchResult> results;

  for (size_t sz : sizes) {
    double t_direct = ping_pong(sz, NITERS, /*use_device_ptr=*/true);
    double t_staged = ping_pong(sz, NITERS, /*use_device_ptr=*/false);
    results.push_back({sz, t_direct, t_staged, t_direct / t_staged});
  }

  if (rank == 0) {
    printf("\n  %12s  %14s  %14s  %10s\n", "Size", "GPU-ptr (us)", "Staged (us)",
           "Ratio");
    printf("  %12s  %14s  %14s  %10s\n", "----", "------------", "----------",
           "-----");
    for (auto &r : results) {
      const char *unit = "B";
      double display = r.bytes;
      if (r.bytes >= 1024 * 1024) {
        display = r.bytes / (1024.0 * 1024.0);
        unit = "MB";
      } else if (r.bytes >= 1024) {
        display = r.bytes / 1024.0;
        unit = "KB";
      }
      printf("  %8.0f %-3s  %12.1f    %12.1f    %8.2fx\n", display, unit,
             r.gpu_direct_us, r.staged_us, r.ratio);
    }

    // Heuristic: if direct is consistently faster (or comparable) to staged,
    // RDMA is likely active. If direct is slower, MPI is probably doing an
    // internal stage through host memory (extra copy overhead).
    double avg_ratio = 0;
    for (auto &r : results) avg_ratio += r.ratio;
    avg_ratio /= results.size();

    printf("\n===== VERDICT =====\n");
    printf("  GPU-aware MPI: YES\n");

    if (avg_ratio < 0.85) {
      printf("  GPUDirect RDMA: LIKELY YES (device-ptr path is %.0f%% faster "
             "than manual staging)\n",
             (1.0 - avg_ratio) * 100);
      printf("  The MPI implementation appears to use a fast zero-copy path.\n");
    } else if (avg_ratio < 1.15) {
      printf("  GPUDirect RDMA: UNCERTAIN (device-ptr and staged times are "
             "similar, ratio=%.2f)\n",
             avg_ratio);
      printf("  MPI may be GPU-aware but staging internally (no RDMA),\n");
      printf("  or RDMA is active but the extra D2H copy in staging is hidden "
             "by MPI overhead.\n");
      printf("  Check: nvidia-smi nvlink -s, or ibv_devinfo for GDR support.\n");
    } else {
      printf("  GPUDirect RDMA: LIKELY NO (device-ptr path is %.0f%% slower "
             "than manual staging)\n",
             (avg_ratio - 1.0) * 100);
      printf("  MPI accepts device pointers but likely stages through host "
             "internally.\n");
      printf("  To enable GPUDirect RDMA:\n");
      printf("    - Load nvidia_peermem kernel module: modprobe nvidia_peermem\n");
      printf("    - UCX: export UCX_TLS=rc,cuda_copy,cuda_ipc,gdr_copy\n");
      printf("    - UCX: export UCX_IB_GPU_DIRECT_RDMA=yes\n");
      printf("    - Verify: ucx_info -d | grep gdr\n");
      printf("    - MVAPICH2: export MV2_USE_GPUDIRECT=1\n");
    }
  }

  MPI_Finalize();
  return 0;
}
