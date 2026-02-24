# ##############################################################################
# Get Kokkos-FFT
# ##############################################################################

message(
  STATUS
    "${BoldYellow}--------------------- Fetching KokkosFFT ---------------------${ColorReset}"
)

include(FetchContent)
FetchContent_Declare(
  KokkosFFT
  DOWNLOAD_EXTRACT_TIMESTAMP FALSE
  URL https://github.com/kokkos/kokkos-fft/archive/refs/tags/v1.0.0.zip
  URL_HASH
    SHA256=80e9c1abdf71df2342ae713c845ba2aeabf1c1a0c2e116795534a8e67734c177
  SYSTEM)
FetchContent_MakeAvailable(KokkosFFT)

message(
  STATUS
    "${BoldYellow}------------------ Fetching KokkosFFT DONE -------------------${ColorReset}\n"
)
