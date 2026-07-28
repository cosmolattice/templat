# ##############################################################################
# Get Kokkos━FFT
# ##############################################################################

message("")
message(
  "${BoldYellow}┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Fetching KokkosFFT ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓${ColorReset}"
)

FetchContent_Declare(
  KokkosFFT
  DOWNLOAD_EXTRACT_TIMESTAMP FALSE
  URL https://github.com/kokkos/kokkos-fft/archive/refs/tags/v1.1.0.zip
  URL_HASH
    SHA256=01000aca33231ea18006a822a321bd737545b2156ad2caffa3b380b3e13f1c30
  ${SYSTEM_MARKER})
FetchContent_MakeAvailable(KokkosFFT)

message(
  "${BoldYellow}┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Fetching KokkosFFT DONE ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛${ColorReset}\n"
)
