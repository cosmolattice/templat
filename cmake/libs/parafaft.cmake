message("") # Blank line for better readability
message(
  "${BoldYellow}┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Fetching ParaFaFT ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓${ColorReset}"
)

# enable FetchContent
set(PARAFAFT_CUDA
    ${CUDA}
    CACHE BOOL "")
set(PARAFAFT_HIP
    ${HIP}
    CACHE BOOL "")
set(PARAFAFT_FFTW_INCLUDE_DIR
    ${FFTW_INCLUDES}
    CACHE PATH "")
set(PARAFAFT_FFTW_LIB_DIR
    ${FFTW_LIB_DIR}
    CACHE PATH "")
# declare ParaFaFT
FetchContent_Declare(
  ParaFaFT # name of the content
  GIT_REPOSITORY https://github.com/cosmolattice/parafaft.git # the repository
  GIT_TAG main # the tag
)
# make available
FetchContent_MakeAvailable(ParaFaFT)

message(
  STATUS "ParaFaFT available with CUDA=${PARAFAFT_CUDA} HIP=${PARAFAFT_HIP}")

message(
  "${BoldYellow}┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ Fetching ParaFaFT DONE ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛${ColorReset}\n"
)
