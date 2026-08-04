if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColorReset "${Esc}[m")
  set(ColorBold "${Esc}[1m")
  set(Red "${Esc}[31m")
  set(Green "${Esc}[32m")
  set(Yellow "${Esc}[33m")
  set(Blue "${Esc}[34m")
  set(Magenta "${Esc}[35m")
  set(Cyan "${Esc}[36m")
  set(White "${Esc}[37m")
  set(BoldRed "${Esc}[1;31m")
  set(BoldGreen "${Esc}[1;32m")
  set(BoldYellow "${Esc}[1;33m")
  set(BoldBlue "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan "${Esc}[1;36m")
  set(BoldWhite "${Esc}[1;37m")
endif()

include(FetchContent)
set(FETCHCONTENT_QUIET
    ON
    CACHE BOOL "Suppress output from FetchContent during configuration")
mark_as_advanced(FETCHCONTENT_QUIET)
# FetchContent_Declare() gained the SYSTEM keyword in CMake 3.25, but 3.25.0 and
# 3.25.1 forwarded it to ExternalProject_Add(), which does not know it and so
# appends it to the preceding keyword's value list (e.g. the declarations below
# end up with URL_HASH "SHA256=...;SYSTEM" and the configure fails). Fixed
# upstream in 3.25.2, see CMake issue #24201. Debian 12 ships 3.25.1, so this
# guard has to be 3.25.2 rather than 3.25.
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.25.2")
  set(SYSTEM_MARKER
      SYSTEM
      CACHE INTERNAL "")
else()
  set(SYSTEM_MARKER
      ""
      CACHE INTERNAL "")
endif()
