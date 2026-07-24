
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

/** \file Tests the umbrella header TempLat.h.
 *
 * Compiling this translation unit at all is the first half of the test: it proves that
 * <TempLat.h> pulls in every public header of the library without collisions, in whichever
 * configuration TempLat was built.
 *
 * The second half runs at test time and guards against the umbrella going stale: it walks
 * include/TempLat, applies the same exclusion policy that TempLat.h documents, and verifies
 * that every remaining header is listed. A newly added header that nobody remembered to add
 * to TempLat.h fails here rather than silently going missing from the public surface.
 */

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "TempLat.h"

#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct TempLatUmbrellaTester {
    static void Test(TDDAssertion &tdd);
  };

  namespace
  {
    /** @brief The header groups that TempLat.h deliberately leaves out; see its file comment. */
    bool isExcludedFromUmbrella(const std::string &relativePath)
    {
      static const std::vector<std::string> excludedPrefixes = {
          "TempLat/fft/external/",      // backend internals, selected by fft/fftlibraryselector.h
          "TempLat/parallel/devices/",  // backend internals, selected by parallel/device.h
          "TempLat/util/tdd/",          // the unit-test harness itself
          "TempLat/util/hash/libkeccak" // private, refuses to be included on its own
      };
      for (const auto &prefix : excludedPrefixes)
        if (relativePath.rfind(prefix, 0) == 0) return true;

      // Empty placeholder headers that only exist to give a test a name.
      return relativePath.ends_with("tester.h") || relativePath.ends_with("log_test.h");
    }

    /** @brief The repository's include directory, derived from this file's own location. */
    std::filesystem::path includeDirectory()
    {
      // <repo>/tests/TempLat/templat.cpp -> <repo>/include
      return std::filesystem::path(__FILE__).parent_path().parent_path().parent_path() / "include";
    }
  } // namespace

  // A handful of spot checks that the umbrella really did define the library, spread over the
  // modules it claims to cover. Merely naming these types is the test.
  static_assert(std::is_class_v<SessionGuard>);
  static_assert(std::is_class_v<MPICommReference>);
  static_assert(std::is_class_v<MemoryToolBox<3>>);
  static_assert(std::is_class_v<Field<double, 3>>);
  static_assert(std::is_class_v<Averager<Field<double, 3>>>);
  static_assert(std::is_class_v<FFTLibrarySelector<3>>);
  static_assert(std::is_class_v<ParameterParser>);

  void TempLatUmbrellaTester::Test(TDDAssertion &tdd)
  {
    const std::filesystem::path includeDir = includeDirectory();
    const std::filesystem::path umbrella = includeDir / "TempLat.h";

    // Out-of-tree copies of the tests (installed, packaged) cannot run the coverage half.
    if (!std::filesystem::exists(umbrella)) {
      sayShort << "Skipping the umbrella coverage check: " << umbrella << " is not reachable from the test binary.\n";
      return;
    }

    std::string umbrellaText;
    {
      std::ifstream in(umbrella);
      umbrellaText.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    }

    std::vector<std::string> missing;
    for (const auto &entry : std::filesystem::recursive_directory_iterator(includeDir / "TempLat")) {
      if (!entry.is_regular_file() || entry.path().extension() != ".h") continue;

      const std::string relativePath = std::filesystem::relative(entry.path(), includeDir).generic_string();
      if (isExcludedFromUmbrella(relativePath)) continue;

      if (umbrellaText.find("\"" + relativePath + "\"") == std::string::npos) missing.push_back(relativePath);
    }

    std::sort(missing.begin(), missing.end());
    for (const auto &header : missing)
      sayShort << "TempLat.h does not include " << header << "\n";

    tdd.verify(missing.empty(), "TempLat.h must include every public header of the library.");
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::TempLatUmbrellaTester> test;
}
