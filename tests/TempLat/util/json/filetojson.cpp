
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
#include "TempLat/util/json/filetojson.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/namedtmpfile.h"
#include "TempLat/util/log/log.h"

namespace TempLat
{

  struct FileToJSONTester {
    static void Test(TDDAssertion &tdd);
  };

  void FileToJSONTester::Test(TDDAssertion &tdd)
  {

    NamedTmpFile ntf;

    ntf << "{\"physics\" : {\"q\" : 1e5, \"g2\" : 4e-5, \"startNorm\" : 1e-10, \"fStar\" : 0.1, \"phi\" : { "
           "\"randomSeed\" : \"Helllooooo!!!\"}, \"chi\" : { \"randomSeed\" : \"Hello to you too. Entropy "
           "perturbations, "
           "aye?\"} }, \"integrator\" : { \"tEnd\" : 200, \"dt\" : 0.0002, \"outputFrequency\" : 1, "
           "\"conformalPowerAlpha\" : 3, \"scaleFactor_initial\" : 1 }, \"output\" : { \"path\" : \"../demoOutput\" } "
           "}";

    std::string fName = ntf.getName();

    ntf.close();

    FileToJSON ftj(fName);

    // say << ftj;

    ntf.remove();

    tdd.verify(ftj.size() != 0);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::FileToJSONTester> test;
}
