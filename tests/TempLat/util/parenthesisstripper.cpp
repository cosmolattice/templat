
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2020
#include "TempLat/util/parenthesisstripper.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct ParenthesisStripperTester {
    static void Test(TDDAssertion &tdd);
  };

  void ParenthesisStripperTester::Test(TDDAssertion &tdd)
  {
    std::string tmp("(a*2)");
    ParenthesisStripper::strip(tmp);
    tdd.verify(tmp == "a*2");
    tmp = "[a*2]";
    ParenthesisStripper::strip(tmp);
    tdd.verify(tmp == "a*2");
    tmp = "{a*2}";
    ParenthesisStripper::strip(tmp);
    tdd.verify(tmp == "a*2");
    tmp = "{a*2}+3";
    ParenthesisStripper::strip(tmp);
    tdd.verify(tmp == "{a*2}+3");
    tmp = "({a*2}+3)";
    ParenthesisStripper::strip(tmp);
    tdd.verify(tmp == "{a*2}+3");

    tmp = "(3*(2+2/2)-2*8-5)*9/(2+4*4)";
    ParenthesisStripper::strip(tmp);
    tdd.verify(tmp == "(3*(2+2/2)-2*8-5)*9/(2+4*4)");
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::ParenthesisStripperTester> test;
}
