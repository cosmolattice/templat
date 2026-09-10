
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019
#include "TempLat/parameters/stringconverter.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{
  MakeException(MockKeywordUnknown);

  /** A stand-in for the enum-like parameter types that live in CosmoInterface (EvolverType,
   * InitialConditionsType). They all follow one pattern: pull a whitespace-delimited token with
   * `in >> tmp`, map it, and throw on anything unrecognised. StringConverter's end-of-input check
   * has to stay correct for that pattern, and those headers are not reachable from here, so the
   * pattern is reproduced rather than imported. */
  struct MockKeyword {
    int value = -1;
  };

  inline std::istream &operator>>(std::istream &in, MockKeyword &k)
  {
    std::string tmp;
    in >> tmp;
    if (tmp == "alpha")
      k.value = 0;
    else if (tmp == "beta")
      k.value = 1;
    else if (tmp.empty()) {
    } // leaves failbit set: the loop must stop without pushing
    else
      throw(MockKeywordUnknown(tmp + " is not a MockKeyword."));
    return in;
  }

  struct StringConverterTester {
    static void Test(TDDAssertion &tdd);
  };

  void StringConverterTester::Test(TDDAssertion &tdd)
  {
    // Does `conv(v, out, "p")` throw the given exception type?
    auto throwsAs = [](auto &conv, auto &out, const std::string &v, auto tag) {
      using E = decltype(tag);
      try {
        conv(v, out, "p");
      } catch (const E &) {
        return true;
      } catch (...) {
        return false;
      }
      return false;
    };
    ParameterParserUnparseableValue bad("");
    MockKeywordUnknown unknown("");

    // ---- doubles: the original case, plus the numeric spellings that must keep working --------
    StringConverter<double> strd;
    MultipleParameterGetter<double> dbles;
    strd("54 67 89", dbles, "test");
    tdd.verify(dbles.size() == 3 && dbles[0]() == 54 && dbles[1]() == 67 && dbles[2]() == 89);
    strd("1e-3", dbles, "test");
    tdd.verify(dbles.size() == 1 && dbles[0]() == 1e-3);
    strd("-2.5 +3 .5", dbles, "test");
    tdd.verify(dbles.size() == 3 && dbles[0]() == -2.5 && dbles[1]() == 3 && dbles[2]() == .5);
    strd("  7.5  ", dbles, "test"); // surrounding whitespace is not trailing junk
    tdd.verify(dbles.size() == 1 && dbles[0]() == 7.5);

    // ---- doubles: malformed input must be refused, not silently dropped ----------------------
    tdd.verify(throwsAs(strd, dbles, "not_a_number", bad));
    tdd.verify(throwsAs(strd, dbles, "0.005 leftover", bad));
    tdd.verify(throwsAs(strd, dbles, "1.0 2.0x", bad));
    tdd.verify(throwsAs(strd, dbles, "abc 1.0", bad));

    // ---- bools: only the boolalpha literals. 0/1 parse into nothing at all, which used to let
    // the caller's default stand while the parameter file still showed what the user wrote ------
    StringConverter<bool> strb;
    MultipleParameterGetter<bool> bools;
    strb("true", bools, "flag");
    tdd.verify(bools.size() == 1 && bools[0]() == true);
    strb("false", bools, "flag");
    tdd.verify(bools.size() == 1 && bools[0]() == false);
    strb("true false true", bools, "flag");
    tdd.verify(bools.size() == 3 && bools[0]() == true && bools[1]() == false && bools[2]() == true);
    for (const auto &v : {"0", "1", "2", "yes", "no", "True", "TRUE", "FALSE", "tru", "truex"})
      tdd.verify(throwsAs(strb, bools, v, bad));

    // ---- integers: "3.5" for an int is the same class of silent mistake ----------------------
    StringConverter<int> stri;
    MultipleParameterGetter<int> ints;
    stri("42", ints, "n");
    tdd.verify(ints.size() == 1 && ints[0]() == 42);
    stri("-7 8", ints, "n");
    tdd.verify(ints.size() == 2 && ints[0]() == -7 && ints[1]() == 8);
    tdd.verify(throwsAs(stri, ints, "abc", bad));
    tdd.verify(throwsAs(stri, ints, "3.5", bad));

    // ---- strings: a path is one token; several words stay several values ---------------------
    StringConverter<std::string> strs;
    MultipleParameterGetter<std::string> strings;
    strs("./diag2/d96/", strings, "outputfile");
    tdd.verify(strings.size() == 1 && strings[0]() == "./diag2/d96/");
    strs("a b c", strings, "words");
    tdd.verify(strings.size() == 3);

    // ---- blank input is the caller's business (an absent value already throws upstream) ------
    tdd.verify(!throwsAs(strd, dbles, "   ", bad));
    tdd.verify(!throwsAs(strd, dbles, "", bad));
    strd("   ", dbles, "test");
    tdd.verify(dbles.size() == 0);

    // ---- enum-like types keep working, and keep their own exception ---------------------------
    StringConverter<MockKeyword> strk;
    MultipleParameterGetter<MockKeyword> keys;
    strk("alpha", keys, "kw");
    tdd.verify(keys.size() == 1 && keys[0]().value == 0);
    strk("alpha beta", keys, "kw");
    tdd.verify(keys.size() == 2 && keys[0]().value == 0 && keys[1]().value == 1);
    tdd.verify(throwsAs(strk, keys, "gamma", unknown)); // not converted into an Unparseable
    tdd.verify(!throwsAs(strk, keys, "  ", bad));

    // ---- a throwing call must not leave stale values behind for the next one ------------------
    try {
      strd("1.0 2.0 junk", dbles, "test");
    } catch (const ParameterParserUnparseableValue &) {
    }
    strd("9", dbles, "test");
    tdd.verify(dbles.size() == 1 && dbles[0]() == 9);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::StringConverterTester> test;
}
