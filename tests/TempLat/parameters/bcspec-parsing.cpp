
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

#include "TempLat/lattice/ghostcells/boundaryconditions.h"
#include "TempLat/parameters/stringconverter.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{
  struct BCSpecParsingTester {
    static void Test(TDDAssertion &tdd);
  };

  static bool throwsBCKeyword(const std::string &s)
  {
    try {
      stringToBCType(s);
    } catch (const BCKeywordUnknown &) {
      return true;
    }
    return false;
  }

  template <size_t NDim> static bool throwsArity(const std::string &s)
  {
    try {
      parseBCSpec<NDim>(s);
    } catch (const BCSpecArityMismatch &) {
      return true;
    }
    return false;
  }

  template <size_t NDim> static bool throwsKeyword(const std::string &s)
  {
    try {
      parseBCSpec<NDim>(s);
    } catch (const BCKeywordUnknown &) {
      return true;
    }
    return false;
  }

  void BCSpecParsingTester::Test(TDDAssertion &tdd)
  {
    // ----- stringToBCType: every valid keyword in every casing variant -----
    tdd.verify(stringToBCType("p") == BCType::Periodic);
    tdd.verify(stringToBCType("pbc") == BCType::Periodic);
    tdd.verify(stringToBCType("periodic") == BCType::Periodic);
    tdd.verify(stringToBCType("PBC") == BCType::Periodic);
    tdd.verify(stringToBCType("Periodic") == BCType::Periodic);
    tdd.verify(stringToBCType("PERIODIC") == BCType::Periodic);

    tdd.verify(stringToBCType("ap") == BCType::Antiperiodic);
    tdd.verify(stringToBCType("apbc") == BCType::Antiperiodic);
    tdd.verify(stringToBCType("antiperiodic") == BCType::Antiperiodic);
    tdd.verify(stringToBCType("APBC") == BCType::Antiperiodic);
    tdd.verify(stringToBCType("AntiPeriodic") == BCType::Antiperiodic);

    tdd.verify(stringToBCType("d") == BCType::Dirichlet);
    tdd.verify(stringToBCType("dbc") == BCType::Dirichlet);
    tdd.verify(stringToBCType("dirichlet") == BCType::Dirichlet);
    tdd.verify(stringToBCType("open") == BCType::Dirichlet);
    tdd.verify(stringToBCType("OPEN") == BCType::Dirichlet);
    tdd.verify(stringToBCType("Dirichlet") == BCType::Dirichlet);

    tdd.verify(stringToBCType("n") == BCType::Neumann);
    tdd.verify(stringToBCType("nbc") == BCType::Neumann);
    tdd.verify(stringToBCType("neumann") == BCType::Neumann);
    tdd.verify(stringToBCType("NEUMANN") == BCType::Neumann);

    // whitespace tolerance on a single keyword
    tdd.verify(stringToBCType("  pbc  ") == BCType::Periodic);

    // ----- stringToBCType: unknown keyword throws -----
    tdd.verify(throwsBCKeyword(""));
    tdd.verify(throwsBCKeyword("foo"));
    tdd.verify(throwsBCKeyword("free"));
    tdd.verify(throwsBCKeyword("twisted"));

    // ----- parseBCSpec: arity errors -----
    tdd.verify(throwsArity<3>("pbc,pbc"));     // too few
    tdd.verify(throwsArity<3>("pbc,pbc,pbc,pbc")); // too many
    tdd.verify(throwsArity<2>("pbc"));         // single token, expect 2
    tdd.verify(throwsArity<1>(""));            // empty, expect 1

    // ----- parseBCSpec: unknown-token error names what's wrong -----
    tdd.verify(throwsKeyword<3>("pbc,quasi,pbc"));

    // ----- parseBCSpec: every BC kind round-trips -----
    {
      auto bc = parseBCSpec<4>("pbc,apbc,dbc,nbc");
      tdd.verify(bc[0] == BCType::Periodic);
      tdd.verify(bc[1] == BCType::Antiperiodic);
      tdd.verify(bc[2] == BCType::Dirichlet);
      tdd.verify(bc[3] == BCType::Neumann);
    }

    // ----- parseBCSpec: whitespace and casing tolerated inside a list -----
    {
      auto bc = parseBCSpec<3>(" Periodic , APBC ,  open ");
      tdd.verify(bc[0] == BCType::Periodic);
      tdd.verify(bc[1] == BCType::Antiperiodic);
      tdd.verify(bc[2] == BCType::Dirichlet);
    }

    // ----- "all-periodic" is exactly equal to allPeriodic<3>() — no silent default -----
    {
      auto parsed = parseBCSpec<3>("pbc,pbc,pbc");
      auto canonical = allPeriodic<3>();
      tdd.verify(parsed == canonical);
    }

    // ----- bcTypeToString is inverse of stringToBCType on canonical names -----
    tdd.verify(bcTypeToString(BCType::Periodic) == "periodic");
    tdd.verify(bcTypeToString(BCType::Antiperiodic) == "antiperiodic");
    tdd.verify(bcTypeToString(BCType::Dirichlet) == "dirichlet");
    tdd.verify(bcTypeToString(BCType::Neumann) == "neumann");

    // ----- serialize → parse → serialize round-trip -----
    {
      BCSpec<3> original{BCType::Antiperiodic, BCType::Periodic, BCType::Neumann};
      std::string s1 = serializeBCSpec<3>(original);
      auto reparsed = parseBCSpec<3>(s1);
      tdd.verify(reparsed == original);
      std::string s2 = serializeBCSpec<3>(reparsed);
      tdd.verify(s1 == s2);
    }

    // ----- parse → serialize → parse round-trip across all four kinds -----
    {
      const std::string in = "apbc,pbc,dbc,nbc";
      auto bc = parseBCSpec<4>(in);
      std::string ser = serializeBCSpec<4>(bc);
      auto bc2 = parseBCSpec<4>(ser);
      tdd.verify(bc == bc2);
    }
  }
} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::BCSpecParsingTester> test;
}
