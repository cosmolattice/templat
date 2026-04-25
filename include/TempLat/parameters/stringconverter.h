#ifndef TEMPLAT_PARAMETERS_STRINGCONVERTER_H
#define TEMPLAT_PARAMETERS_STRINGCONVERTER_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2019

#include "TempLat/lattice/ghostcells/boundaryconditions.h"
#include "TempLat/parameters/multipleparametergetter.h"
#include "TempLat/util/exception.h"
#include <cctype>
#include <sstream>
#include <string>
#include <vector>

namespace TempLat
{
  MakeException(BCKeywordUnknown);
  MakeException(BCSpecArityMismatch);

  /** @brief A class which wraps a single function, splitting it in lines and passing each
   *  line to ParameterGetter, returning the result in your provided MultipleParameterGetter<T>.
   *
   * Unit test: ctest -R test-stringconverter
   **/
  template <class T> class StringConverter
  {
  public:
    // Put public methods here. These should change very little over time.
    StringConverter() = default;
    void operator()(const std::string &str, MultipleParameterGetter<T> &arr, const std::string &name)
    {
      arr.clear();
      T tmp;
      std::istringstream iss(str);

      while (iss >> std::boolalpha >> std::skipws >> tmp) {
        arr.push_back(ParameterGetter<T>(tmp, name));
      }
    }
  };

  inline BCType stringToBCType(const std::string &s)
  {
    std::string k;
    k.reserve(s.size());
    for (char c : s) {
      if (!std::isspace(static_cast<unsigned char>(c)))
        k.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }
    if (k == "p" || k == "pbc" || k == "periodic") return BCType::Periodic;
    if (k == "ap" || k == "apbc" || k == "antiperiodic") return BCType::Antiperiodic;
    if (k == "d" || k == "dbc" || k == "dirichlet" || k == "open") return BCType::Dirichlet;
    if (k == "n" || k == "nbc" || k == "neumann") return BCType::Neumann;
    throw(BCKeywordUnknown("Unknown boundary-condition keyword: '" + s +
                           "'. Expected one of: p/pbc/periodic, ap/apbc/antiperiodic, "
                           "d/dbc/dirichlet/open, n/nbc/neumann."));
  }

  inline std::string bcTypeToString(BCType bc)
  {
    switch (bc) {
    case BCType::Periodic: return "periodic";
    case BCType::Antiperiodic: return "antiperiodic";
    case BCType::Dirichlet: return "dirichlet";
    case BCType::Neumann: return "neumann";
    }
    return "unknown";
  }

  /** @brief Parse a comma-separated BC keyword list (e.g. "apbc,pbc,pbc") into a `BCSpec<NDim>`.
   *
   * The list must have exactly NDim entries. Whitespace around tokens is allowed.
   * Unknown keywords throw `BCKeywordUnknown`; wrong arity throws `BCSpecArityMismatch`.
   * Recognized keywords (case-insensitive): p/pbc/periodic, ap/apbc/antiperiodic,
   * d/dbc/dirichlet/open, n/nbc/neumann.
   */
  template <size_t NDim> BCSpec<NDim> parseBCSpec(const std::string &s)
  {
    std::vector<std::string> tokens;
    std::string item;
    std::stringstream ss(s);
    while (std::getline(ss, item, ',')) tokens.push_back(item);
    if (tokens.size() != NDim) {
      throw(BCSpecArityMismatch("BC spec '" + s + "' has " + std::to_string(tokens.size()) +
                                " entries; expected exactly " + std::to_string(NDim) + "."));
    }
    BCSpec<NDim> out{};
    for (size_t i = 0; i < NDim; ++i) out[i] = stringToBCType(tokens[i]);
    return out;
  }

  /** @brief Serialize a `BCSpec<NDim>` back to a canonical comma-separated keyword list. */
  template <size_t NDim> std::string serializeBCSpec(const BCSpec<NDim> &bc)
  {
    std::string out;
    for (size_t i = 0; i < NDim; ++i) {
      if (i) out.push_back(',');
      out += bcTypeToString(bc[i]);
    }
    return out;
  }
} // namespace TempLat

#endif
