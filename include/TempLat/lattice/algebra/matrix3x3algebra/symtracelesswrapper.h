#ifndef COSMOINTERFACE_MATRIX3X3ALGEBRA_SYMTRACELESSWRAPPER_H
#define COSMOINTERFACE_MATRIX3X3ALGEBRA_SYMTRACELESSWRAPPER_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026

#include "TempLat/lattice/algebra/matrix3x3algebra/allmatrixoperator.h"
#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/helpers/getstring.h"
#include "TempLat/util/rangeiteration/tagliteral.h"

#include "TempLat/lattice/algebra/helpers/preget.h"
#include "TempLat/lattice/algebra/helpers/postget.h"

namespace TempLat
{
  /** @brief A class which wraps two objects as a symmetric traceless .
   *
   *
   * Unit test: ctest -R test-symtracelesswrapperwrapper
   **/
  template <class R0, class R1, class R2, class R3, class R4, class R5>
  class SymTracelessWrapper : public SymTracelessOperator
  {
  public:
    // Put public methods here. These should change very little over time.

    SymTracelessWrapper(const R0 &pR0, const R1 &pR1, const R2 &pR2, const R3 &pR3, const R4 &pR4, const R5 &pR5)
        : mR0(pR0), mR1(pR1), mR2(pR2), mR3(pR3), mR4(pR4), mR5(pR5)
    {
    }

    auto SymTracelessGet(Tag<0> t) const { return (2. / 3.) * mR0 - (1. / 3.) * mR3 - (1. / 3.) * mR5; }
    auto SymTracelessGet(Tag<1> t) const { return mR1; }
    auto SymTracelessGet(Tag<2> t) const { return mR2; }
    auto SymTracelessGet(Tag<3> t) const { return -(1. / 3.) * mR0 + (2. / 3.) * mR3 - (1. / 3.) * mR5; }
    auto SymTracelessGet(Tag<4> t) const { return mR4; }
    auto SymTracelessGet(Tag<5> t) const { return -(1. / 3.) * mR0 - (1. / 3.) * mR3 + (2. / 3.) * mR5; }

    auto SymTracelessGet(Tag<1> t1, Tag<1> t2) const { return SymTracelessGet(0_c); }
    auto SymTracelessGet(Tag<1> t1, Tag<2> t2) const { return SymTracelessGet(1_c); }
    auto SymTracelessGet(Tag<1> t1, Tag<3> t2) const { return SymTracelessGet(2_c); }
    auto SymTracelessGet(Tag<2> t1, Tag<1> t2) const { return SymTracelessGet(1_c); }
    auto SymTracelessGet(Tag<2> t1, Tag<2> t2) const { return SymTracelessGet(3_c); }
    auto SymTracelessGet(Tag<2> t1, Tag<3> t2) const { return SymTracelessGet(4_c); }
    auto SymTracelessGet(Tag<3> t1, Tag<1> t2) const { return SymTracelessGet(2_c); }
    auto SymTracelessGet(Tag<3> t1, Tag<2> t2) const { return SymTracelessGet(4_c); }
    auto SymTracelessGet(Tag<3> t1, Tag<3> t2) const { return SymTracelessGet(5_c); }

    template <int N> auto operator()(Tag<N> t) const
    {
      static_assert(N >= 0 && N <= 4, "Operator(): N must be 0 or 4 for SymTracelessWrapper");
      return SymTracelessGet(t);
    }

    template <typename... IDX>
      requires requires(std::decay_t<R0> r0, std::decay_t<R1> r1, std::decay_t<R2> r2, std::decay_t<R3> r3,
                        std::decay_t<R4> r4, std::decay_t<R5> r5, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(r0, idx...);
        DoEval::eval(r1, idx...);
        DoEval::eval(r2, idx...);
        DoEval::eval(r3, idx...);
        DoEval::eval(r4, idx...);
        DoEval::eval(r5, idx...);
      }
    DEVICE_INLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      device::array<decltype(DoEval::eval(mR0, idx...)), 5> result;
      result[0] = DoEval::eval(mR0, idx...);
      result[1] = DoEval::eval(mR1, idx...);
      result[2] = DoEval::eval(mR2, idx...);
      result[3] = DoEval::eval(mR3, idx...);
      result[4] = DoEval::eval(mR4, idx...);
      const auto trace = result[0] + result[3] + DoEval::eval(mR5, idx...);
      result[0] -= 1. / 3. * trace;
      result[3] -= 1. / 3. * trace;
      return result;
    }

    void preGet()
    {
      PreGet::apply(mR0);
      PreGet::apply(mR1);
      PreGet::apply(mR2);
      PreGet::apply(mR3);
      PreGet::apply(mR4);
      PreGet::apply(mR5);
    }

    void postGet()
    {
      PostGet::apply(mR0);
      PostGet::apply(mR1);
      PostGet::apply(mR2);
      PostGet::apply(mR3);
      PostGet::apply(mR4);
      PostGet::apply(mR5);
    }

    std::string toString() const
    {
      return "SymTraceless(" + GetString::get(mR0) + "," + GetString::get(mR1) + "," + GetString::get(mR2) + "," +
             GetString::get(mR3) + "," + GetString::get(mR4) + ")";
    }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
    R0 mR0;
    R1 mR1;
    R2 mR2;
    R3 mR3;
    R4 mR4;
    R5 mR5;
  };

  template <typename R0, typename R1, typename R2, typename R3, typename R4, typename R5>
  SymTracelessWrapper<R0, R1, R2, R3, R4, R5> ConstructSymTraceless(const R0 &r0, const R1 &r1, const R2 &r2,
                                                                    const R3 &r3, const R4 &r4, const R5 &r5)
  {
    return {r0, r1, r2, r3, r4, r5};
  }
} // namespace TempLat

#endif
