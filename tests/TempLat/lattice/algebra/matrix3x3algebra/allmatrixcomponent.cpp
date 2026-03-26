/* This file is part of CosmoLattice, available at www.cosmolattice.net .
 C opyright Daniel G. Figueroa,* Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
 Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026
#include "TempLat/lattice/algebra/matrix3x3algebra/allmatrixcomponents.h"
#include "TempLat/util/rangeiteration/tagliteral.h"
#include "TempLat/util/rangeiteration/tag.h"

#include "TempLat/util/ndloop.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/symtracelessfield.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/algebra/operators/operators.h"

#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

    struct AllMatrixComponentsTester {
        static void Test(TDDAssertion &tdd);
    };

    void AllMatrixComponentsTester::Test(TDDAssertion &tdd)
    {
        struct MatrixStruct {
            DEVICE_FORCEINLINE_FUNCTION
            double MatrixGet(Tag<0> t) { return 6; };
            DEVICE_FORCEINLINE_FUNCTION
            int MatrixGet(Tag<1> t1, Tag<2> t2) { return 64; };
        };

        MatrixStruct matrix;
        tdd.verify(getComponent(matrix, 0_c) == 6);
        tdd.verify(getComponent(matrix, 1_c, 2_c) == 64);

        struct SymStruct {
            DEVICE_FORCEINLINE_FUNCTION
            double SymGet(Tag<0> t) { return 6; };
            DEVICE_FORCEINLINE_FUNCTION
            int SymGet(Tag<1> t1, Tag<1> t2) { return 64; };
        };

        SymStruct sym;
        tdd.verify(getComponent(sym, 0_c) == 6);
        tdd.verify(getComponent(sym, 1_c, 1_c) == 64);


        // Test whether it works with complex fields.
        constexpr size_t NDim = 3;
        using T = double;
        ptrdiff_t nGrid = 16, nGhost = 2;
        auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);
        toolBox->setVerbose();

        SpatialCoordinate<NDim> x(toolBox);

        SymTracelessField<T, NDim> fs("fs", toolBox);

        Field<T, NDim> x1("x1", toolBox);
        Field<T, NDim> x2("x2", toolBox);
        Field<T, NDim> x3("x3", toolBox);

        x1 = x(1_c);
        x2 = x(2_c);
        x3 = x(3_c);

        fs = ConstructSymTraceless(x1*x1, x2*x1, x3*x1, x2*x2, x3*x2, x3*x3);

        Field<T, NDim> f11("f11", toolBox);
        Field<T, NDim> f12("f12", toolBox);
        Field<T, NDim> f13("f13", toolBox);
        Field<T, NDim> f22("f22", toolBox);
        Field<T, NDim> f23("f23", toolBox);

        f11 = fs(1_c, 1_c);
        f12 = fs(1_c, 2_c);
        f13 = fs(1_c, 3_c);
        f22 = fs(2_c, 2_c);
        f23 = fs(2_c, 3_c);

        {
            auto view1 = x1.getLocalNDHostView();
            auto view2 = x2.getLocalNDHostView();
            auto view3 = x3.getLocalNDHostView();
            auto view11 = f11.getLocalNDHostView();
            auto view12 = f12.getLocalNDHostView();
            auto view13 = f13.getLocalNDHostView();
            auto view22 = f22.getLocalNDHostView();
            auto view23 = f23.getLocalNDHostView();
            // auto view33 = f33.getLocalNDHostView();
            bool all_correct = true;
            NDLoop<NDim>(view11, [&](const auto... idx) {
                bool this_correct = true;
                this_correct &= AlmostEqual(view11(idx...), 2. / 3. * view1(idx...)*view1(idx...) - 1. / 3. * view2(idx...)*view2(idx...) - 1. / 3. * view3(idx...)*view3(idx...) );
                this_correct &= AlmostEqual(view12(idx...), view1(idx...)*view2(idx...));
                this_correct &= AlmostEqual(view13(idx...), view1(idx...)*view3(idx...));
                this_correct &= AlmostEqual(view22(idx...), - 1. / 3. * view1(idx...)*view1(idx...) + 2. / 3. * view2(idx...)*view2(idx...) - 1. / 3. * view3(idx...)*view3(idx...));
                this_correct &= AlmostEqual(view23(idx...), view2(idx...)*view3(idx...));
                // this_correct &= AlmostEqual(view33(idx...), - 1. / 3. * view1(idx...)*view1(idx...) - 1. / 3. * view2(idx...)*view2(idx...) + 2. / 3. * view3(idx...)*view3(idx...));
                if (!this_correct) {
                    std::cout << "Matrix operation test failed at index (";
                    ((std::cout << idx << ", "), ...);
                    std::cout << ") got (" << view11(idx...) << ", "
                                           << view12(idx...) << ", "
                                           << view13(idx...) << ", "
                                           << view22(idx...) << ", "
                                           << view23(idx...) << "), expected ("
                                           << 2. / 3. * view1(idx...)*view1(idx...) - 1. / 3. * view2(idx...)*view2(idx...) - 1. / 3. * view3(idx...)*view3(idx...) << ", "
                                           << view1(idx...)*view2(idx...) << ", "
                                           << view1(idx...)*view3(idx...) << ", "
                                           << - 1. / 3. * view1(idx...)*view1(idx...) + 2. / 3. * view2(idx...)*view2(idx...) - 1. / 3. * view3(idx...)*view3(idx...) << ", "
                                           << view2(idx...)*view3(idx...) << ")\n";
                }
                all_correct = all_correct && this_correct;
            });
            tdd.verify(all_correct);
        }

    }

} // namespace TempLat

namespace
{
    TempLat::TDDContainer<TempLat::AllMatrixComponentsTester> test;
}
