#ifndef TEMPLAT_LATTICE_MEASUREMENTS_PROJECTIONHELPERS_UNBINNEDRADIALPROJECTIONRESULT_H
#define TEMPLAT_LATTICE_MEASUREMENTS_PROJECTIONHELPERS_UNBINNEDRADIALPROJECTIONRESULT_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026
//            Based on: Wessel Valkenburg, Year: 2019

#include "TempLat/util/exception.h"
#include "TempLat/lattice/algebra/helpers/getfloattype.h"
#include "TempLat/lattice/measuringtools/projectionhelpers/radialprojectionsinglebinandvalue.h"
#include "TempLat/lattice/measuringtools/projectionhelpers/radialprojectionsinglequantity.h"
#include "TempLat/lattice/measuringtools/projectionhelpers/radialprojectionrebinner.h"

#include "TempLat/parallel/device.h"

namespace TempLat
{

  /** @brief A class which holds the result of a radial projection without binning.
   *  For each possible value of the Fourier coordinate it saves:
   *  the average value,
   *  the sample variance (= <f^2> - <f>^2 ) of the values,
   *  and the number of entries in a bin (multiplicity).
   *
   *  Upon call to finalize, constructs the parent vector, which you can use as you wish.
   *  This contains pairs for each value of the momentum, with the momentum (first) and
   *  the value result (second), as a `RadialProjectionSingleDatum'
   *
   *  Unit test: ctest -R test-unbinnedradialprojectionresult
   **/

  template <typename T> class UnbinnedRadialProjectionResult : public std::vector<std::tuple<T, RadialProjectionSingleDatum<T>>>
  {
  public:
    using floatType = typename GetFloatType<T>::type;

    // Put public methods here. These should change very little over time.
    UnbinnedRadialProjectionResult(size_t nBins, bool pIsInFourier = false) :
    std::vector<std::tuple<T, RadialProjectionSingleDatum<T>>>(),
    finalizedOnce(false),
    mNBins(nBins),
    mValues(mNBins),
    unset(false),
    mIsInFourier(pIsInFourier)
    {
      mMultiplicitiesDevice = DeviceView("RadialProjectionResult::mMultiplicitiesDevice", mNBins);
      mMultiplicities = device::memory::createMirrorView(mMultiplicitiesDevice);
    }

    auto getNBins()
    {
      return mNBins;
    }

    /** \brief Rescale the results with a function of x or k (bin location),
     *  using for now a simple lambda function of single float parameter, which
     *  is your x in f(x). In other words, you give f(x).
     */
    template <typename LL>
    UnbinnedRadialProjectionResult& rescale(LL rescaler) {
      //for (auto&& it : *this) {
      for (size_t i = 0; i< this->size(); ++i)
      {
        std::get<1>((*this)[i]) *= rescaler( std::get<0>((*this)[i]) );
      }
      return *this;
    }

    /** @brief Rescale the bin positions with a normalization (for example dimensionful).
     */
    UnbinnedRadialProjectionResult& rescaleBins(T scale) {
      for (auto&& it : *this) {
        std::get<0>(it) *= scale;
      }
      return *this;
    }

    auto integrate(bool dummy)
    {
      auto total = 0.;
      for (size_t i = 0; i < this->size(); ++i) {
        total +=  std::get<1>((*this)[i]).average / std::get<0>((*this)[i]) * std::get<1>((*this)[i]).multiplicity;
      }
      return total;
    }

    std::string toString(int verbosity = 0) const
    {
      if ((ptrdiff_t)this->size() < 1) return "";
      std::stringstream sstream;
      sstream << this->front().getHeader(verbosity) << "\n";
      for (auto &&it : *this) {
        sstream << it.toString(verbosity) << "\n";
      }
      return sstream.str();
    }

    friend std::ostream &operator<<(std::ostream &ostream, const UnbinnedRadialProjectionResult &rpr)
    {
      ostream << rpr.toString() << "\n";
      return ostream;
    }

    template <typename S> friend class RadialProjector;

    template <typename S>
    friend UnbinnedRadialProjectionResult<S> operator+(const UnbinnedRadialProjectionResult<S> &a, const UnbinnedRadialProjectionResult<S> &b);

    DEVICE_FORCEINLINE_FUNCTION
    void add_device(ptrdiff_t i, const T &value, const T &weight = (T)1) const
    {
      mValues.add_device(i, value, weight);
      device::atomic_add(&mMultiplicitiesDevice(i), weight);
    }

    /** @brief RadialProjector calls this as the last step, does the transposition of the result vecotrs into one vector
     * of RadialProjectionSingleBinAndValue<T>. */
    UnbinnedRadialProjectionResult &finalize(MPICommReference comm)
    {
      if (finalizedOnce) throw RadialProjectionResultFinalizationException("Can only finalize once per instance.");

      finalizedOnce = true;
      pull();

      comm.Allreduce(mMultiplicities, MPI_SUM);

      mValues.finalize(comm);

      auto &mFullResult = *this;

      mFullResult.clear();
      for (size_t i = 0; i < mNBins; ++i) {
        if (mMultiplicities[i] > 0.) {
          std::tuple<T, RadialProjectionSingleDatum<T>> next(std::sqrt(i), mValues.getFinal(i, mMultiplicities[i]));
          this->push_back(next);
        }
      }

      return *this;
    }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
    bool finalizedOnce;
    size_t mNBins;
    RadialProjectionSingleQuantity<T> mValues;

    using DeviceView = device::memory::NDView<floatType, 1>;
    using HostMirror = typename DeviceView::host_mirror_type;

    HostMirror mMultiplicities;
    DeviceView mMultiplicitiesDevice;

    bool unset;
    bool mUseBinCentralValues;
    bool mIsInFourier;

    void pull()
    {
      mValues.pull();
      device::memory::copyDeviceToHost(mMultiplicitiesDevice, mMultiplicities.data());
    }

    void push()
    {
      mValues.push();
      device::memory::copyHostToDevice(mMultiplicities.data(), mMultiplicitiesDevice);
    }
  };

  template <typename T>
  UnbinnedRadialProjectionResult<T> operator+(const UnbinnedRadialProjectionResult<T> &a, const UnbinnedRadialProjectionResult<T> &b)
  {
    UnbinnedRadialProjectionResult<T> res(a);

    for (size_t i = 0; i < a.size(); ++i) {
      std::get<1>(res[i]).average += std::get<1>(b[i]).average;
    }
    return res;
  }

  template <typename T, class R> UnbinnedRadialProjectionResult<T> operator*(R &&func, const UnbinnedRadialProjectionResult<T> &obj)
  {
    UnbinnedRadialProjectionResult<T> res(obj);
    return res.rescale(func);
  }

  template <typename T> UnbinnedRadialProjectionResult<T> operator*(double scale, const UnbinnedRadialProjectionResult<T> &obj)
  {
    auto func = [&](auto x) { return scale; };
    return func * obj;
  }

  template <typename T> UnbinnedRadialProjectionResult<T> operator*(float scale, const UnbinnedRadialProjectionResult<T> &obj)
  {
    auto func = [&](auto x) { return scale; };
    return func * obj;
  }

  template <typename T> UnbinnedRadialProjectionResult<T> operator*(int scale, const UnbinnedRadialProjectionResult<T> &obj)
  {
    auto func = [&](auto x) { return scale; };
    return func * obj;
  }
} // namespace TempLat

#endif
