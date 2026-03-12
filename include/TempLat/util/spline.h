/*
 * spline.h
 *
 * GPU-compatible cubic spline interpolation library.
 * Based on the cubic spline library by Tino Kluge (ttk448 at gmail.com).
 *
 * Two-layer design:
 *   - SplineData: lightweight GPU-evaluatable handle (raw pointers, no STL)
 *   - Spline: host-side owner that fits coefficients and manages device memory
 *
 * ---------------------------------------------------------------------
 * Copyright (C) 2011, 2014, 2016, 2021 Tino Kluge (ttk448 at gmail.com)
 * Modified for GPU compatibility, 2025.
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------
 */

#ifndef TEMPLAT_UTIL_SPLINE_H
#define TEMPLAT_UTIL_SPLINE_H

#include <cassert>
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>

#include "TempLat/parallel/device_memory.h"

namespace TempLat
{

  // =====================================================================
  // SplineData — GPU-evaluatable handle (lightweight, copyable by value)
  // =====================================================================

  template <typename T> struct SplineData {
    const T *x = nullptr;
    const T *y = nullptr;
    const T *b = nullptr;
    const T *c = nullptr;
    const T *d = nullptr;
    T c0 = T(0);
    size_t n = 0;

    DEVICE_FUNCTION
    size_t find_closest(T val) const
    {
      // Manual binary search replacing std::upper_bound.
      // Returns largest idx such that x[idx] <= val (0 if val < x[0]).
      if (n == 0) return 0;
      size_t lo = 0;
      size_t hi = n;
      while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (x[mid] <= val)
          lo = mid + 1;
        else
          hi = mid;
      }
      return lo > 0 ? lo - 1 : 0;
    }

    DEVICE_FUNCTION
    T operator()(T val) const
    {
      size_t idx = find_closest(val);
      T h = val - x[idx];
      if (val < x[0]) {
        return (c0 * h + b[0]) * h + y[0];
      } else if (val > x[n - 1]) {
        return (c[n - 1] * h + b[n - 1]) * h + y[n - 1];
      } else {
        return ((d[idx] * h + c[idx]) * h + b[idx]) * h + y[idx];
      }
    }

    DEVICE_FUNCTION
    T deriv(int order, T val) const
    {
      size_t idx = find_closest(val);
      T h = val - x[idx];
      if (val < x[0]) {
        switch (order) {
        case 1:
          return T(2) * c0 * h + b[0];
        case 2:
          return T(2) * c0;
        default:
          return T(0);
        }
      } else if (val > x[n - 1]) {
        switch (order) {
        case 1:
          return T(2) * c[n - 1] * h + b[n - 1];
        case 2:
          return T(2) * c[n - 1];
        default:
          return T(0);
        }
      } else {
        switch (order) {
        case 1:
          return (T(3) * d[idx] * h + T(2) * c[idx]) * h + b[idx];
        case 2:
          return T(6) * d[idx] * h + T(2) * c[idx];
        case 3:
          return T(6) * d[idx];
        default:
          return T(0);
        }
      }
    }
  };

  // =====================================================================
  // detail — band matrix solver and cubic root helpers (host-only)
  // =====================================================================

  namespace detail
  {

    template <typename T> class band_matrix
    {
    private:
      std::vector<std::vector<T>> m_upper;
      std::vector<std::vector<T>> m_lower;

    public:
      band_matrix() {}
      band_matrix(int dim, int n_u, int n_l) { resize(dim, n_u, n_l); }

      void resize(int dim, int n_u, int n_l)
      {
        assert(dim > 0);
        assert(n_u >= 0);
        assert(n_l >= 0);
        m_upper.resize(n_u + 1);
        m_lower.resize(n_l + 1);
        for (size_t i = 0; i < m_upper.size(); i++)
          m_upper[i].resize(dim);
        for (size_t i = 0; i < m_lower.size(); i++)
          m_lower[i].resize(dim);
      }

      int dim() const { return m_upper.size() > 0 ? (int)m_upper[0].size() : 0; }
      int num_upper() const { return (int)m_upper.size() - 1; }
      int num_lower() const { return (int)m_lower.size() - 1; }

      T &operator()(int i, int j)
      {
        int k = j - i;
        assert((i >= 0) && (i < dim()) && (j >= 0) && (j < dim()));
        assert((-num_lower() <= k) && (k <= num_upper()));
        return k >= 0 ? m_upper[k][i] : m_lower[-k][i];
      }

      T operator()(int i, int j) const
      {
        int k = j - i;
        assert((i >= 0) && (i < dim()) && (j >= 0) && (j < dim()));
        assert((-num_lower() <= k) && (k <= num_upper()));
        return k >= 0 ? m_upper[k][i] : m_lower[-k][i];
      }

      T &saved_diag(int i)
      {
        assert((i >= 0) && (i < dim()));
        return m_lower[0][i];
      }

      T saved_diag(int i) const
      {
        assert((i >= 0) && (i < dim()));
        return m_lower[0][i];
      }

      void lu_decompose()
      {
        int i_max, j_max;
        int j_min;
        T x;
        for (int i = 0; i < this->dim(); i++) {
          assert(this->operator()(i, i) != T(0));
          this->saved_diag(i) = T(1) / this->operator()(i, i);
          j_min = std::max(0, i - this->num_lower());
          j_max = std::min(this->dim() - 1, i + this->num_upper());
          for (int j = j_min; j <= j_max; j++) {
            this->operator()(i, j) *= this->saved_diag(i);
          }
          this->operator()(i, i) = T(1);
        }
        for (int k = 0; k < this->dim(); k++) {
          i_max = std::min(this->dim() - 1, k + this->num_lower());
          for (int i = k + 1; i <= i_max; i++) {
            assert(this->operator()(k, k) != T(0));
            x = -this->operator()(i, k) / this->operator()(k, k);
            this->operator()(i, k) = -x;
            j_max = std::min(this->dim() - 1, k + this->num_upper());
            for (int j = k + 1; j <= j_max; j++) {
              this->operator()(i, j) = this->operator()(i, j) + x * this->operator()(k, j);
            }
          }
        }
      }

      std::vector<T> l_solve(const std::vector<T> &b) const
      {
        assert(this->dim() == (int)b.size());
        std::vector<T> x(this->dim());
        int j_start;
        T sum;
        for (int i = 0; i < this->dim(); i++) {
          sum = T(0);
          j_start = std::max(0, i - this->num_lower());
          for (int j = j_start; j < i; j++)
            sum += this->operator()(i, j) * x[j];
          x[i] = (b[i] * this->saved_diag(i)) - sum;
        }
        return x;
      }

      std::vector<T> r_solve(const std::vector<T> &b) const
      {
        assert(this->dim() == (int)b.size());
        std::vector<T> x(this->dim());
        int j_stop;
        T sum;
        for (int i = this->dim() - 1; i >= 0; i--) {
          sum = T(0);
          j_stop = std::min(this->dim() - 1, i + this->num_upper());
          for (int j = i + 1; j <= j_stop; j++)
            sum += this->operator()(i, j) * x[j];
          x[i] = (b[i] - sum) / this->operator()(i, i);
        }
        return x;
      }

      std::vector<T> lu_solve(const std::vector<T> &b, bool is_lu_decomposed = false)
      {
        assert(this->dim() == (int)b.size());
        if (!is_lu_decomposed) this->lu_decompose();
        auto y = this->l_solve(b);
        return this->r_solve(y);
      }
    };

    template <typename T> inline T get_eps() { return std::numeric_limits<T>::epsilon(); }

    template <typename T> inline std::vector<T> solve_linear(T a, T b)
    {
      std::vector<T> x;
      if (b == T(0)) {
        if (a == T(0)) {
          x.push_back(T(0));
        }
        return x;
      }
      x.push_back(-a / b);
      return x;
    }

    template <typename T> inline std::vector<T> solve_quadratic(T a, T b, T c, int newton_iter = 0)
    {
      if (c == T(0)) return solve_linear(a, b);
      T p = T(0.5) * b / c;
      T q = a / c;
      T discr = p * p - q;
      const T eps = T(0.5) * get_eps<T>();
      T discr_err = (T(6) * (p * p) + T(3) * std::fabs(q) + std::fabs(discr)) * eps;

      std::vector<T> x;
      if (std::fabs(discr) <= discr_err) {
        x.push_back(-p);
      } else if (discr > T(0)) {
        x.push_back(-p - std::sqrt(discr));
        x.push_back(-p + std::sqrt(discr));
      }
      for (size_t i = 0; i < x.size(); i++) {
        for (int k = 0; k < newton_iter; k++) {
          T f = (c * x[i] + b) * x[i] + a;
          T f1 = T(2) * c * x[i] + b;
          if (std::fabs(f1) > T(1e-8)) x[i] -= f / f1;
        }
      }
      return x;
    }

    template <typename T> inline std::vector<T> solve_cubic(T a, T b, T c, T d, int newton_iter = 0)
    {
      if (d == T(0)) return solve_quadratic(a, b, c, newton_iter);
      if (d != T(1)) {
        a /= d;
        b /= d;
        c /= d;
      }
      std::vector<T> z;
      T p = -(T(1) / T(3)) * b + (T(1) / T(9)) * (c * c);
      T r = T(2) * (c * c) - T(9) * b;
      T q = -T(0.5) * a - (T(1) / T(54)) * (c * r);
      T discr = p * p * p - q * q;
      const T eps = get_eps<T>();
      T p_err = eps * (std::fabs(b) + (T(4) / T(9)) * (c * c) + std::fabs(p));
      T r_err = eps * (T(6) * (c * c) + T(18) * std::fabs(b) + std::fabs(r));
      T q_err =
          T(0.5) * std::fabs(a) * eps + (T(1) / T(54)) * std::fabs(c) * (r_err + std::fabs(r) * T(3) * eps) + std::fabs(q) * eps;
      T discr_err =
          (p * p) * (T(3) * p_err + std::fabs(p) * T(2) * eps) + std::fabs(q) * (T(2) * q_err + std::fabs(q) * eps) + std::fabs(discr) * eps;

      if (std::fabs(discr) <= discr_err) {
        if (std::fabs(p) <= p_err) {
          z.push_back(T(0));
        } else {
          z.push_back(T(2) * q / p);
          z.push_back(-q / p);
        }
      } else if (discr > T(0)) {
        T ac = (T(1) / T(3)) * std::acos(q / (p * std::sqrt(p)));
        T sq = T(2) * std::sqrt(p);
        z.push_back(sq * std::cos(ac));
        z.push_back(sq * std::cos(ac - T(2) * T(M_PI) / T(3)));
        z.push_back(sq * std::cos(ac - T(4) * T(M_PI) / T(3)));
      } else {
        T sgnq = (q >= T(0) ? T(1) : T(-1));
        T basis = std::fabs(q) + std::sqrt(-discr);
        T C = sgnq * std::pow(basis, T(1) / T(3));
        z.push_back(C + p / C);
      }
      for (size_t i = 0; i < z.size(); i++) {
        z[i] -= (T(1) / T(3)) * c;
        for (int k = 0; k < newton_iter; k++) {
          T f = ((z[i] + c) * z[i] + b) * z[i] + a;
          T f1 = (T(3) * z[i] + T(2) * c) * z[i] + b;
          if (std::fabs(f1) > T(1e-8)) z[i] -= f / f1;
        }
      }
      if (a == T(0)) {
        assert(z.size() > 0);
        T xmin = std::fabs(z[0]);
        size_t imin = 0;
        for (size_t i = 1; i < z.size(); i++) {
          if (xmin > std::fabs(z[i])) {
            xmin = std::fabs(z[i]);
            imin = i;
          }
        }
        z[imin] = T(0);
      }
      std::sort(z.begin(), z.end());
      return z;
    }

  } // namespace detail

  // =====================================================================
  // Spline — host-side owner with device memory
  // =====================================================================

  template <typename T = double> class Spline
  {
  public:
    enum SplineType { linear = 10, cspline = 30, cspline_hermite = 31 };

    enum BdType { first_deriv = 1, second_deriv = 2, not_a_knot = 3 };

    Spline()
        : m_type(cspline), m_left(second_deriv), m_right(second_deriv), m_left_value(T(0)), m_right_value(T(0)),
          m_made_monotonic(false)
    {
    }

    Spline(const std::vector<T> &X, const std::vector<T> &Y, SplineType type = cspline,
           bool make_monotonic = false, BdType left = second_deriv, T left_value = T(0),
           BdType right = second_deriv, T right_value = T(0))
        : m_type(type), m_left(left), m_right(right), m_left_value(left_value), m_right_value(right_value),
          m_made_monotonic(false)
    {
      this->set_points(X, Y, m_type);
      if (make_monotonic) this->make_monotonic();
    }

    void set_boundary(BdType left, T left_value, BdType right, T right_value)
    {
      assert(m_x.size() == 0);
      m_left = left;
      m_right = right;
      m_left_value = left_value;
      m_right_value = right_value;
    }

    void set_points(const std::vector<T> &x, const std::vector<T> &y, SplineType type = cspline)
    {
      assert(x.size() == y.size());
      assert(x.size() >= 3);
      if (m_left == not_a_knot || m_right == not_a_knot) assert(x.size() >= 4);
      m_type = type;
      m_made_monotonic = false;
      m_x = x;
      m_y = y;
      int n = (int)x.size();
      for (int i = 0; i < n - 1; i++)
        assert(m_x[i] < m_x[i + 1]);

      if (type == linear) {
        m_d.resize(n);
        m_c.resize(n);
        m_b.resize(n);
        for (int i = 0; i < n - 1; i++) {
          m_d[i] = T(0);
          m_c[i] = T(0);
          m_b[i] = (m_y[i + 1] - m_y[i]) / (m_x[i + 1] - m_x[i]);
        }
        m_b[n - 1] = m_b[n - 2];
        m_c[n - 1] = T(0);
        m_d[n - 1] = T(0);
      } else if (type == cspline) {
        int n_upper = (m_left == not_a_knot) ? 2 : 1;
        int n_lower = (m_right == not_a_knot) ? 2 : 1;
        detail::band_matrix<T> A(n, n_upper, n_lower);
        std::vector<T> rhs(n);
        for (int i = 1; i < n - 1; i++) {
          A(i, i - 1) = (T(1) / T(3)) * (x[i] - x[i - 1]);
          A(i, i) = (T(2) / T(3)) * (x[i + 1] - x[i - 1]);
          A(i, i + 1) = (T(1) / T(3)) * (x[i + 1] - x[i]);
          rhs[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
        }
        if (m_left == second_deriv) {
          A(0, 0) = T(2);
          A(0, 1) = T(0);
          rhs[0] = m_left_value;
        } else if (m_left == first_deriv) {
          A(0, 0) = T(2) * (x[1] - x[0]);
          A(0, 1) = x[1] - x[0];
          rhs[0] = T(3) * ((y[1] - y[0]) / (x[1] - x[0]) - m_left_value);
        } else if (m_left == not_a_knot) {
          A(0, 0) = -(x[2] - x[1]);
          A(0, 1) = x[2] - x[0];
          A(0, 2) = -(x[1] - x[0]);
          rhs[0] = T(0);
        } else {
          assert(false);
        }
        if (m_right == second_deriv) {
          A(n - 1, n - 1) = T(2);
          A(n - 1, n - 2) = T(0);
          rhs[n - 1] = m_right_value;
        } else if (m_right == first_deriv) {
          A(n - 1, n - 1) = T(2) * (x[n - 1] - x[n - 2]);
          A(n - 1, n - 2) = x[n - 1] - x[n - 2];
          rhs[n - 1] = T(3) * (m_right_value - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
        } else if (m_right == not_a_knot) {
          A(n - 1, n - 3) = -(x[n - 1] - x[n - 2]);
          A(n - 1, n - 2) = x[n - 1] - x[n - 3];
          A(n - 1, n - 1) = -(x[n - 2] - x[n - 3]);
          rhs[0] = T(0);
        } else {
          assert(false);
        }
        m_c = A.lu_solve(rhs);
        m_d.resize(n);
        m_b.resize(n);
        for (int i = 0; i < n - 1; i++) {
          m_d[i] = (T(1) / T(3)) * (m_c[i + 1] - m_c[i]) / (x[i + 1] - x[i]);
          m_b[i] =
              (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (T(1) / T(3)) * (T(2) * m_c[i] + m_c[i + 1]) * (x[i + 1] - x[i]);
        }
        T h = x[n - 1] - x[n - 2];
        m_d[n - 1] = T(0);
        m_b[n - 1] = T(3) * m_d[n - 2] * h * h + T(2) * m_c[n - 2] * h + m_b[n - 2];
        if (m_right == first_deriv) m_c[n - 1] = T(0);
      } else if (type == cspline_hermite) {
        m_b.resize(n);
        m_c.resize(n);
        m_d.resize(n);
        for (int i = 1; i < n - 1; i++) {
          const T h = m_x[i + 1] - m_x[i];
          const T hl = m_x[i] - m_x[i - 1];
          m_b[i] = -h / (hl * (hl + h)) * m_y[i - 1] + (h - hl) / (hl * h) * m_y[i] + hl / (h * (hl + h)) * m_y[i + 1];
        }
        if (m_left == first_deriv) {
          m_b[0] = m_left_value;
        } else if (m_left == second_deriv) {
          const T h = m_x[1] - m_x[0];
          m_b[0] = T(0.5) * (-m_b[1] - T(0.5) * m_left_value * h + T(3) * (m_y[1] - m_y[0]) / h);
        } else if (m_left == not_a_knot) {
          const T h0 = m_x[1] - m_x[0];
          const T h1 = m_x[2] - m_x[1];
          m_b[0] = -m_b[1] + T(2) * (m_y[1] - m_y[0]) / h0 +
                   h0 * h0 / (h1 * h1) * (m_b[1] + m_b[2] - T(2) * (m_y[2] - m_y[1]) / h1);
        } else {
          assert(false);
        }
        if (m_right == first_deriv) {
          m_b[n - 1] = m_right_value;
          m_c[n - 1] = T(0);
        } else if (m_right == second_deriv) {
          const T h = m_x[n - 1] - m_x[n - 2];
          m_b[n - 1] = T(0.5) * (-m_b[n - 2] + T(0.5) * m_right_value * h + T(3) * (m_y[n - 1] - m_y[n - 2]) / h);
          m_c[n - 1] = T(0.5) * m_right_value;
        } else if (m_right == not_a_knot) {
          const T h0 = m_x[n - 2] - m_x[n - 3];
          const T h1 = m_x[n - 1] - m_x[n - 2];
          m_b[n - 1] = -m_b[n - 2] + T(2) * (m_y[n - 1] - m_y[n - 2]) / h1 +
                       h1 * h1 / (h0 * h0) * (m_b[n - 3] + m_b[n - 2] - T(2) * (m_y[n - 2] - m_y[n - 3]) / h0);
          m_c[n - 1] = (m_b[n - 2] + T(2) * m_b[n - 1]) / h1 - T(3) * (m_y[n - 1] - m_y[n - 2]) / (h1 * h1);
        } else {
          assert(false);
        }
        m_d[n - 1] = T(0);
        set_coeffs_from_b();
      } else {
        assert(false);
      }

      m_c0 = (m_left == first_deriv) ? T(0) : m_c[0];

      // Push data to device memory
      push();
    }

    bool make_monotonic()
    {
      assert(m_x.size() == m_y.size());
      assert(m_x.size() == m_b.size());
      assert(m_x.size() > 2);
      bool modified = false;
      const int n = (int)m_x.size();
      for (int i = 0; i < n; i++) {
        int im1 = std::max(i - 1, 0);
        int ip1 = std::min(i + 1, n - 1);
        if (((m_y[im1] <= m_y[i]) && (m_y[i] <= m_y[ip1]) && m_b[i] < T(0)) ||
            ((m_y[im1] >= m_y[i]) && (m_y[i] >= m_y[ip1]) && m_b[i] > T(0))) {
          modified = true;
          m_b[i] = T(0);
        }
      }
      for (int i = 0; i < n - 1; i++) {
        T h = m_x[i + 1] - m_x[i];
        T avg = (m_y[i + 1] - m_y[i]) / h;
        if (avg == T(0) && (m_b[i] != T(0) || m_b[i + 1] != T(0))) {
          modified = true;
          m_b[i] = T(0);
          m_b[i + 1] = T(0);
        } else if ((m_b[i] >= T(0) && m_b[i + 1] >= T(0) && avg > T(0)) ||
                   (m_b[i] <= T(0) && m_b[i + 1] <= T(0) && avg < T(0))) {
          T r = std::sqrt(m_b[i] * m_b[i] + m_b[i + 1] * m_b[i + 1]) / std::fabs(avg);
          if (r > T(3)) {
            modified = true;
            m_b[i] *= (T(3) / r);
            m_b[i + 1] *= (T(3) / r);
          }
        }
      }
      if (modified) {
        set_coeffs_from_b();
        m_made_monotonic = true;
      }
      return modified;
    }

    // Evaluate on device (requires push() to have been called)
    DEVICE_INLINE_FUNCTION
    T operator()(T x) const { return data()(x); }

    // Host-side evaluation (uses std::vector data directly, no push() needed)
    T eval_host(T x) const
    {
      size_t n = m_x.size();
      size_t idx = find_closest(x);
      T h = x - m_x[idx];
      if (x < m_x[0]) {
        return (m_c0 * h + m_b[0]) * h + m_y[0];
      } else if (x > m_x[n - 1]) {
        return (m_c[n - 1] * h + m_b[n - 1]) * h + m_y[n - 1];
      } else {
        return ((m_d[idx] * h + m_c[idx]) * h + m_b[idx]) * h + m_y[idx];
      }
    }

    // Evaluate derivative on device (requires push() to have been called)
    DEVICE_INLINE_FUNCTION
    T deriv(int order, T x) const { return data().deriv(order, x); }

    // Host-side derivative evaluation (uses std::vector data directly, no push() needed)
    T deriv_host(int order, T x) const
    {
      assert(order > 0);
      size_t n = m_x.size();
      size_t idx = find_closest(x);
      T h = x - m_x[idx];
      T interpol;
      if (x < m_x[0]) {
        switch (order) {
        case 1:
          interpol = T(2) * m_c0 * h + m_b[0];
          break;
        case 2:
          interpol = T(2) * m_c0;
          break;
        default:
          interpol = T(0);
          break;
        }
      } else if (x > m_x[n - 1]) {
        switch (order) {
        case 1:
          interpol = T(2) * m_c[n - 1] * h + m_b[n - 1];
          break;
        case 2:
          interpol = T(2) * m_c[n - 1];
          break;
        default:
          interpol = T(0);
          break;
        }
      } else {
        switch (order) {
        case 1:
          interpol = (T(3) * m_d[idx] * h + T(2) * m_c[idx]) * h + m_b[idx];
          break;
        case 2:
          interpol = T(6) * m_d[idx] * h + T(2) * m_c[idx];
          break;
        case 3:
          interpol = T(6) * m_d[idx];
          break;
        default:
          interpol = T(0);
          break;
        }
      }
      return interpol;
    }

    std::vector<T> solve(T y, bool ignore_extrapolation = true) const
    {
      std::vector<T> x;
      std::vector<T> root;
      const size_t n = m_x.size();
      if (!ignore_extrapolation) {
        root = detail::solve_cubic<T>(m_y[0] - y, m_b[0], m_c0, T(0), 1);
        for (size_t j = 0; j < root.size(); j++) {
          if (root[j] < T(0)) x.push_back(m_x[0] + root[j]);
        }
      }
      for (size_t i = 0; i < n - 1; i++) {
        root = detail::solve_cubic<T>(m_y[i] - y, m_b[i], m_c[i], m_d[i], 1);
        for (size_t j = 0; j < root.size(); j++) {
          T h = (i > 0) ? (m_x[i] - m_x[i - 1]) : T(0);
          T eps = detail::get_eps<T>() * T(512) * std::min(h, T(1));
          if ((-eps <= root[j]) && (root[j] < m_x[i + 1] - m_x[i])) {
            T new_root = m_x[i] + root[j];
            if (x.size() > 0 && x.back() + eps > new_root) {
              x.back() = new_root;
            } else {
              x.push_back(new_root);
            }
          }
        }
      }
      if (!ignore_extrapolation) {
        root = detail::solve_cubic<T>(m_y[n - 1] - y, m_b[n - 1], m_c[n - 1], T(0), 1);
        for (size_t j = 0; j < root.size(); j++) {
          if (T(0) <= root[j]) x.push_back(m_x[n - 1] + root[j]);
        }
      }
      return x;
    }

    // Copy coefficient arrays to device and cache the SplineData handle
    void push()
    {
      size_t n = m_x.size();
      m_x_device = device::memory::NDView<T, 1>("spline_x", n);
      m_y_device = device::memory::NDView<T, 1>("spline_y", n);
      m_b_device = device::memory::NDView<T, 1>("spline_b", n);
      m_c_device = device::memory::NDView<T, 1>("spline_c", n);
      m_d_device = device::memory::NDView<T, 1>("spline_d", n);
      device::memory::copyHostToDevice(m_x.data(), m_x_device);
      device::memory::copyHostToDevice(m_y.data(), m_y_device);
      device::memory::copyHostToDevice(m_b.data(), m_b_device);
      device::memory::copyHostToDevice(m_c.data(), m_c_device);
      device::memory::copyHostToDevice(m_d.data(), m_d_device);
      m_data = SplineData<T>{
          m_x_device.data(), m_y_device.data(), m_b_device.data(), m_c_device.data(), m_d_device.data(), m_c0, n};
    }

    // Return the cached GPU-evaluatable handle (requires previous push())
    DEVICE_INLINE_FUNCTION
    const SplineData<T> &data() const { return m_data; }

    std::vector<T> get_x() const { return m_x; }
    std::vector<T> get_y() const { return m_y; }
    T get_x_min() const
    {
      assert(!m_x.empty());
      return m_x.front();
    }
    T get_x_max() const
    {
      assert(!m_x.empty());
      return m_x.back();
    }

  private:
    // Host-side coefficient storage
    std::vector<T> m_x, m_y;
    std::vector<T> m_b, m_c, m_d;
    T m_c0 = T(0);
    SplineType m_type;
    BdType m_left, m_right;
    T m_left_value, m_right_value;
    bool m_made_monotonic;

    // Cached device handle
    SplineData<T> m_data;

    // Device-side views (Kokkos::View is reference-counted)
    device::memory::NDView<T, 1> m_x_device;
    device::memory::NDView<T, 1> m_y_device;
    device::memory::NDView<T, 1> m_b_device;
    device::memory::NDView<T, 1> m_c_device;
    device::memory::NDView<T, 1> m_d_device;

    void set_coeffs_from_b()
    {
      assert(m_x.size() == m_y.size());
      assert(m_x.size() == m_b.size());
      assert(m_x.size() > 2);
      size_t n = m_b.size();
      if (m_c.size() != n) m_c.resize(n);
      if (m_d.size() != n) m_d.resize(n);
      for (size_t i = 0; i < n - 1; i++) {
        const T h = m_x[i + 1] - m_x[i];
        m_c[i] = (T(3) * (m_y[i + 1] - m_y[i]) / h - (T(2) * m_b[i] + m_b[i + 1])) / h;
        m_d[i] = ((m_b[i + 1] - m_b[i]) / (T(3) * h) - (T(2) / T(3)) * m_c[i]) / h;
      }
      m_c0 = (m_left == first_deriv) ? T(0) : m_c[0];
    }

    size_t find_closest(T x) const
    {
      auto it = std::upper_bound(m_x.begin(), m_x.end(), x);
      return std::max(int(it - m_x.begin()) - 1, 0);
    }
  };

} // namespace TempLat

#endif /* TEMPLAT_UTIL_SPLINE_H */
