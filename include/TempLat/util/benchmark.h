#ifndef TEMPLAT_UTIL_TempLatBENCHMARK_H
#define TEMPLAT_UTIL_TempLatBENCHMARK_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler,  Year: 2025

#include <cstddef>
#include <functional>
#include <limits>
#include <sstream>
#include <iomanip>
#include <list>
#include <fstream>

#include "TempLat/util/log/saycomplete.h"
#include "TempLat/util/timer.h"
#include "TempLat/parallel/device_iteration.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace TempLat
{
  /**
   * @brief A class to benchmark your code.
   *
   */
  class Benchmark
  {
  public:
    class Measurer
    {
    public:
      template <typename F> void measure(const std::string tag, F &&f)
      {
        device::iteration::fence();
#ifdef HAVE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        const Timer timer;
        f();
        device::iteration::fence();
#ifdef HAVE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        const size_t elapsed = timer.nanoseconds();
        measurements.emplace_back(tag, elapsed);
      }

      void collect()
      {
#ifdef HAVE_MPI
        int world_size, world_rank;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

        // Serialize measurements: flatten to tag lengths, tags, and elapsed times
        std::vector<int> tag_lengths;
        std::vector<char> tag_chars;
        std::vector<double> elapsed_times;
        for (const auto &m : measurements) {
          tag_lengths.push_back(m.first.size());
          tag_chars.insert(tag_chars.end(), m.first.begin(), m.first.end());
          elapsed_times.push_back(m.second);
        }
        int num_measurements = measurements.size();

        // Gather sizes at root
        std::vector<int> recv_counts(world_size);
        MPI_Gather(&num_measurements, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

        // Gather tag lengths
        std::vector<int> tag_len_recv;
        if (world_rank == 0) tag_len_recv.resize(world_size * num_measurements); // over-allocate
        MPI_Gather(tag_lengths.data(), num_measurements, MPI_INT, world_rank == 0 ? tag_len_recv.data() : nullptr,
                   num_measurements, MPI_INT, 0, MPI_COMM_WORLD);

        // Gather tags
        int tag_chars_len = tag_chars.size();
        std::vector<int> tag_chars_counts(world_size);
        MPI_Gather(&tag_chars_len, 1, MPI_INT, tag_chars_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
        std::vector<char> tag_chars_recv;
        int total_tag_chars = 0;
        if (world_rank == 0) {
          for (int c : tag_chars_counts)
            total_tag_chars += c;
          tag_chars_recv.resize(total_tag_chars);
        }
        std::vector<int> tag_chars_displs(world_size);
        if (world_rank == 0) {
          int offset = 0;
          for (int i = 0; i < world_size; ++i) {
            tag_chars_displs[i] = offset;
            offset += tag_chars_counts[i];
          }
        }
        MPI_Gatherv(tag_chars.data(), tag_chars_len, MPI_CHAR, world_rank == 0 ? tag_chars_recv.data() : nullptr,
                    world_rank == 0 ? tag_chars_counts.data() : nullptr,
                    world_rank == 0 ? tag_chars_displs.data() : nullptr, MPI_CHAR, 0, MPI_COMM_WORLD);

        // Gather elapsed times
        std::vector<double> elapsed_recv;
        int total_measurements = 0;
        if (world_rank == 0) {
          for (int c : recv_counts)
            total_measurements += c;
          elapsed_recv.resize(total_measurements);
        }
        std::vector<int> elapsed_displs(world_size);
        if (world_rank == 0) {
          int offset = 0;
          for (int i = 0; i < world_size; ++i) {
            elapsed_displs[i] = offset;
            offset += recv_counts[i];
          }
        }
        MPI_Gatherv(elapsed_times.data(), num_measurements, MPI_DOUBLE, world_rank == 0 ? elapsed_recv.data() : nullptr,
                    world_rank == 0 ? recv_counts.data() : nullptr, world_rank == 0 ? elapsed_displs.data() : nullptr,
                    MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Reconstruct all_measurements at root
        if (world_rank == 0) {
          std::vector<std::pair<std::string, double>> all_measurements;
          int tag_pos = 0;
          int elapsed_pos = 0;
          for (int i = 0; i < world_size; ++i) {
            for (int j = 0; j < recv_counts[i]; ++j) {
              int len = tag_len_recv[i * num_measurements + j];
              std::string tag(tag_chars_recv.begin() + tag_pos, tag_chars_recv.begin() + tag_pos + len);
              tag_pos += len;
              double elapsed = elapsed_recv[elapsed_pos++];
              all_measurements.emplace_back(tag, elapsed);
            }
          }
          measurements = all_measurements;
        }
#else
        // No MPI: nothing to do
#endif
      }

      double getAverage(const std::string &tag) const
      {
        double total = 0;
        size_t count = 0;
        for (const auto &[measurementTag, elapsed] : measurements) {
          if (measurementTag == tag) {
            total += elapsed;
            ++count;
          }
        }
        return count > 0 ? total / count : 0;
      }

      auto getMeasurement(const std::string tag) const
      {
        double average = getAverage(tag);
        double totalSquaredDiff = 0;
        double count = 0;
        double min = std::numeric_limits<double>::max();
        double max = 0;
        for (const auto &[measurementTag, elapsed] : measurements) {
          if (measurementTag == tag) {
            double diff = elapsed > average ? elapsed - average : average - elapsed;
            totalSquaredDiff += diff * diff;
            ++count;
            min = std::min(min, elapsed);
            max = std::max(max, elapsed);
          }
        }
        double stdD = count > 0 ? std::sqrt(totalSquaredDiff) / count : 0;

        return std::make_tuple(tag, average, stdD, min, max, count);
      }

      std::vector<std::string> getTags() const
      {
        std::vector<std::string> tags;
        for (const auto &[tag, _] : measurements) {
          if (std::find(tags.begin(), tags.end(), tag) == tags.end()) {
            tags.push_back(tag);
          }
        }
        return tags;
      }

    private:
      std::vector<std::pair<std::string, double>> measurements; // tag, elapsed time in nanoseconds
      friend class Benchmark;
    };

    template <typename F> Benchmark(F &&function) : mFunction(std::forward<F>(function)) {}

    void run(size_t n = 0)
    {
      int rank = 0;
#ifdef HAVE_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

      Measurer dead_measurer;
      if (n == 0) {
        if (rank == 0) sayMPI << "Estimating number of iterations to run for the benchmark.\n";
        // If n is 0, we run the function once and check how long it takes.
        Timer timer;
        mFunction(dead_measurer);
        size_t elapsed = timer.nanoseconds();
        // We wish to run the benchmark for no longer than 30s
        n = 30'000'000'000 / elapsed;
        n = std::max(n, size_t(1));
        n = std::min(n, size_t(1000));
      }

      // warmup
      if (n >= 10) {
        if (rank == 0) sayMPI << "Running warmup for " << std::max((size_t)10, n / 2) << " iterations.\n";
        for (size_t i = 0; i < std::max((size_t)10, n / 2); ++i)
          mFunction(dead_measurer);
      }

#ifdef HAVE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif

      if (rank == 0) sayMPI << "Running benchmark for " << n << " iterations.\n";

      Measurer measurer;
      for (size_t i = 0; i < n; ++i)
        mFunction(measurer);

      measurer.collect();

      for (const auto &tag : measurer.getTags())
        mMeasurements[tag] = measurer.getMeasurement(tag);
    }

    /**
     * @brief Creates a unique log file in the current directory with the benchmark results.
     *
     */
    void log(std::string name) const
    {
      std::string formatted_time = std::to_string(std::time(nullptr));
      std::string filename = name + "_bench_" + formatted_time + ".csv";

      std::ofstream logFile(filename, std::ios::app);
      if (logFile.is_open()) {
        logFile << "Tag,Average[s],StdDev[s],Count\n";
        for (const auto &el : mMeasurements) {
          auto [measurementTag_, measurementData] = el;
          std::string measurementTag = "\"" + measurementTag_ + "\"";
          // replace whitespace in the tag with underscores
          for (auto &c : measurementTag)
            if (std::isspace(c)) c = '_';

          const auto &[tag, average, stdDev, minT, maxT, count] = measurementData;

          const double averageInSeconds = average / 1e9;
          const double stdDevInSeconds = stdDev / 1e9;

          logFile << measurementTag << "," << averageInSeconds << "," << stdDevInSeconds << "," << count << "\n";
        }
        logFile.close();
      } else {
        sayMPI << "Could not open " << filename << " for writing.\n";
      }
    }

    void print() const
    {
#ifdef HAVE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
      int world_rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
      if (world_rank != 0) return; // Only print from the root process
#endif
      const size_t tagWidth = 28;
      std::list<std::string> outputs;
      for (const auto &el : mMeasurements) {
        const auto &[measurementTag, measurementData] = el;
        const auto &[tag, average, stdDev, minT, maxT, count] = measurementData;

        const auto [averageStr, averageLevel] = formatTime(average, 0);
        const auto [stdDevStr, stdDevLevel] = formatTime(stdDev, averageLevel);
        const auto [minStr, minLevel] = formatTime(minT, averageLevel);
        const auto [maxStr, maxLevel] = formatTime(maxT, averageLevel);
        const std::string countStr = std::to_string(count);

        std::stringstream ss;
        ss << "\033[32mTag:\033[1;4;34m" << std::setw(tagWidth - 4) << measurementTag << "\033[0m" // tag
           << "\n    \033[1;34m|\033[0m Average : " << std::setw(tagWidth - 16 - averageStr.size()) << ""
           << averageStr // average
           << "\n    \033[1;34m|\033[0m Std Dev : " << std::setw(tagWidth - 16 - stdDevStr.size()) << ""
           << stdDevStr                                                                                         // std
           << "\n    \033[1;34m|\033[0m     Min : " << std::setw(tagWidth - 16 - minStr.size()) << "" << minStr // std
           << "\n    \033[1;34m|\033[0m     Max : " << std::setw(tagWidth - 16 - maxStr.size()) << "" << maxStr // std
           << "\n    \033[1;34m|\033[0m Count :   " << std::setw(tagWidth - 16 - countStr.size()) << ""
           << countStr; // count
        outputs.push_back(ss.str());
      }

      const auto terminal_width(WEXITSTATUS(std::system("exit `tput cols`")));
      const size_t nextTo = (terminal_width - 8) / (tagWidth + 8); // 8 for the spaces in between

      // Now glue the outputs together
      std::string output;
      std::vector<std::string> curLines;
      while (!outputs.empty()) {
        curLines.clear();
        // Take the first nextTo outputs, split them into lines, and glue them together.
        for (size_t i = 0; i < nextTo && !outputs.empty(); ++i) {

          const std::string curOutput = outputs.front();
          outputs.pop_front();

          std::istringstream iss(curOutput);
          std::string line;
          size_t lineCount = 0;
          while (std::getline(iss, line)) {
            if (curLines.size() <= lineCount)
              curLines.push_back("    " + line);
            else
              curLines[lineCount] += "        " + line;
            ++lineCount;
          }
        }

        for (size_t i = 0; i < curLines.size(); ++i)
          output += curLines[i] + "\n";

        output += "\n";
      }

      std::cout << "\nBenchmark results:\n\n" << output;
    }

  private:
    std::function<void(Measurer &)> mFunction;

    using Measurement =
        std::tuple<std::string, double, double, double, double, size_t>; // tag, average, std, min, max, count
    std::map<std::string, Measurement> mMeasurements;

    static std::pair<std::string, int> formatTime(size_t time_ns, int lv = 0)
    {
      size_t total = time_ns;
      const size_t nanoseconds = total % 1000;
      total /= 1000;

      const size_t micro = total % 1000;
      total /= 1000;

      const size_t milli = total % 1000;
      total /= 1000;

      const size_t sec = total % 60;
      total /= 60;

      const size_t min = total % 60;
      total /= 60;

      const size_t hours = total;

      if ((lv == 0 && hours > 0) || (lv == 1 && hours > 0))
        return make_pair(std::to_string(hours) + "h " + std::to_string(min) + "min", 1);
      else if ((lv == 0 && min > 0) || (lv == 2 && min > 0))
        return make_pair(std::to_string(min) + "min " + std::to_string(sec) + "s", 2);
      else if ((lv == 0 && sec > 0) || (lv == 3 && sec > 0))
        return make_pair(std::to_string(sec) + "s " + std::to_string(milli) + "ms", 3);
      else if ((lv == 0 && milli > 0) || (lv == 4 && milli > 0))
        return make_pair(std::to_string(milli) + "ms " + std::to_string(micro) + "us", 4);
      else if ((lv == 0 && micro > 0) || (lv == 5 && micro > 0))
        return make_pair(std::to_string(micro) + "us " + std::to_string(nanoseconds) + "ns", 5);
      else if (lv == 1)
        return make_pair(std::to_string(min) + "min", -2);
      else if (lv == 2)
        return make_pair(std::to_string(sec) + "s", -3);
      else if (lv == 3)
        return make_pair(std::to_string(milli) + "ms", -4);
      else if (lv == 4)
        return make_pair(std::to_string(micro) + "us", -5);
      else if (lv == 5 || lv == 0)
        return make_pair(std::to_string(nanoseconds) + "ns", 6);
      throw std::runtime_error("Invalid level for formatting time: " + std::to_string(lv));
    }
  };
} // namespace TempLat

#endif