#include <iostream>
#include "find_peaks.hpp"

int main() {
    std::vector<double> signal = {0, 2, 1, 3, 0, 1, 2, 1, 5, 7, 8, 1, 3, 0};

    signal_processing::FindPeaksOptions options;
    options.height = {1.5, std::numeric_limits<double>::infinity()};
    options.distance = 1;
    options.prominence = {0.5, std::numeric_limits<double>::infinity()};
    options.width = {0.0, std::numeric_limits<double>::infinity()};
    options.rel_height = 0.5;

    std::vector<int> peaks;
    signal_processing::PeakProperties properties;
    signal_processing::find_peaks(signal, peaks, properties, options);

    std::cout << "Peaks at indices: ";
    for (int idx : peaks) {
        std::cout << idx << " ";
    }
    std::cout << std::endl;

    // Output peak properties
    for (size_t i = 0; i < peaks.size(); ++i) {
        std::cout << "Peak at index " << peaks[i]
                  << " has height " << properties.peak_heights[i]
                  << ", prominence " << properties.prominences[i]
                  << ", width " << properties.widths[i]
                  << std::endl;
    }

    return 0;
}

