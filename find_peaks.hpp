#ifndef FIND_PEAKS_HPP
#define FIND_PEAKS_HPP

#include <vector>
#include <limits>
#include <algorithm>
#include <cmath>
#include <numeric> // For std::iota

namespace signal_processing {

/**
 * @brief Options for peak finding, mirroring scipy.signal.find_peaks parameters.
 */
struct FindPeaksOptions {
    // Height of peaks: (min_height, max_height)
    std::pair<double, double> height = { -std::numeric_limits<double>::infinity(),
                                          std::numeric_limits<double>::infinity() };

    // Threshold: (min_threshold, max_threshold)
    std::pair<double, double> threshold = { -std::numeric_limits<double>::infinity(),
                                             std::numeric_limits<double>::infinity() };

    // Minimum distance between peaks
    int distance = 1;

    // Prominence: (min_prominence, max_prominence)
    std::pair<double, double> prominence = { -std::numeric_limits<double>::infinity(),
                                              std::numeric_limits<double>::infinity() };

    // Width: (min_width, max_width)
    std::pair<double, double> width = { -std::numeric_limits<double>::infinity(),
                                         std::numeric_limits<double>::infinity() };

    // Window length for prominence calculation
    int wlen = 0; // 0 indicates using the entire signal

    // Relative height for width calculation
    double rel_height = 0.5;

    // Plateau size: (min_plateau_size, max_plateau_size)
    std::pair<int, int> plateau_size = { 0, std::numeric_limits<int>::max() };
};

/**
 * @brief Structure to hold properties of detected peaks.
 */
struct PeakProperties {
    std::vector<double> peak_heights;
    std::vector<double> prominences;
    std::vector<int> left_bases;
    std::vector<int> right_bases;
    std::vector<double> widths;
    std::vector<double> width_heights;
    std::vector<double> left_ips;
    std::vector<double> right_ips;
    std::vector<int> plateau_sizes;
    std::vector<int> left_edges;
    std::vector<int> right_edges;
    std::vector<double> left_thresholds;
    std::vector<double> right_thresholds;
};

/**
 * @brief Finds peaks in a 1D signal array with properties similar to scipy.signal.find_peaks.
 *
 * @param x The input signal as a vector of doubles.
 * @param peaks Output vector of indices where peaks are located in the input signal.
 * @param properties Output PeakProperties containing various peak attributes.
 * @param options A struct containing optional parameters for peak detection.
 */
void find_peaks(const std::vector<double>& x, std::vector<int>& peaks, PeakProperties& properties,
                const FindPeaksOptions& options = FindPeaksOptions()) {
    const size_t N = x.size();

    if (N < 3) {
        return; // Not enough data points to have peaks
    }

    std::vector<int> all_peaks;
    std::vector<int> left_edges;
    std::vector<int> right_edges;

    // Step 1: Find all local maxima and plateaus
    for (size_t i = 1; i < N - 1; ++i) {
        if (x[i - 1] < x[i] && x[i] > x[i + 1]) {
            // Sharp peak
            all_peaks.push_back(static_cast<int>(i));
            left_edges.push_back(static_cast<int>(i));
            right_edges.push_back(static_cast<int>(i));
        } else if (x[i - 1] < x[i] && x[i] == x[i + 1]) {
            // Start of plateau
            size_t j = i + 1;
            while (j + 1 < N && x[j] == x[j + 1]) {
                ++j;
            }
            if (j + 1 < N && x[j] > x[j + 1]) {
                // Plateau peak
                int peak_idx = static_cast<int>((i + j) / 2);
                all_peaks.push_back(peak_idx);
                left_edges.push_back(static_cast<int>(i));
                right_edges.push_back(static_cast<int>(j));
                i = j;
            }
        }
    }

    // Initialize properties
    PeakProperties props;

    // Store initial properties
    for (size_t idx = 0; idx < all_peaks.size(); ++idx) {
        int peak_idx = all_peaks[idx];
        props.peak_heights.push_back(x[peak_idx]);
        props.left_edges.push_back(left_edges[idx]);
        props.right_edges.push_back(right_edges[idx]);
        int plateau_size = right_edges[idx] - left_edges[idx] + 1;
        props.plateau_sizes.push_back(plateau_size);
    }

    // At each step, we'll maintain indices of the peaks that satisfy the conditions
    std::vector<int> keep_indices(all_peaks.size());
    std::iota(keep_indices.begin(), keep_indices.end(), 0); // Initialize with all indices

    // Helper lambda to filter peaks and properties based on keep_indices
    auto filter_peaks_and_properties = [&](const std::vector<int>& indices) {
        std::vector<int> filtered_peaks;
        PeakProperties filtered_props;

        for (int idx : indices) {
            filtered_peaks.push_back(all_peaks[idx]);
            filtered_props.peak_heights.push_back(props.peak_heights[idx]);
            filtered_props.left_edges.push_back(props.left_edges[idx]);
            filtered_props.right_edges.push_back(props.right_edges[idx]);
            filtered_props.plateau_sizes.push_back(props.plateau_sizes[idx]);

            if (!props.left_thresholds.empty())
                filtered_props.left_thresholds.push_back(props.left_thresholds[idx]);
            if (!props.right_thresholds.empty())
                filtered_props.right_thresholds.push_back(props.right_thresholds[idx]);
            if (!props.prominences.empty())
                filtered_props.prominences.push_back(props.prominences[idx]);
            if (!props.left_bases.empty())
                filtered_props.left_bases.push_back(props.left_bases[idx]);
            if (!props.right_bases.empty())
                filtered_props.right_bases.push_back(props.right_bases[idx]);
            if (!props.widths.empty())
                filtered_props.widths.push_back(props.widths[idx]);
            if (!props.width_heights.empty())
                filtered_props.width_heights.push_back(props.width_heights[idx]);
            if (!props.left_ips.empty())
                filtered_props.left_ips.push_back(props.left_ips[idx]);
            if (!props.right_ips.empty())
                filtered_props.right_ips.push_back(props.right_ips[idx]);
        }

        // Update all_peaks and properties
        all_peaks = filtered_peaks;
        props = filtered_props;
    };

    // Step 2: Apply plateau_size condition
    if (options.plateau_size.first > 0 || options.plateau_size.second < std::numeric_limits<int>::max()) {
        std::vector<int> new_keep_indices;
        for (int idx : keep_indices) {
            int plateau_size = props.plateau_sizes[idx];
            if (plateau_size >= options.plateau_size.first && plateau_size <= options.plateau_size.second) {
                new_keep_indices.push_back(idx);
            }
        }
        keep_indices = new_keep_indices;
        filter_peaks_and_properties(keep_indices);
    }

    // Step 3: Apply height condition
    if (options.height.first > -std::numeric_limits<double>::infinity() ||
        options.height.second < std::numeric_limits<double>::infinity()) {
        std::vector<int> new_keep_indices;
        for (int idx : keep_indices) {
            double peak_height = props.peak_heights[idx];
            if (peak_height >= options.height.first && peak_height <= options.height.second) {
                new_keep_indices.push_back(idx);
            }
        }
        keep_indices = new_keep_indices;
        filter_peaks_and_properties(keep_indices);
    }

    // Step 4: Apply threshold condition
    if (options.threshold.first > -std::numeric_limits<double>::infinity() ||
        options.threshold.second < std::numeric_limits<double>::infinity()) {
        // Compute thresholds
        props.left_thresholds.clear();
        props.right_thresholds.clear();
        for (int idx : keep_indices) {
            int peak_idx = all_peaks[idx];
            double left_thresh = peak_idx > 0 ? x[peak_idx] - x[peak_idx - 1] : std::numeric_limits<double>::infinity();
            double right_thresh = peak_idx + 1 < N ? x[peak_idx] - x[peak_idx + 1] : std::numeric_limits<double>::infinity();
            props.left_thresholds.push_back(left_thresh);
            props.right_thresholds.push_back(right_thresh);
        }

        std::vector<int> new_keep_indices;
        for (size_t i = 0; i < keep_indices.size(); ++i) {
            double left_thresh = props.left_thresholds[i];
            double right_thresh = props.right_thresholds[i];
            if (left_thresh >= options.threshold.first && left_thresh <= options.threshold.second &&
                right_thresh >= options.threshold.first && right_thresh <= options.threshold.second) {
                new_keep_indices.push_back(keep_indices[i]);
            }
        }
        keep_indices = new_keep_indices;
        filter_peaks_and_properties(keep_indices);
    }

    // Step 5: Apply distance condition
    if (options.distance > 1 && all_peaks.size() > 1) {
        // Sort peaks by their amplitude in descending order
        std::vector<std::pair<double, int>> peak_heights;
        for (size_t i = 0; i < all_peaks.size(); ++i) {
            peak_heights.emplace_back(props.peak_heights[i], static_cast<int>(i));
        }
        std::sort(peak_heights.begin(), peak_heights.end(),
                  [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                      return a.first > b.first;
                  });

        std::vector<bool> occupied(N, false);
        std::vector<int> new_keep_indices;

        for (const auto& p : peak_heights) {
            int idx = p.second;
            int peak_idx = all_peaks[idx];
            int left = std::max(0, peak_idx - options.distance);
            int right = std::min(static_cast<int>(N) - 1, peak_idx + options.distance);

            bool is_occupied = false;
            for (int i = left; i <= right; ++i) {
                if (occupied[i]) {
                    is_occupied = true;
                    break;
                }
            }
            if (!is_occupied) {
                new_keep_indices.push_back(idx);
                occupied[peak_idx] = true;
            }
        }

        // Sort new_keep_indices to maintain original order
        std::sort(new_keep_indices.begin(), new_keep_indices.end());
        keep_indices = new_keep_indices;
        filter_peaks_and_properties(keep_indices);
    }

    // Step 6: Calculate prominences if required
    if (options.prominence.first > -std::numeric_limits<double>::infinity() ||
        options.prominence.second < std::numeric_limits<double>::infinity() ||
        options.width.first > -std::numeric_limits<double>::infinity() ||
        options.width.second < std::numeric_limits<double>::infinity()) {
        // Compute prominences
        props.prominences.clear();
        props.left_bases.clear();
        props.right_bases.clear();

        for (size_t i = 0; i < all_peaks.size(); ++i) {
            int peak_idx = all_peaks[i];

            // Left base
            int left_base = peak_idx;
            double left_min = x[peak_idx];
            for (int j = peak_idx - 1; j >= 0; --j) {
                if (x[j] < left_min) {
                    left_min = x[j];
                    left_base = j;
                }
                if (x[j] > x[peak_idx]) {
                    break;
                }
            }

            // Right base
            int right_base = peak_idx;
            double right_min = x[peak_idx];
            for (int j = peak_idx + 1; j < static_cast<int>(N); ++j) {
                if (x[j] < right_min) {
                    right_min = x[j];
                    right_base = j;
                }
                if (x[j] > x[peak_idx]) {
                    break;
                }
            }

            double prominence = x[peak_idx] - std::max(left_min, right_min);
            props.prominences.push_back(prominence);
            props.left_bases.push_back(left_base);
            props.right_bases.push_back(right_base);
        }

        // Apply prominence condition
        if (options.prominence.first > -std::numeric_limits<double>::infinity() ||
            options.prominence.second < std::numeric_limits<double>::infinity()) {
            std::vector<int> new_keep_indices;
            for (size_t i = 0; i < all_peaks.size(); ++i) {
                double prominence = props.prominences[i];
                if (prominence >= options.prominence.first && prominence <= options.prominence.second) {
                    new_keep_indices.push_back(static_cast<int>(i));
                }
            }
            keep_indices = new_keep_indices;
            filter_peaks_and_properties(keep_indices);
        }
    }

    // Step 7: Calculate widths if required
    if (options.width.first > -std::numeric_limits<double>::infinity() ||
        options.width.second < std::numeric_limits<double>::infinity()) {
        props.widths.clear();
        props.width_heights.clear();
        props.left_ips.clear();
        props.right_ips.clear();

        for (size_t i = 0; i < all_peaks.size(); ++i) {
            int peak_idx = all_peaks[i];
            double prominence = props.prominences[i];
            double width_height = x[peak_idx] - options.rel_height * prominence;

            // Left interpolation
            double left_ip = static_cast<double>(peak_idx);
            for (int j = peak_idx; j >= props.left_bases[i]; --j) {
                if (x[j] <= width_height) {
                    if (j < peak_idx) {
                        left_ip = j + (width_height - x[j]) / (x[j + 1] - x[j]);
                    }
                    break;
                }
            }

            // Right interpolation
            double right_ip = static_cast<double>(peak_idx);
            for (int j = peak_idx; j <= props.right_bases[i]; ++j) {
                if (x[j] <= width_height) {
                    if (j > peak_idx) {
                        right_ip = j - (width_height - x[j]) / (x[j - 1] - x[j]);
                    }
                    break;
                }
            }

            double width = right_ip - left_ip;
            props.widths.push_back(width);
            props.width_heights.push_back(width_height);
            props.left_ips.push_back(left_ip);
            props.right_ips.push_back(right_ip);
        }

        // Apply width condition
        std::vector<int> new_keep_indices;
        for (size_t i = 0; i < all_peaks.size(); ++i) {
            double width = props.widths[i];
            if (width >= options.width.first && width <= options.width.second) {
                new_keep_indices.push_back(static_cast<int>(i));
            }
        }
        keep_indices = new_keep_indices;
        filter_peaks_and_properties(keep_indices);
    }

    // Final assignment
    peaks = all_peaks;
    properties = props;
}

} // namespace signal_processing

#endif // FIND_PEAKS_HPP

