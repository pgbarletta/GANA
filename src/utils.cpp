#include "GANA/utils.hpp"
namespace GANA {

// Helper function for getting the indices that sort a vector.
template<typename T> std::vector<unsigned int> sort_indices(
	const std::vector<T> &v) {
	// initialize original index locations
	std::vector<unsigned int> idx(v.size());
	std::iota(idx.begin(), idx.end(), 0);

	// sort indices based on comparing values in v
	sort(idx.begin(), idx.end(),
			[&v](unsigned int i1, unsigned int i2) {return v[i1] < v[i2];});

	return idx;
}

// Helper function to get the indices of the true elements of a bool array.
// Optimized for large (>500) and sparse bool arrays.
void getIndicesFromSparseBoolArray(
	bool *in_array, const int n_in, std::vector<int> &indices) {

	auto sz_lo = sizeof(long);
	long *cin_array = (long *) in_array;
	int divi = n_in / sz_lo, resto = n_in % sz_lo;

	for (size_t i = 0; i < divi; ++i) {
		if (cin_array[i] != 0) {
			size_t lo = i * sz_lo, hi = lo + sz_lo;
			for (size_t j = lo; j < hi; ++j) {
				if (in_array[j] != 0) {
					indices.push_back(j);
				}
			}
		}
	}

	for (size_t i = n_in - resto; i < n_in; ++i) {
		if (in_array[i] != 0) {
			indices.push_back(i);
		}
	}

	return;
}
} // namespace GANA
