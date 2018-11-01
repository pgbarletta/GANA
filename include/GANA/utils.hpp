#ifndef GANA_UTILS
#define GANA_UTILS

extern float rsltion;
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <thrust/device_vector.h>

namespace GANA {

    // Turn a coordinate in the matching matrix grid index.
    inline int32_t cont_to_grid(float x) {
        return static_cast<int32_t>(fabs(x - fmod(x, rsltion)) / rsltion);
    }

    // Turn a grid index into a xyz coordinate.
    inline float grid_to_cont(int32_t idx) {
        return static_cast<float>(idx * rsltion);
    }

	// Helper function for getting the indices that sort a cont_vector.
	template <typename T>
	std::vector<unsigned int> sort_indices(const std::vector<T> &v);

	// Helper function to get the indices of the true elements of a bool array.
	// Optimized for large (>500) and sparse bool arrays.
	void get_indices_from_sparse_bool_array(bool *in_array, const int n_in,
	std::vector<int> &indices);

	// Template structure to pass to kernel
	template <typename T>
	struct Karray
	{
	    T*  array_;
	    int size_;
	};

	// Function to convert thrust::device_vector to Karray struct
	template <typename T>
	Karray<T> to_kernel(thrust::device_vector<T>& Darray)
	{
		Karray<T> KDarray;
	    KDarray.array_ = thrust::raw_pointer_cast(Darray.data());
	    KDarray.size_ = static_cast<int>(Darray.size());

	    return KDarray;
	}

} // namespace GANA

#endif // GANA_UTILS
