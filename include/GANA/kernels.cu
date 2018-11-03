inline void gpuAssert(cudaError_t code, char const *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}
#define GEC(ans) { gpuAssert((ans), __FILE__, __LINE__); } // GPU Error Check


// Kernels
__global__ void empiezo(GANA::Tetrahedron *in_IA, int const n) {

	int ti = threadIdx.x + blockIdx.x * blockDim.x;

	if (ti < n) {
        in_IA[ti]._p[0][0] += 10.;
        in_IA[ti]._p[0][1] += 10.;
        in_IA[ti]._p[0][2] += 10.;
        in_IA[ti]._p[1][0] += 10.;
        in_IA[ti]._p[1][1] += 10.;
        in_IA[ti]._p[1][2] += 10.;
        in_IA[ti]._p[2][0] += 10.;
        in_IA[ti]._p[2][1] += 10.;
        in_IA[ti]._p[2][2] += 10.;
        in_IA[ti]._p[3][0] += 10.;
        in_IA[ti]._p[3][1] += 10.;
        in_IA[ti]._p[3][2] += 10.;
	}

	return;
}

__global__ void init_grilla(float* grilla, int grilla_size, float x_min,
	float y_min, float z_min, float x_max, float y_max, float z_max,
	int x_cnt, int y_cnt, int z_cnt, float resolution) {

	const int nproc = gridDim.x * blockDim.x;
	int ti = threadIdx.x + blockIdx.x * blockDim.x;

	while (ti < grilla_size) {
		if (ti % 3 == 0) {
			float x_step = x_min + resolution * ((ti % (x_cnt * 3)) / 3);
			grilla[ti] = x_step;
		} else if (ti % 3 == 1) {
			float y_step = y_min
					+ resolution * ((ti % (x_cnt * y_cnt * 3)) / (3 * x_cnt));
			grilla[ti] = y_step;

		} else { // ti % 3 == 2
			float z_step = z_min + resolution * (ti / (3 * x_cnt * y_cnt));
			grilla[ti] = z_step;
		}
		ti = ti + nproc;
	}

	return;
}


__global__ void in_bbox(float* molecule_points, int x_min, int y_min, int z_min,
                int x_max, int y_max, int z_max, int natoms, bool* d_atoms_in_x,
                bool* d_atoms_in_y, bool* d_atoms_in_z, bool* atoms_in_bbox) {

        int ti = threadIdx.x + blockIdx.x * blockDim.x;

        // Get the atoms that lie inside the planes delimited by [xmin, xmax],
        // [ymin, ymax] and [zmin, zmax].
        if (ti < natoms) {
                d_atoms_in_x[ti] = ( (molecule_points[ti*3] > x_min) &&
                                (molecule_points[ti*3] < x_max) ) ? true : false;

                d_atoms_in_y[ti] = ( (molecule_points[ti*3 + 1] > y_min) &&
                                                (molecule_points[ti*3 + 1] < y_max) ) ? true : false;

                d_atoms_in_z[ti] = ( (molecule_points[ti*3 + 2] > z_min) &&
                                                (molecule_points[ti*3 + 2] < z_max) ) ? true : false;
        }
        __syncthreads();
        // Now get the joint set of the previos atoms to get those that lie inside
        // the cube
        if (ti < natoms) {
                atoms_in_bbox[ti] = d_atoms_in_x[ti] && d_atoms_in_y[ti] &&
                                d_atoms_in_z[ti];
        }

        // Now, index the "molecule_points" array to get those atoms
        return;
}
