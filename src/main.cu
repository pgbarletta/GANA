float rsltion = .1;
#include <cuda.h>
#include <cuda_runtime.h>

#include <vector>
#include <stdexcept>
#include <string>
#include <typeinfo>

#include "GANA/utils.hpp"
#include "GANA/continuous.hpp"
#include "GANA/grid.hpp"
#include "GANA/kernels.cu"

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}
#define GEC(ans) { gpuAssert((ans), __FILE__, __LINE__); } // GPU Error Check

int main(int argc, char **argv) {
	// Get positions and associated variables ready.
	if (argc != 4) {
		std::cerr << "Usage: GANA in_pdb resolution ou_pdb" << '\n';
		return 0;
	}

	std::vector<unsigned int> indices = {300, 600, 900, 1200, 1500, 1800, 1240,
		400, 500, 700, 800, 1000, 1100};
	///////////

	try {
		rsltion = std::stof(argv[2]);
	} catch(...) {
		std::cerr << "Bad input resolution. Please specify a "
			<< "number between .01 and 1" << '\n';
	}

	/////////////////////////
	GANA::molecule prote(argv[1]), *Dprote;
	// Paso molécula a GPU. Falta pasar los arrays.
	GEC( cudaMalloc((void **) &Dprote, sizeof(GANA::molecule)) );
	GEC( cudaMemcpyAsync(Dprote, &prote, sizeof(GANA::molecule), cudaMemcpyHostToDevice) );

	// Paso los arrays.
	float *Dradii, *Din_radii;
	GANA::point *Dxyz, *Din_xyz;
	const auto xyz_sz = sizeof(float) * prote._natoms * 3,
		rad_sz = sizeof(float) * prote._natoms;
	GEC( cudaMalloc((void **) &Dxyz, xyz_sz) );
	GEC( cudaMalloc((void **) &Din_xyz, xyz_sz) );
	GEC( cudaMalloc((void **) &Dradii, rad_sz) );
	GEC( cudaMalloc((void **) &Din_radii, rad_sz) );
	GEC( cudaMemcpyAsync(Dxyz, prote._xyz, xyz_sz, cudaMemcpyHostToDevice) );
	GEC( cudaMemcpyAsync(Din_xyz, prote._in_xyz, xyz_sz, cudaMemcpyHostToDevice) );
	GEC( cudaMemcpyAsync(Dradii, prote._radii, rad_sz, cudaMemcpyHostToDevice) );
	GEC( cudaMemcpyAsync(Din_radii, prote._in_radii, rad_sz, cudaMemcpyHostToDevice) );

	// Apunto los pointers de la molécula (en GPU) a los arrays (en GPU).
	GEC( cudaMemcpyAsync(&(Dprote->_xyz), &Dxyz, sizeof(GANA::point *),
		cudaMemcpyHostToDevice) );
	GEC( cudaMemcpyAsync(&(Dprote->_in_xyz), &Din_xyz, sizeof(GANA::point *),
		cudaMemcpyHostToDevice) );
	GEC( cudaMemcpyAsync(&(Dprote->_radii), &Dradii, sizeof(float *),
		cudaMemcpyHostToDevice) );
	GEC( cudaMemcpyAsync(&(Dprote->_in_radii), &Din_radii, sizeof(float *),
		cudaMemcpyHostToDevice) );
	/////////////////////////

	GANA::triangulation incl_area(prote, indices);
	incl_area.draw("aux/ia.pdb");

	const dim3 dB0(1024, 1, 1);
	const dim3 dG0(10, 1, 1);
	empiezo<<<dG0, dB0>>>(incl_area._Dtetrahedrons, incl_area._ntetrahedrons);
	cudaDeviceSynchronize();
	cudaMemcpy(incl_area._tetrahedrons, incl_area._Dtetrahedrons,
			sizeof(GANA::tetrahedron) * incl_area._ntetrahedrons, cudaMemcpyDeviceToHost);

	incl_area.draw("aux/after.pdb");


	return 0;
}
