#include "GANA/primitives.hpp"
#include <cuda.h>
#include <cuda_runtime.h>

namespace GANA {
// Draw triangle.
void Triangle::draw(FILE *out_file, unsigned int start_idx, unsigned int resid) {
	const auto i = start_idx++;
	const auto j = start_idx++;
	const auto k = start_idx;

	_p0.draw(out_file, i, resid);
	_p1.draw(out_file, j, resid);
    _p2.draw(out_file, k, resid);

	fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", i, j, k);
	fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", j, k, i);
	return;
}

// Draw tetrahedron.
void Tetrahedron::draw (FILE *out_file, unsigned int start_idx, unsigned int resid) {
	const auto i = start_idx++;
	const auto j = start_idx++;
	const auto k = start_idx++;
	const auto l = start_idx;

	_p0.draw(out_file, i, resid);
	_p1.draw(out_file, j, resid);
    _p2.draw(out_file, k, resid);
	_p3.draw(out_file, l, resid);

	fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", i, j, k, l);
	fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", j, k, l, i);
	fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", k, l, i, j);
	fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", l, i, j, k);
	return;
}

// Draw cube.
void Cube::draw(FILE *out_file, unsigned int start_idx, unsigned int resid) {
	const auto i = start_idx++;
	const auto j = start_idx++;
	const auto k = start_idx++;
	const auto l = start_idx++;
	const auto ii = start_idx++;
	const auto jj = start_idx++;
	const auto kk = start_idx++;
	const auto ll = start_idx;

	_p0.draw(out_file, i, resid);
	_p1.draw(out_file, j, resid);
	_p2.draw(out_file, k, resid);
	_p3.draw(out_file, l, resid);
	_p4.draw(out_file, ii, resid);
	_p5.draw(out_file, jj, resid);
	_p6.draw(out_file, kk, resid);
	_p7.draw(out_file, ll, resid);

	fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", i, j, l, ii);
	fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", k, j, l, kk);
	fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", jj, ii, kk, j);
	fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", ll, ii, kk, l);

	return;
}

// Draw prism. Can't draw connectivity properly if the prism wasn't constructed
// with proper Point ordering. SO this class is kind of useless.
void Prism::draw(FILE *out_file, unsigned int start_idx, unsigned int resid) {
	const auto i = start_idx++;
	const auto j = start_idx++;
	const auto k = start_idx++;
	const auto l = start_idx++;
	const auto ii = start_idx++;
	const auto jj = start_idx++;
	const auto kk = start_idx++;
	const auto ll = start_idx;

	_p0.draw(out_file, i, resid);
	_p1.draw(out_file, j, resid);
	_p2.draw(out_file, k, resid);
	_p3.draw(out_file, l, resid);
	_p4.draw(out_file, ii, resid);
	_p5.draw(out_file, jj, resid);
	_p6.draw(out_file, kk, resid);
	_p7.draw(out_file, ll, resid);

	fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", i, j, l, ii);
	fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", k, j, l, kk);
	fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", jj, ii, kk, j);
	fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", ll, ii, kk, l);

	return;
}

}
