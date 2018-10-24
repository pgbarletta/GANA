#include "GANA/grid.hpp"
namespace GANA {

// Draw Grid_point as atom.
void Grid_point::draw(FILE *ou_fil, unsigned int idx, unsigned int resid) {
    const float fx = grid_to_cont(_xyz[0]);
    const float fy = grid_to_cont(_xyz[1]);
    const float fz = grid_to_cont(_xyz[2]);
    fmt::print(
        ou_fil,
        "{: <6}{: >5} {: <4s} {:3} {:1}{: >4}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {: >2s}\n",
        "HETATM", idx, "H", "GPU", "A", resid, fx, fy, fz, 1.0, 0.0, "H");
    return;
}

std::ostream& operator<<(std::ostream &stream, const Grid_point &p) {
	stream << p[0] << " " << p[1] << " " << p[2];
	return stream;
}

// Turn a grid point into a continuous point, given the resolution.
point Grid_point_to_point(const Grid_point &in_point) {
	return point(grid_to_cont(in_point[0]), grid_to_cont(in_point[1]),
		grid_to_cont(in_point[2]));
}

// Turn a grid point into a continuous point, given the resolution.
Grid_point point_to_Grid_point(const point &in_point) {
	return Grid_point(grid_to_cont(in_point[0]), grid_to_cont(in_point[1]),
		grid_to_cont(in_point[2]));
}

Grid_molecule::Grid_molecule(const molecule &in_mol,
    const point &orig_point) {

    _orig_vtor = point_to_vector(orig_point);
    _natoms = in_mol._natoms;

    _xyz = (Grid_point*) malloc(sizeof(Grid_point) * _natoms * 3);
	_in_xyz = (Grid_point*) malloc(sizeof(Grid_point) * _natoms * 3);
	_radii = (float*) malloc(sizeof(float) * _natoms);
	_in_radii = (float*) malloc(sizeof(float) * _natoms);

    for (size_t i = 0 ; i < _natoms ; ++i) {
        _xyz[i] = point_to_Grid_point(in_mol._xyz[i] - _orig_vtor);
        _radii[i] = in_mol._radii[i];
	}
}

void Grid_molecule::draw(const std::string &ou_fil) {
	
	FILE *file = std::fopen(ou_fil.c_str(), "w");
	if(file) {
		for (size_t i = 1 ; i <= _natoms ; ++i) {
			_xyz[i-1].draw(file, i, i);
		}
	} else {
		std::cerr << "Could not open " << ou_fil << ". " << '\n';
	}
	std::fclose(file);
return;
}

}
