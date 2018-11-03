#include "GANA/continuous.hpp"

#include <CGAL/Delaunay_triangulation_3.h>
using Delaunay = CGAL::Delaunay_triangulation_3<EPIC>;
using Finite_cells_iterator = Delaunay::Finite_cells_iterator;
#include <cuda.h>
#include <cuda_runtime.h>

namespace GANA {

Molecule::Molecule(std::string const &in_filename) {
	chemfiles::Trajectory in_trj(in_filename);
	auto in_frm = in_trj.read();
	auto in_top = in_frm.topology();
	auto in_xyz = in_frm.positions();
	_natoms = in_xyz.size();

	_xyz = (Point*) malloc(sizeof(Point) * _natoms * 3);
	_in_xyz = (Point*) malloc(sizeof(Point) * _natoms * 3);
	_radii = (float*) malloc(sizeof(float) * _natoms);
	_in_radii = (float*) malloc(sizeof(float) * _natoms);

	// Get atoms positions and VdW radii.
	size_t j = 0;
	for (const auto &residuo : in_top.residues()) {
		for (const auto &i : residuo) {
			const auto atom = in_xyz[i];
			_xyz[i] = Point(atom[0], atom[1], atom[2]);
			_radii[j] = in_top[i].vdw_radius().value_or(1.5);
			++j;
		}
	}
}

// Draw the molecule in the **out_file** path in PDB format.
void Molecule::draw(std::string const &out_file) {
	
	FILE *file = std::fopen(out_file.c_str(), "w");
	if(file) {
		for (size_t i = 1 ; i <= _natoms ; ++i) {
			_xyz[i-1].draw(file, i, i);
		}
	} else {
		std::cerr << "Could not open " << out_file << ". " << '\n';
	}
	std::fclose(file);
	return;
}

ConvexHull::ConvexHull(
	Molecule const &prote, std::vector<int> const &indices) {

	// Get the coordinates of the input indices atoms.
	auto point_vtor = atom_indices_to_points(prote, indices);
	
	// Get the convex hull.
	Polyhedron con_hul;
	CGAL::convex_hull_3(point_vtor.begin(), point_vtor.end(), con_hul);

	// Turn CGAL's polyhedron holding the convex hull into GANA triangles.
	P_Facet_const_iterator f_ite = con_hul.facets_begin();
	const P_Facet_const_iterator f_end = con_hul.facets_end();
	_ntriangles = std::distance(f_ite, f_end);
	_triangles = (Triangle *) malloc(sizeof(Triangle) * _ntriangles);

	for (size_t i = 0 ; f_ite != f_end ; ++f_ite) {
	    P_Halfedge_around_facet_const_circulator he_ite = f_ite->facet_begin();
		_triangles[i] = Triangle(he_ite->vertex()->point(),
		(he_ite++)->vertex()->point(), (he_ite++)->vertex()->point());
		++i;
	}
}

void ConvexHull::draw(std::string const &out_file) {
	
	FILE *file = std::fopen(out_file.c_str(), "w");
	if(file) {
		int resid = 1;
		for (size_t i = 0 ; i < _ntriangles ; ++i) {
			const auto start_idx = i * 3 + 1;
			_triangles[i].draw(file, start_idx, resid++);
		}
		for (size_t i = 0 ; i < _ntriangles ; ++i) {
			const auto j = i * 3 + 1;
			fmt::print(file, "CONECT {:>4} {:>4} {:>4}\n", j, j+1, j+2);
			fmt::print(file, "CONECT {:>4} {:>4} {:>4}\n", j+1, j, j+2);
		}
	} else {
		std::cerr << "Could not open " << out_file << ". " << '\n';
	}
	std::fclose(file);
	
	return;
}

Triangulation::Triangulation(
	Molecule const &prote, std::vector<int> const &indices) {

	// Get the coordinates of the input indices atoms.
	auto point_vtor = atom_indices_to_points(prote, indices);
	
	// Get the convex hull.
	Polyhedron con_hul;
	CGAL::convex_hull_3(point_vtor.begin(), point_vtor.end(), con_hul);
	
	// Get the triangulation from the convex hull.
	Delaunay T(con_hul.points_begin(), con_hul.points_end());
	
	// Turn CGAL's Delaunay triangulation into GANA triangulation.
	auto cell_ite = T.finite_cells_begin();
	const auto cell_end = T.finite_cells_end();

	_ntetrahedrons = std::distance(cell_ite, cell_end);
	const int sz_t = sizeof(Tetrahedron) * _ntetrahedrons;
	const int sz_c = sizeof(Cube) * _ntetrahedrons;

	_tetrahedrons = (Tetrahedron *) malloc(sz_t);
	cudaMalloc((void **)&_Dtetrahedrons, sz_t);
	_bboxes = (Cube *) malloc(sz_c);
	cudaMalloc((void **)&_Dbboxes, sz_c);

	// // Initialize bounding box.
	// constexpr float x = -999.0f;
	// constexpr float y = -999.0f;
	// constexpr float z = -999.0f;
	// bbo_x = cube(point(-x, -y, -z), point(-x, -y, z), point(-x, y, z),
	// 	point(-x, y, -z), point(x, -y, -z), point(x, -y, z),
	// 	point(x, y, z), point(x, y, -z));

	for (size_t i = 0 ; i < _ntetrahedrons ; ++i) {
		const auto p0 = cell_ite->vertex(0)->point();
		const auto p1 = cell_ite->vertex(1)->point();
		const auto p2 = cell_ite->vertex(2)->point();
		const auto p3 = cell_ite->vertex(3)->point();
		cell_ite++;
        _tetrahedrons[i] = Tetrahedron(p0, p1, p2, p3);
	}
	cudaMemcpyAsync(_Dtetrahedrons, _tetrahedrons, sz_t, cudaMemcpyHostToDevice);

	return;
}

void Triangulation::draw(std::string const &out_file) {
	
	FILE *file = std::fopen(out_file.c_str(), "w");
	if(file) {
		int resid = 1;
		for (size_t i = 0 ; i < _ntetrahedrons ; ++i) {
			const auto start_idx = i * 4 + 1;
			_tetrahedrons[i].draw(file, start_idx, resid++);
		}
		for (size_t i = 0 ; i < _ntetrahedrons ; ++i) {
			const auto j = i * 4 + 1;
			fmt::print(file, "CONECT {:>4} {:>4} {:>4}\n",
				j, j+1, j+2, j+3);
			fmt::print(file, "CONECT {:>4} {:>4} {:>4}\n",
				j+1, j+2, j+3, j);
			fmt::print(file, "CONECT {:>4} {:>4} {:>4}\n",
				j+2, j+3, j, j+1);
			fmt::print(file, "CONECT {:>4} {:>4} {:>4}\n",
				j+3, j, j+1, j+2);
		}
	} else {
		std::cerr << "Could not open " << out_file << ". " << '\n';
	}
	std::fclose(file);
	
	return;
}

BoundingBox::BoundingBox(Tetrahedron const &in_T) noexcept {
	 //in_T._p0.

	return;
}



} // namespace GANA
