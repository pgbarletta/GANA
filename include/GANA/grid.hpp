#ifndef GANA_GRID
#define GANA_GRID

extern float rsltion;
#include <string>
#include <fstream>
#include <iostream>
#include <cuda.h>

#include <fmt/format.h>
#include "GANA/continuous.hpp"
#include "GANA/utils.hpp"

namespace GANA {

    class GridPoint {
    public:
        using grid_idx_t = int32_t;

        GridPoint() = default;

        GridPoint(const grid_idx_t x, const grid_idx_t y, const grid_idx_t z)
            noexcept : _xyz{x, y, z} {}

        GridPoint(const Point &p) noexcept :
            _xyz{cont_to_grid(p[0]), cont_to_grid(p[1]), cont_to_grid(p[2])} {}

        // Draw GridPoint as atom.
        void draw(FILE *ou_fil, unsigned int idx, unsigned int resid);

        // // Returns a vector starting on this point coordinates.
        // vector operator-(const point& p) const {
        //     return vector(_x - p[0], _y - p[1], _z - p[2],
        //         _x, _y, _z);
        // }
        // bool operator==(const point &p) const {
        //         return (_x == p[0] && _y == p[1] && _z == p[2]);
        // }
        grid_idx_t &operator[] (int idx) {
        	return _xyz[idx];
        }
        grid_idx_t operator[] (int idx) const {
        	return _xyz[idx];
        }
    private:
        grid_idx_t _xyz[3];
    };
    std::ostream& operator<<(std::ostream &stream, const GridPoint& t);
    // Turn a grid point into a continuous point, given the resolution.
    Point GridPoint_to_point(const GridPoint &in_point);

    // Turn a grid point into a continuous point, given the resolution.
    GridPoint point_to_GridPoint(const Point &in_point);

    class GridMolecule {
    public:
        GridMolecule() = default;
        GridMolecule(const Molecule &in_mol, const Point &orig_point);

        ~GridMolecule() {
            free(_xyz);
            free(_in_xyz);
            free(_radii);
            free(_in_radii);
        }

        void draw(const std::string &ou_fil);

        unsigned int _natoms;
        Vector _orig_vtor;
        GridPoint *_xyz, *_in_xyz;
        float *_radii, *_in_radii;
    };

    class GridConvexHull {
    public:
        GridConvexHull() = default;
        GridConvexHull(const ConvexHull& CH);

        ~GridConvexHull() {
            free(dots);
        }

        void draw(const std::string &ou_fil);
        
        float *dots;
    };

    class GridBool {
    public:
        GridBool() = default;
        GridBool(const ConvexHull& CH);

    };
} //namespace GANA

#endif // GANA_GRID
