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

    class Grid_point {
    public:
        using grid_idx_t = int32_t;

        Grid_point() = default;

        Grid_point(const grid_idx_t x, const grid_idx_t y, const grid_idx_t z)
            noexcept : _xyz{x, y, z} {}

        Grid_point(const point &p) noexcept :
            _xyz{cont_to_grid(p[0]), cont_to_grid(p[1]), cont_to_grid(p[2])} {}

        // Draw Grid_point as atom.
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
    std::ostream& operator<<(std::ostream &stream, const Grid_point& t);
    // Turn a grid point into a continuous point, given the resolution.
    point Grid_point_to_point(const Grid_point &in_point);

    // Turn a grid point into a continuous point, given the resolution.
    Grid_point point_to_Grid_point(const point &in_point);

    class Grid_molecule {
    public:
        Grid_molecule() = default;
        Grid_molecule(const molecule &in_filename, const point &orig_point);

        ~Grid_molecule() {
            free(_xyz);
            free(_in_xyz);
            free(_radii);
            free(_in_radii);
        }

        void draw(const std::string &ou_fil);

        unsigned int _natoms;
        vector _orig_vtor;
        Grid_point *_xyz, *_in_xyz;
        float *_radii, *_in_radii;
    };

    class Grid_convex_hull {
    public:
        Grid_convex_hull() = default;
        Grid_convex_hull(const convex_hull& CH);

        ~Grid_convex_hull() {
            free(dots);
        }

        void draw(const std::string &ou_fil);
        
        float *dots;
    };

    class Grid_bool {
    public:
        Grid_bool() = default;
        Grid_bool(const convex_hull& CH);

    };
} //namespace GANA

#endif // GANA_GRID
