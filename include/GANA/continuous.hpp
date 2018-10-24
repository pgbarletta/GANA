#ifndef GANA_CONTINUOUS_H
#define GANA_CONTINUOUS_H

extern float rsltion;
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
using EPIC = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = EPIC::Point_3;
using Polyhedron = CGAL::Polyhedron_3<EPIC>;

#include <cuda.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include "chemfiles.hpp"
#include <fmt/format.h>
#include "GANA/utils.hpp"

using P_Facet_iterator = Polyhedron::Facet_iterator;
using P_Facet_const_iterator = Polyhedron::Facet_const_iterator;
using P_Edge_iterator = Polyhedron::Edge_iterator;
using P_Edge_const_iterator = Polyhedron::Edge_const_iterator;
using P_Halfedge_around_facet_circulator = Polyhedron::Halfedge_around_facet_circulator;
using P_Halfedge_around_facet_const_circulator = Polyhedron::Halfedge_around_facet_const_circulator;
using P_Vertex_iterator = Polyhedron::Vertex_iterator;
using P_Vertex_const_iterator = Polyhedron::Vertex_const_iterator;


namespace GANA {

    class vector {
    public:
        vector() = default;
        vector(float x, float y, float z) : _vxyz{x, y, z},
            _origin{0., 0., 0.} {}
        vector(float x, float y, float z, float ox, float oy, float oz) :
            _vxyz{x, y, z},  _origin{ox, oy, oz} {}

        // Vector components access and modification.
        __host__ __device__ float &operator[] (int idx) {
        	return _vxyz[idx];
        }
        __host__ __device__ float operator[] (int idx) const {
        	return _vxyz[idx];
        }
        
        // Vector origin access and modification.
        float get_ox() const {
            return _origin[0];
        }
        float get_oy() const {
            return _origin[1];
        }
        float get_oz() const {
            return _origin[2];
        }
        void set_ox(float ox) {
            _origin[0] = ox;
        }
        void set_oy(float oy) {
            _origin[1] = oy;
        }
        void set_oz(float oz) {
            _origin[2] = oz;
        }
    private:
        float _vxyz[3], _origin[3];
    };
    
    inline std::ostream& operator<<(std::ostream &stream, const vector &v) {
	    stream << v[0] << " " << v[1] << " " << v[2];
	    return stream;
    }

    // Returns a vector starting on this same vector coordinates.
    inline vector operator+(const vector &lhs, const vector &rhs) {
        return vector(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2],
            lhs.get_ox(), lhs.get_oy(), lhs.get_oz());
    }

    // Returns a vector starting on this same vector coordinates.
    inline vector operator-(const vector &lhs, const vector &rhs) {
        return vector(lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2],
            lhs.get_ox(), lhs.get_oy(), lhs.get_oz());
    }

    inline bool operator==(const vector &lhs, const vector &rhs) {
        return (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2] &&
            lhs.get_ox() == rhs.get_ox() && lhs.get_oy() == rhs.get_oy() &&
            lhs.get_oz() == rhs.get_oz());
    }

    // Get the magnitude of the vector.
     inline float norm(const vector &v) {
    	 const float dx = v[0] - v.get_ox();
    	 const float dy = v[1] - v.get_oy();
    	 const float dz = v[2] - v.get_oz();
    	 return std::sqrt(dx*dx + dy*dy + dz*dz);
     }

    class point {
    public:
        point() = default;
        point(float x, float y, float z) : _xyz{x, y, z} {}
        point(Point p) : _xyz{static_cast<float>(CGAL::to_double(p.x())),
        	static_cast<float>(CGAL::to_double(p.y())),
        	static_cast<float>(CGAL::to_double(p.z()))} {}

        // Draw atom.
        void draw(FILE *out_file, unsigned int idx, unsigned int resid) {
		    fmt::print(
                out_file,
                "{: <6}{: >5} {: <4s} {:3} {:1}{: >4}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {: >2s}\n",
                "HETATM", idx, "H", "GPU", "A", resid, _xyz[0], _xyz[1], _xyz[2], 1.0, 0.0, "H"
			    );
		    return;
        }

        __host__ __device__ float &operator[] (int idx) {
        	return _xyz[idx];
        }
        __host__ __device__ float operator[] (int idx) const {
        	return _xyz[idx];
        }

    private:
        float _xyz[3];
    };
    
    inline std::ostream& operator<<(std::ostream &stream, const point &p) {
    	    stream << p[0] << " " << p[1] << " " << p[2];
    	    return stream;
    }

    // Returns a vector starting on this point coordinates.
    inline vector operator-(const point& lhs, const point& rhs) {
        return vector(lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2],
            lhs[0], lhs[1], lhs[2]);
    }
    // Displaces the point along the vector. 
    inline point operator+(const point& p, const vector& v) {
        return point(p[0] + (v[0] - v.get_ox()), p[1] + (v[1] - v.get_oy()),
            p[2] + (v[2] - v.get_oz()));
    }
    // Displaces the point along the vector. 
    inline point operator-(const point& p, const vector& v) {
        return point(p[0] - (v[0] - v.get_ox()), p[1] - (v[1] - v.get_oy()),
            p[2] - (v[2] - v.get_oz()));
    }
    
    inline bool operator==(const point& lhs, const point& rhs) {
        return (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2]);
    }

    inline vector point_to_vector(const point &in_point) {
        return vector(in_point[0], in_point[1], in_point[2]);
    }

    // Get the distance between 2 points
    inline float distance(const point &p0, const point &p1) {
   	 const float dx = p0[0] - p1[0];
   	 const float dy = p0[1] - p1[1];
   	 const float dz = p0[2] - p1[2];
   	 return std::sqrt(dx*dx + dy*dy + dz*dz);
    }

    class molecule {
    public:
        molecule() = default;
        molecule(const std::string &in_filename);

        ~molecule() {
            free(_xyz);
            free(_in_xyz);
            free(_radii);
            free(_in_radii);
        }

        void draw(const std::string &out_file);

        unsigned int _natoms;
        point *_xyz, *_in_xyz;
        float *_radii, *_in_radii;
    };

    class triangle {
    public:
        triangle() = default;
        triangle(float x0, float y0, float z0,
            float x1, float y1, float z1,
            float x2, float y2, float z2) : _p0(point(x0, y0, z0)),
            _p1(point(x1, y1, z1)), _p2(point(x2, y2, z2)) {}
        // From GANA::point
        triangle(point p0, point p1, point p2) : _p0(p0), _p1(p1), _p2(p2) {}
        // From CGAL Point
        triangle(Point p0, Point p1, Point p2) : _p0(point(p0)), _p1(point(p1)),
        _p2(point(p2)) {}
        // Draw triangle.
        void draw(FILE *out_file, unsigned int start_idx, unsigned int resid);

        point _p0, _p1, _p2;
    };
    
    inline std::ostream& operator<<(std::ostream &stream, const triangle& t) {
	    stream << t._p0 << "\t" << t._p1 << "\t" << t._p2;
	    return stream;
    }

    class convex_hull {
    public:
        convex_hull() = default;
            convex_hull(const molecule &prote,
            const std::vector<unsigned int> &indices);
        
        ~convex_hull() {
            free(_triangles);
        }
            
        void draw(const std::string &out_file);

        triangle *_triangles;
        unsigned int _ntriangles;
    };

    class tetrahedron {
    public:
        tetrahedron() = default;
        tetrahedron(float x0, float y0, float z0,
            float x1, float y1, float z1, float x2, float y2, float z2,
            float x3, float y3, float z3) : _p0(point(x0, y0, z0)),
            _p1(point(x1, y1, z1)), _p2(point(x2, y2, z2)),
            _p3(point(x3, y3, z3)) {}
        // From GANA::point
        tetrahedron(point p0, point p1, point p2, point p3) : 
            _p0(p0), _p1(p1), _p2(p2), _p3(p3) {}
        // From CGAL Point
        tetrahedron(Point p0, Point p1, Point p2, Point p3) : 
            _p0(point(p0)), _p1(point(p1)), _p2(point(p2)), _p3(point(p3)) {}
        // Draw tetrahedron.
        void draw (FILE *out_file, unsigned int start_idx, unsigned int resid);

        point _p0, _p1, _p2, _p3;
    };
    
    inline std::ostream& operator<<(std::ostream &stream, const tetrahedron& t) {
	    stream << t._p0 << "\t" << t._p1 << "\t" << t._p2 << "\t" << t._p3;
	    return stream;
    }

    class prism {
    public:
        prism() = default;
        prism(float x0, float y0, float z0,
            float x1, float y1, float z1, float x2, float y2, float z2,
            float x3, float y3, float z3, float x4, float y4, float z4,
            float x5, float y5, float z5, float x6, float y6, float z6,
            float x7, float y7, float z7) : _p0(point(x0, y0, z0)),
            _p1(point(x1, y1, z1)), _p2(point(x2, y2, z2)),
            _p3(point(x3, y3, z3)), _p4(point(x4, y4, z4)),
            _p5(point(x5, y5, z5)), _p6(point(x6, y6, z6)),
            _p7(point(x7, y7, z7)) {}
        // From GANA::point
        prism(point p0, point p1, point p2, point p3, point p4, point p5,
        point p6, point p7) : 
            _p0(p0), _p1(p1), _p2(p2), _p3(p3), _p4(p4), _p5(p5), _p6(p6),
            _p7(p7) {}
        // From CGAL Point
        prism(Point p0, Point p1, Point p2, Point p3, Point p4, Point p5,
            Point p6,  Point p7) :
            _p0(point(p0)), _p1(point(p1)), _p2(point(p2)), _p3(point(p3)),
            _p4(point(p4)), _p5(point(p5)), _p6(point(p6)), _p7(point(p7)) {}
        // Draw prism.
        void draw(FILE *out_file, unsigned int start_idx, unsigned int resid);
        
        point _p0, _p1, _p2, _p3, _p4, _p5, _p6, _p7;
    };
    std::ostream& operator<<(std::ostream &stream, const prism& t);

    class cube {
    public:
        cube() = default;
        cube(float p0x, float p0y, float p0z, float dim) :
            _p0(point(p0x, p0y, p0z)), _dim(dim) {
            _p1 = _p0 + vector(0.f, 0.f, _dim);
            _p2 = _p0 + vector(0.f, _dim, _dim);
            _p3 = _p0 + vector(0.f, _dim, 0.f);
            _p4 = _p0 + vector(_dim, 0.f, 0.f);
            _p5 = _p0 + vector(_dim, 0.f, _dim);
            _p6 = _p0 + vector(_dim, _dim, _dim);
            _p7 = _p0 + vector(_dim, _dim, 0.f);
        }
        // From GANA::point.
        cube(point p0, float dim) : _p0(p0), _dim(dim) {
            _p1 = _p0 + vector(0.f, 0.f, _dim);
            _p2 = _p0 + vector(0.f, _dim, _dim);
            _p3 = _p0 + vector(0.f, _dim, 0.f);
            _p4 = _p0 + vector(_dim, 0.f, 0.f);
            _p5 = _p0 + vector(_dim, 0.f, _dim);
            _p6 = _p0 + vector(_dim, _dim, _dim);
            _p7 = _p0 + vector(_dim, _dim, 0.f);
        }
        // From CGAL Point.
        cube(Point p0, float dim) : _p0(point(p0)), _dim(dim) {
            _p1 = _p0 + vector(0.f, 0.f, _dim);
            _p2 = _p0 + vector(0.f, _dim, _dim);
            _p3 = _p0 + vector(0.f, _dim, 0.f);
            _p4 = _p0 + vector(_dim, 0.f, 0.f);
            _p5 = _p0 + vector(_dim, 0.f, _dim);
            _p6 = _p0 + vector(_dim, _dim, _dim);
            _p7 = _p0 + vector(_dim, _dim, 0.f);
        }
        // Draw cube.
        void draw(FILE *out_file, unsigned int start_idx, unsigned int resid);

        point _p0, _p1, _p2, _p3, _p4, _p5, _p6, _p7;
        float _dim;
    };
    std::ostream& operator<<(std::ostream &stream, const cube& t);

    class triangulation {
    public:
        triangulation() = default;
        triangulation(const molecule &prote,
            const std::vector<unsigned int> &indices);
        
        ~triangulation() {
        	 free(_tetrahedrons);
         	 free(_bboxes);
         	 cudaFree(_Dtetrahedrons);
         	 cudaFree(_Dbboxes);
         }
            
        void draw(const std::string &out_file);

        unsigned int _ntetrahedrons;
        tetrahedron *_tetrahedrons, *_Dtetrahedrons;
        cube *_bboxes, *_Dbboxes;
    };

    class Bounding_box {
    public:
    	Bounding_box() = default;

    	Bounding_box(float xmin, float ymin, float zmin, float xmax,
    		float ymax, float zmax) noexcept : _xmin(xmin), _ymin(ymin),
    		_zmin(zmin), _xmax(xmax), _ymax(ymax), _zmax(zmax) {};

    	Bounding_box(tetrahedron const &in_tet) noexcept;

    	float _xmin, _ymin, _zmin, _xmax, _ymax, _zmax;
    };

    // Get the coordinates of the "indices" atoms as CGAL points.
inline void get_point_set(const molecule &prote,
	    const std::vector<unsigned int> &indices, std::vector<Point> &point_set) {

		point_set.reserve(indices.size());
		for(const auto &idx : indices) {
			const auto xyz = prote._xyz[idx - 1];
			point_set.push_back(Point(xyz[0], xyz[1], xyz[2]));
		}
		return;
	}
} // namespace 

#endif // _H
