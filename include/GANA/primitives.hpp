#ifndef GANA_PRIMITIVES_H
#define GANA_PRIMITIVES_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
using EPIC = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = EPIC::Point_3;

#include <cuda.h>
#include <fmt/format.h>

namespace GANA {

	class Vector {
    public:
        Vector() = default;
        Vector(float x, float y, float z) noexcept : _vxyz{x, y, z},
            _origin{0., 0., 0.} {}
        Vector(float x, float y, float z, float ox, float oy, float oz) noexcept :
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

    inline std::ostream& operator<<(std::ostream &stream, const Vector &v) {
	    stream << v[0] << " " << v[1] << " " << v[2];
	    return stream;
    }

    // Returns a Vector starting on this same Vector coordinates.
    inline Vector operator+(const Vector &lhs, const Vector &rhs) {
        return Vector(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2],
            lhs.get_ox(), lhs.get_oy(), lhs.get_oz());
    }

    // Returns a Vector starting on this same Vector coordinates.
    inline Vector operator-(const Vector &lhs, const Vector &rhs) {
        return Vector(lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2],
            lhs.get_ox(), lhs.get_oy(), lhs.get_oz());
    }

    inline bool operator==(const Vector &lhs, const Vector &rhs) {
        return (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2] &&
            lhs.get_ox() == rhs.get_ox() && lhs.get_oy() == rhs.get_oy() &&
            lhs.get_oz() == rhs.get_oz());
    }

    // Get the magnitude of the Vector.
     inline float norm(const Vector &v) {
    	 const float dx = v[0] - v.get_ox();
    	 const float dy = v[1] - v.get_oy();
    	 const float dz = v[2] - v.get_oz();
    	 return std::sqrt(dx*dx + dy*dy + dz*dz);
     }

    class Point {
    public:
        Point() = default;
        Point(float x, float y, float z) : _xyz{x, y, z} {}
        Point(Point_3 p) : _xyz{static_cast<float>(CGAL::to_double(p.x())),
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

    inline std::ostream& operator<<(std::ostream &stream, const Point &p) {
    	    stream << p[0] << " " << p[1] << " " << p[2];
    	    return stream;
    }

    // Returns a Vector starting on this Point coordinates.
    inline Vector operator-(const Point& lhs, const Point& rhs) {
        return Vector(lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2],
            lhs[0], lhs[1], lhs[2]);
    }
    // Displaces the Point along the Vector.
    inline Point operator+(const Point& p, const Vector& v) {
        return Point(p[0] + (v[0] - v.get_ox()), p[1] + (v[1] - v.get_oy()),
            p[2] + (v[2] - v.get_oz()));
    }
    // Displaces the Point along the Vector.
    inline Point operator-(const Point& p, const Vector& v) {
        return Point(p[0] - (v[0] - v.get_ox()), p[1] - (v[1] - v.get_oy()),
            p[2] - (v[2] - v.get_oz()));
    }

    inline bool operator==(const Point& lhs, const Point& rhs) {
        return (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2]);
    }

    inline Vector point_to_vector(const Point &in_point) {
        return Vector(in_point[0], in_point[1], in_point[2]);
    }

    // Get the distance between 2 points
    inline float distance(const Point &p0, const Point &p1) {
   	 const float dx = p0[0] - p1[0];
   	 const float dy = p0[1] - p1[1];
   	 const float dz = p0[2] - p1[2];
   	 return std::sqrt(dx*dx + dy*dy + dz*dz);
    }

    class Triangle {
    public:
        Triangle() = default;
        Triangle(float x0, float y0, float z0,
            float x1, float y1, float z1,
            float x2, float y2, float z2) : _p0(Point(x0, y0, z0)),
            _p1(Point(x1, y1, z1)), _p2(Point(x2, y2, z2)) {}
        // From GANA::Point
        Triangle(Point p0, Point p1, Point p2) : _p0(p0), _p1(p1), _p2(p2) {}
        // From CGAL Point
        Triangle(Point_3 p0, Point_3 p1, Point_3 p2) : _p0(Point(p0)), _p1(Point(p1)),
        _p2(Point(p2)) {}
        // Draw triangle.
        void draw(FILE *out_file, unsigned int start_idx, unsigned int resid);

        Point _p0, _p1, _p2;
    };

    inline std::ostream& operator<<(std::ostream &stream, const Triangle& t) {
	    stream << t._p0 << "\t" << t._p1 << "\t" << t._p2;
	    return stream;
    }

    class Tetrahedron {
    public:
        Tetrahedron() = default;
        Tetrahedron(float x0, float y0, float z0,
            float x1, float y1, float z1, float x2, float y2, float z2,
            float x3, float y3, float z3) : _p0(Point(x0, y0, z0)),
            _p1(Point(x1, y1, z1)), _p2(Point(x2, y2, z2)),
            _p3(Point(x3, y3, z3)) {}
        // From GANA::Point
        Tetrahedron(Point p0, Point p1, Point p2, Point p3) :
            _p0(p0), _p1(p1), _p2(p2), _p3(p3) {}
        // From CGAL Point
        Tetrahedron(Point_3 p0, Point_3 p1, Point_3 p2, Point_3 p3) :
            _p0(Point(p0)), _p1(Point(p1)), _p2(Point(p2)), _p3(Point(p3)) {}
        // Draw tetrahedron.
        void draw (FILE *out_file, unsigned int start_idx, unsigned int resid);

        Point _p0, _p1, _p2, _p3;
    };

    inline std::ostream& operator<<(std::ostream &stream, const Tetrahedron& t) {
	    stream << t._p0 << "\t" << t._p1 << "\t" << t._p2 << "\t" << t._p3;
	    return stream;
    }

    class Prism {
    public:
        Prism() = default;
        Prism(float x0, float y0, float z0,
            float x1, float y1, float z1, float x2, float y2, float z2,
            float x3, float y3, float z3, float x4, float y4, float z4,
            float x5, float y5, float z5, float x6, float y6, float z6,
            float x7, float y7, float z7) : _p0(Point(x0, y0, z0)),
            _p1(Point(x1, y1, z1)), _p2(Point(x2, y2, z2)),
            _p3(Point(x3, y3, z3)), _p4(Point(x4, y4, z4)),
            _p5(Point(x5, y5, z5)), _p6(Point(x6, y6, z6)),
            _p7(Point(x7, y7, z7)) {}
        // From GANA::Point
        Prism(Point p0, Point p1, Point p2, Point p3, Point p4, Point p5,
        Point p6, Point p7) :
            _p0(p0), _p1(p1), _p2(p2), _p3(p3), _p4(p4), _p5(p5), _p6(p6),
            _p7(p7) {}
        // From CGAL Point
        Prism(Point_3 p0, Point_3 p1, Point_3 p2, Point_3 p3, Point_3 p4,
        	Point_3 p5, Point_3 p6,  Point_3 p7) :
            _p0(Point(p0)), _p1(Point(p1)), _p2(Point(p2)), _p3(Point(p3)),
            _p4(Point(p4)), _p5(Point(p5)), _p6(Point(p6)), _p7(Point(p7)) {}
        // Draw prism.
        void draw(FILE *out_file, unsigned int start_idx, unsigned int resid);

        Point _p0, _p1, _p2, _p3, _p4, _p5, _p6, _p7;
    };
    std::ostream& operator<<(std::ostream &stream, const Prism& t);

    class Cube {
    public:
        Cube() = default;
        Cube(float p0x, float p0y, float p0z, float dim) :
            _p0(Point(p0x, p0y, p0z)), _dim(dim) {
            _p1 = _p0 + Vector(0.f, 0.f, _dim);
            _p2 = _p0 + Vector(0.f, _dim, _dim);
            _p3 = _p0 + Vector(0.f, _dim, 0.f);
            _p4 = _p0 + Vector(_dim, 0.f, 0.f);
            _p5 = _p0 + Vector(_dim, 0.f, _dim);
            _p6 = _p0 + Vector(_dim, _dim, _dim);
            _p7 = _p0 + Vector(_dim, _dim, 0.f);
        }
        // From GANA::Point.
        Cube(Point p0, float dim) : _p0(p0), _dim(dim) {
            _p1 = _p0 + Vector(0.f, 0.f, _dim);
            _p2 = _p0 + Vector(0.f, _dim, _dim);
            _p3 = _p0 + Vector(0.f, _dim, 0.f);
            _p4 = _p0 + Vector(_dim, 0.f, 0.f);
            _p5 = _p0 + Vector(_dim, 0.f, _dim);
            _p6 = _p0 + Vector(_dim, _dim, _dim);
            _p7 = _p0 + Vector(_dim, _dim, 0.f);
        }
        // From CGAL Point.
        Cube(Point_3 p0, float dim) : _p0(Point(p0)), _dim(dim) {
            _p1 = _p0 + Vector(0.f, 0.f, _dim);
            _p2 = _p0 + Vector(0.f, _dim, _dim);
            _p3 = _p0 + Vector(0.f, _dim, 0.f);
            _p4 = _p0 + Vector(_dim, 0.f, 0.f);
            _p5 = _p0 + Vector(_dim, 0.f, _dim);
            _p6 = _p0 + Vector(_dim, _dim, _dim);
            _p7 = _p0 + Vector(_dim, _dim, 0.f);
        }
        // Draw cube.
        void draw(FILE *out_file, unsigned int start_idx, unsigned int resid);

        Point _p0, _p1, _p2, _p3, _p4, _p5, _p6, _p7;
        float _dim;
    };
    std::ostream& operator<<(std::ostream &stream, const Cube& t);

}
#endif // _H
