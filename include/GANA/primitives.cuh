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

    Vector(float const x, float const y, float const z) noexcept :
    	_vxyz{x, y, z}, _origin{0., 0., 0.} {}

    Vector(float const x, float const y, float const z,
    	float const ox, float const oy, float const oz) noexcept :
        _vxyz{x, y, z},  _origin{ox, oy, oz} {}

    // Vector components access and modification.
    __host__ __device__
    float &operator[] (int const idx) {
    	return _vxyz[idx];
    }
    __host__ __device__
    float operator[] (int const idx) const {
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

inline std::ostream& operator<<(std::ostream &stream, Vector const &v) {
    stream << v[0] << " " << v[1] << " " << v[2];
    return stream;
}

// Returns a Vector starting on this same Vector coordinates.
inline Vector operator+(Vector const &lhs, Vector const &rhs) {
    return Vector(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2],
        lhs.get_ox(), lhs.get_oy(), lhs.get_oz());
}

// Returns a Vector starting on this same Vector coordinates.
inline Vector operator-(Vector const &lhs, Vector const &rhs) {
    return Vector(lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2],
        lhs.get_ox(), lhs.get_oy(), lhs.get_oz());
}

inline bool operator==(Vector const &lhs, Vector const &rhs) {
    return (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2] &&
        lhs.get_ox() == rhs.get_ox() && lhs.get_oy() == rhs.get_oy() &&
        lhs.get_oz() == rhs.get_oz());
}

// Get the magnitude of the Vector.
 inline float norm(Vector const &v) {
	 const float dx = v[0] - v.get_ox();
	 const float dy = v[1] - v.get_oy();
	 const float dz = v[2] - v.get_oz();
	 return std::sqrt(dx*dx + dy*dy + dz*dz);
 }

class Point {
public:
    Point() = default;

    Point(float const x, float const y, float const z) : _xyz{x, y, z} {}

    Point(Point_3 const p) : _xyz{static_cast<float>(CGAL::to_double(p.x())),
    	static_cast<float>(CGAL::to_double(p.y())),
    	static_cast<float>(CGAL::to_double(p.z()))} {}

    // Draw atom.
    void draw(FILE *out_file, int const idx, int const resid) {
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

inline std::ostream& operator<<(std::ostream &stream, Point const &p) {
	    stream << p[0] << " " << p[1] << " " << p[2];
	    return stream;
}

// Returns a Vector starting on this Point coordinates.
inline Vector operator-(Point const &lhs, Point const &rhs) {
    return Vector(lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2],
        lhs[0], lhs[1], lhs[2]);
}
// Displaces the Point along the Vector.
inline Point operator+(Point const &p, Vector const &v) {
    return Point(p[0] + (v[0] - v.get_ox()), p[1] + (v[1] - v.get_oy()),
        p[2] + (v[2] - v.get_oz()));
}
// Displaces the Point along the Vector.
inline Point operator-(Point const &p, Vector const &v) {
    return Point(p[0] - (v[0] - v.get_ox()), p[1] - (v[1] - v.get_oy()),
        p[2] - (v[2] - v.get_oz()));
}

inline bool operator==(Point const &lhs, Point const &rhs) {
    return (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2]);
}

inline Vector point_to_vector(Point const &in_point) {
    return Vector(in_point[0], in_point[1], in_point[2]);
}

// Get the distance between 2 points
inline float distance(Point const &p0, Point const &p1) {
 float const dx = p0[0] - p1[0];
 float const dy = p0[1] - p1[1];
 float const dz = p0[2] - p1[2];
 return std::sqrt(dx*dx + dy*dy + dz*dz);
}

class Triangle {
public:
    Triangle() = default;

    Triangle(float const x0, float const y0, float const z0,
        float const x1, float const y1, float const z1,
        float const x2, float const y2, float const z2) :
        _p({Point(x0, y0, z0), Point(x1, y1, z1), Point(x2, y2, z2)}) {}

    // From GANA::Point
    Triangle(Point const &p0, Point const &p1, Point const &p2) :
    	_p({p0, p1, p2}) {}

    // From CGAL Point
    Triangle(Point_3 const &p0, Point_3 const &p1, Point_3 const &p2) :
    	_p({Point(p0), Point(p1), Point(p2)}) {}

    // Draw triangle.
    void draw(FILE *out_file, int const start_idx, int const resid);

    __host__ __device__
    Point operator[](int const idx) const {
    	return _p[idx];
    }
    __host__ __device__
    Point &operator[](int const idx) {
    	return _p[idx];
    }

    Point _p[3];
};

inline std::ostream& operator<<(std::ostream &stream, const Triangle& t) {
    stream << t._p[0] << "\t" << t._p[1] << "\t" << t._p[2];
    return stream;
}

class Tetrahedron {
public:
    Tetrahedron() = default;
    Tetrahedron(float const x0, float const y0, float const z0,
        float const x1, float const y1, float const z1,
        float const x2, float const y2, float const z2,
        float const x3, float const y3, float const z3) :
        _p({Point(x0, y0, z0), Point(x1, y1, z1), Point(x2, y2, z2),
    	Point(x3, y3, z3)}) {}

    // From GANA::Point
    Tetrahedron(Point const &p0, Point const &p1, Point const &p2,
    		Point const &p3) : _p({p0, p1, p2, p3}) {}

    // From CGAL Point
    Tetrahedron(Point_3 const &p0, Point_3 const &p1, Point_3 const &p2,
    	Point_3 const &p3) : _p({Point(p0), Point(p1), Point(p2), Point(p3)}) {}

    // Draw tetrahedron.
    void draw (FILE *out_file, int const start_idx, int const resid);

    __host__ __device__
    Point operator[](int const idx) const {
    	return _p[idx];
    }
    __host__ __device__
    Point &operator[](int const idx) {
    	return _p[idx];
    }

    Point _p[4];
};

inline std::ostream& operator<<(std::ostream &stream, Tetrahedron const &t) {
    stream << t._p[0] << "\t" << t._p[1] << "\t" << t._p[2] << "\t" << t._p[3];
    return stream;
}

class Cube {
public:
    Cube() = default;

    Cube(float const p0x, float const p0y, float const p0z, float const dim);

    // From GANA::Point.
    Cube(Point const p0, float const dim);

    // From CGAL Point.
    Cube(Point_3 const p0, float const dim);

    // Draw cube.
    void draw(FILE *out_file, int const start_idx, int const resid);

    Point _p[8];
    float _dim;
};
std::ostream& operator<<(std::ostream &stream, Cube const &t);

class Prism {
public:
    Prism() = default;

    // From GANA::Point
    Prism(Point const &p0, Point const &p1, Point const &p2, Point const &p3,
    	Point const &p4, Point const &p5, Point const &p6, Point const &p7);

    // From CGAL Point
    Prism(Point_3 const &p0, Point_3 const &p1, Point_3 const &p2,
    	Point_3 const &p3, Point_3 const &p4, Point_3 const &p5,
    	Point_3 const &p6,  Point_3 const &p7);

    __host__ __device__
    Point operator[](int const idx) const {
    	return _p[idx];
    }
    __host__ __device__
    Point &operator[](int const idx) {
    	return _p[idx];
    }

    // Draw prism.
    void draw(FILE *out_file, int const start_idx, int const resid);

    Point _p[8];
};
std::ostream& operator<<(std::ostream &stream, Prism const &t);

auto determinant(Vector const &v0, Vector const &v1, Vector const &v2) -> float;
auto determinant(Tetrahedron const &t) -> float;

}
#endif // _H
