#pragma once

//#include "assert.hh"
#include <cmath>

class Vec3;
typedef const Vec3& Vec3_arg;
class DVec3;
typedef const DVec3& DVec3_arg;

////////////////////////////////////////////////////////////////////////////////
// Float Vec3 
class Vec3
{
public:
	float x,y,z;

	inline Vec3() : x(), y(), z() {}
	inline Vec3(float x, float y, float z) : x(x),y(y),z(z) {}
	Vec3(DVec3_arg other);
	inline Vec3(Vec3_arg other) : x(other.x),y(other.y),z(other.z) {}
	inline Vec3(const float* v) : x(v[0]), y(v[1]), z(v[2]) {}
	inline Vec3& operator=(Vec3_arg  other)
		{
			x = other.x; y = other.y; z = other.z;
			return *this;
		}

	inline bool operator==(Vec3_arg  other) const 
		{
			return x == other.x && y == other.y && z == other.z;
		}
		
	inline bool operator!=(Vec3_arg  other) const
		{
			return !(*this == other);
		}

	inline Vec3 operator-() const
		{
			return Vec3(-x,-y,-z);
		}

	inline Vec3_arg operator-=(Vec3_arg other)
		{
			x-= other.x; y-= other.y; z-= other.z;
			return *this;
		}

	inline Vec3_arg  operator+=(Vec3_arg other)
		{
			x+= other.x; y+= other.y; z+= other.z;
			return *this;
		}

	inline Vec3_arg  operator*=(float s)
		{
			x*=s; y*=s; z*=s;
			return *this;
		}

	inline Vec3_arg  operator/=(float s)
		{
			x/=s;y/=s;z/=s;
			return *this;
		}

	inline float operator[](int idx) const {
		//ASSERT(idx >= 0 && idx < 3);
		return (&x)[idx];
	}
	
	inline float &operator[](int idx) {
		//ASSERT(idx >= 0 && idx < 3);
		return (&x)[idx];
	}

	inline void set(float x, float y, float z) 
		{ 
			this->x = x;
			this->y = y;
			this->z = z;
		}
};

////////////////////////////////////////////////////////////////////////////////
// Vec3 util
inline Vec3 operator+(Vec3_arg lhs, Vec3_arg rhs) 
{
	return Vec3(lhs.x + rhs.x,
				lhs.y + rhs.y,
				lhs.z + rhs.z);
}


inline Vec3 operator-(Vec3_arg lhs, Vec3_arg rhs) 
{
	return Vec3(lhs.x - rhs.x,
				lhs.y - rhs.y,
				lhs.z - rhs.z);
}


inline Vec3 operator*(Vec3_arg lhs, float s) 
{
	return Vec3(lhs.x*s,lhs.y*s,lhs.z*s);
}
		

inline Vec3 operator* (float scalar, Vec3_arg  other) 
{
	return other*scalar;
}

inline Vec3 operator/(Vec3_arg lhs, float s) 
{
	return Vec3(lhs.x/s, lhs.y/s, lhs.z/s);
}

inline Vec3 Lerp(float param,
                Vec3_arg  left,
				Vec3_arg  right)
{
	return (1.f - param) * left + param * right;
}

inline float dot(Vec3_arg lhs, Vec3_arg rhs)
{
	return 
		lhs.x * rhs.x + 
		lhs.y * rhs.y + 
		lhs.z * rhs.z;
}

inline float magnitude_squared(Vec3_arg vec)
{
	return dot(vec,vec);
}

inline float magnitude(Vec3_arg vec)
{
	return sqrt( magnitude_squared(vec) );
}

inline Vec3 normalize(Vec3_arg vec) 
{
	float mag = magnitude(vec);
	return vec/mag;
}

inline Vec3 normalize_safe(Vec3_arg vec, Vec3_arg defaultVec)
{
    float mag = magnitude(vec);
    if(mag > 1e-3f) {
        return vec/mag;
    } else {
        return defaultVec;
    }
}

inline Vec3 cross(Vec3_arg lhs, Vec3_arg rhs)
{
	return Vec3( lhs.y * rhs.z - lhs.z * rhs.y,
				 lhs.z * rhs.x - lhs.x * rhs.z,
				 lhs.x * rhs.y - lhs.y * rhs.x );
}

inline Vec3 vec_xz(Vec3_arg v)
{
    return Vec3(v.x, 0.0, v.z);
}

inline Vec3 line_normal_xz(Vec3_arg v)
{
    return Vec3(-v.z, 0, v.x);
}

inline Vec3 VecMin(Vec3_arg l, Vec3_arg r)
{
	return Vec3( 
		l.x < r.x ? l.x : r.x,
		l.y < r.y ? l.y : r.y,
		l.z < r.z ? l.z : r.z);
}

inline Vec3 VecMax(Vec3_arg l, Vec3_arg r)
{
	return Vec3( 
		l.x > r.x ? l.x : r.x,
		l.y > r.y ? l.y : r.y,
		l.z > r.z ? l.z : r.z);
}

////////////////////////////////////////////////////////////////////////////////
// Double vector3

class DVec3
{
public:
	double x,y,z;

	inline DVec3() : x(), y(), z() {}
	inline DVec3(double x, double y, double z) : x(x),y(y),z(z) {}
	inline DVec3(DVec3_arg other) : x(other.x),y(other.y),z(other.z) {}
	inline DVec3(Vec3_arg other) : x(other.x),y(other.y),z(other.z) {}
	inline DVec3(const double* v) : x(v[0]), y(v[1]), z(v[2]) {}
	inline DVec3& operator=(DVec3_arg  other)
		{
			x = other.x; y = other.y; z = other.z;
			return *this;
		}

	inline bool operator==(DVec3_arg  other) const 
		{
			return x == other.x && y == other.y && z == other.z;
		}
		
	inline bool operator!=(DVec3_arg  other) const
		{
			return !(*this == other);
		}

	inline DVec3 operator-() const
		{
			return DVec3(-x,-y,-z);
		}

	inline DVec3_arg operator-=(DVec3_arg other)
		{
			x-= other.x; y-= other.y; z-= other.z;
			return *this;
		}

	inline DVec3_arg  operator+=(DVec3_arg other)
		{
			x+= other.x; y+= other.y; z+= other.z;
			return *this;
		}

	inline DVec3_arg  operator*=(double s)
		{
			x*=s; y*=s; z*=s;
			return *this;
		}

	inline DVec3_arg  operator/=(double s)
		{
			x/=s;y/=s;z/=s;
			return *this;
		}

	inline double operator[](int idx) const {
		//ASSERT(idx >= 0 && idx < 3);
		return (&x)[idx];
	}
	
	inline double &operator[](int idx) {
		//ASSERT(idx >= 0 && idx < 3);
		return (&x)[idx];
	}

	inline void set(double x, double y, double z) 
		{ 
			this->x = x;
			this->y = y;
			this->z = z;
		}
};

////////////////////////////////////////////////////////////////////////////////
// DVec3 util

inline DVec3 operator+(DVec3_arg lhs, DVec3_arg rhs) 
{
	return DVec3(lhs.x + rhs.x,
				lhs.y + rhs.y,
				lhs.z + rhs.z);
}


inline DVec3 operator-(DVec3_arg lhs, DVec3_arg rhs) 
{
	return DVec3(lhs.x - rhs.x,
				lhs.y - rhs.y,
				lhs.z - rhs.z);
}


inline DVec3 operator*(DVec3_arg lhs, double s) 
{
	return DVec3(lhs.x*s,lhs.y*s,lhs.z*s);
}
		

inline DVec3 operator* (double scalar, DVec3_arg  other) 
{
	return other*scalar;
}

inline DVec3 operator/(DVec3_arg lhs, double s) 
{
	return DVec3(lhs.x/s, lhs.y/s, lhs.z/s);
}

inline DVec3 Lerp(double param,
                DVec3_arg  left,
				DVec3_arg  right)
{
	return (1. - param) * left + param * right;
}

inline double dot(DVec3_arg lhs, DVec3_arg rhs)
{
	return 
		lhs.x * rhs.x + 
		lhs.y * rhs.y + 
		lhs.z * rhs.z;
}

inline double magnitude_squared(DVec3_arg vec)
{
	return dot(vec,vec);
}

inline double magnitude(DVec3_arg vec)
{
	return sqrt( magnitude_squared(vec) );
}

inline DVec3 normalize(DVec3_arg vec) 
{
	double mag = magnitude(vec);
	return vec/mag;
}

inline DVec3 normalize_safe(DVec3_arg vec, DVec3_arg defaultVec)
{
    double mag = magnitude(vec);
    if(mag > 1e-3) {
        return vec/mag;
    } else {
        return defaultVec;
    }
}

inline DVec3 cross(DVec3_arg lhs, DVec3_arg rhs)
{
	return DVec3( lhs.y * rhs.z - lhs.z * rhs.y,
				 lhs.z * rhs.x - lhs.x * rhs.z,
				 lhs.x * rhs.y - lhs.y * rhs.x );
}

inline DVec3 vec_xz(DVec3_arg v)
{
    return DVec3(v.x, 0.0, v.z);
}

inline DVec3 line_normal_xz(DVec3_arg v)
{
    return DVec3(-v.z, 0, v.x);
}

inline DVec3 VecMin(DVec3_arg l, DVec3_arg r)
{
	return DVec3( 
		l.x < r.x ? l.x : r.x,
		l.y < r.y ? l.y : r.y,
		l.z < r.z ? l.z : r.z);
}

inline DVec3 VecMax(DVec3_arg l, DVec3_arg r)
{
	return DVec3( 
		l.x > r.x ? l.x : r.x,
		l.y > r.y ? l.y : r.y,
		l.z > r.z ? l.z : r.z);
}

////////////////////////////////////////////////////////////////////////////////
// Extra constructors
inline Vec3::Vec3(DVec3_arg other)
	: x(other.x), y(other.y), z(other.z)
{
}

