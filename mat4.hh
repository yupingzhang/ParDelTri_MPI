#pragma once

//#include "assert.hh"
#include <cstring>

////////////////////////////////////////////////////////////////////////////////
/* Matrix math conventions

    Multiplication with points and vectors is post-multiply. ie:

    p' = M p
    v' = M v
*/

class Mat4;
typedef const Mat4& Mat4_arg;
class DMat4;
typedef const DMat4& DMat4_arg;

inline Mat4 operator*(Mat4_arg lhs, Mat4_arg rhs) ;
inline DMat4 operator*(DMat4_arg lhs, DMat4_arg rhs) ;
////////////////////////////////////////////////////////////////////////////////
// Fat 4x4 matrix.
class Mat4
{
public:
	float m[16];

	inline Mat4( ) {}		
	Mat4( float v ) {
		for(int i = 0; i < 16; ++i) {
			m[i] = v;
		}
	}

	struct ident_t { };
	Mat4( ident_t ) {	// use like: Mat4( (Mat4::ident_t()) )
		for(int i = 0; i < 16; ++i) {
			m[i] = 0.0;
		}
		m[0] = 1.f; m[5] = 1.f; m[10] = 1.f; m[15] = 1.f;
	}	

	Mat4( float a, float b, float c, float d,
			float e, float f, float g, float h,
			float i, float j, float k, float l,
			float mm, float n, float o, float p)
	{
		m[0] = a; m[1] = b; m[2] = c; m[3] = d;
		m[4] = e; m[5] = f; m[6] = g; m[7] = h;
		m[8] = i; m[9] = j; m[10] = k; m[11] = l;
		m[12] = mm; m[13] = n; m[14] = o; m[15] = p;
	}

	Mat4( float v[16] ) 
	{
		memcpy(m,v,sizeof(float)*16);
	}

	Mat4(Mat4_arg  other) {
		memcpy(m,other.m,sizeof(m));
	}

	Mat4& operator=(Mat4_arg  other) {
		if(this != &other) {
			memcpy(m,other.m,sizeof(m));
		} 
		return *this;
	}

	Mat4& operator*=(Mat4_arg  other) 
	{
		*this = (*this * other);
		return *this;
	}

	Mat4& operator+=(Mat4_arg  other) 
	{
		for(int i = 0; i < 16; ++i) {
			m[i] += other.m[i];
		}
		return *this;
	}

	inline Vec3 Row3(int rownum) const
	{
		//ASSERT(rownum >= 0 && rownum < 4);
		return Vec3(&m[rownum*4]);
	}

	inline Vec3 Col3(int colnum) const
	{
		//ASSERT(colnum >= 0 && colnum < 4);
		return Vec3(m[colnum],m[colnum+4],m[colnum+8]);
	}
};

inline Mat4 operator+(Mat4_arg lhs, Mat4_arg  rhs) 
{
	Mat4 result;
	for(int i = 0; i < 16; ++i) {
		for(int j = 0; j < 16; ++j) {
			result.m[i] = lhs.m[i] + rhs.m[i];
		}
	}
	return result;
}

inline Mat4 operator*(Mat4_arg lhs, Mat4_arg rhs) 
{
	Mat4 result;
				
	int resultIdx = 0;
	for(int rowIdx = 0; rowIdx < 16; rowIdx += 4)
	{
		for(int col = 0; col < 4; ++col)
		{
			result.m[resultIdx] = lhs.m[rowIdx] * rhs.m[col]
				+ lhs.m[rowIdx + 1] * rhs.m[col + 4]
				+ lhs.m[rowIdx + 2] * rhs.m[col + 8]
				+ lhs.m[rowIdx + 3] * rhs.m[col + 12];
			++resultIdx;
		}
	}
	return result;
}

inline float det3x3(float a11, float a12, float a13,
					float a21, float a22, float a23,
					float a31, float a32, float a33) 
{
	float d1 = a22*a33 - a23*a32;
	float d2 = a23*a31 - a21*a33; // negated so orders is reversed
	float d3 = a21*a32 - a22*a31;
	
	float r = a11 * d1 + a12 * d2 + a13 * d3;
	return r;
}

inline float det(Mat4_arg m)
{	
	float d1 = det3x3(m.m[1*4 + 1], m.m[1*4 + 2], m.m[1*4 + 3],
					  m.m[2*4 + 1], m.m[2*4 + 2], m.m[2*4 + 3],
					  m.m[3*4 + 1], m.m[3*4 + 2], m.m[3*4 + 3]);
	float d2 = det3x3(m.m[1*4 + 0], m.m[1*4 + 2], m.m[1*4 + 3],
					  m.m[2*4 + 0], m.m[2*4 + 2], m.m[2*4 + 3],
					  m.m[3*4 + 0], m.m[3*4 + 2], m.m[3*4 + 3]);
	float d3 = det3x3(m.m[1*4 + 0], m.m[1*4 + 1], m.m[1*4 + 3],
					  m.m[2*4 + 0], m.m[2*4 + 1], m.m[2*4 + 3],
					  m.m[3*4 + 0], m.m[3*4 + 1], m.m[3*4 + 3]);
	float d4 = det3x3(m.m[1*4 + 0], m.m[1*4 + 1], m.m[1*4 + 2],
					  m.m[2*4 + 0], m.m[2*4 + 1], m.m[2*4 + 2],
					  m.m[3*4 + 0], m.m[3*4 + 1], m.m[3*4 + 2]);
				
	float r = m.m[0*4 + 0] * d1 - m.m[0*4 + 1] * d2 + m.m[0*4 + 2] * d3 - m.m[0*4 + 3] * d4;
	return r;	
}

inline double ddet3x3(double a11, double a12, double a13,
					double a21, double a22, double a23,
					double a31, double a32, double a33) 
{
	double d1 = a22*a33 - a23*a32;
	double d2 = a23*a31 - a21*a33; // negated so order is reversed
	double d3 = a21*a32 - a22*a31;
	
	double r = a11 * d1 + a12 * d2 + a13 * d3;
	return r;
}

inline double ddet(Mat4_arg m)
{	
	double d1 = ddet3x3(m.m[1*4 + 1], m.m[1*4 + 2], m.m[1*4 + 3],
					  m.m[2*4 + 1], m.m[2*4 + 2], m.m[2*4 + 3],
					  m.m[3*4 + 1], m.m[3*4 + 2], m.m[3*4 + 3]);
	double d2 = ddet3x3(m.m[1*4 + 0], m.m[1*4 + 2], m.m[1*4 + 3],
					  m.m[2*4 + 0], m.m[2*4 + 2], m.m[2*4 + 3],
					  m.m[3*4 + 0], m.m[3*4 + 2], m.m[3*4 + 3]);
	double d3 = ddet3x3(m.m[1*4 + 0], m.m[1*4 + 1], m.m[1*4 + 3],
					  m.m[2*4 + 0], m.m[2*4 + 1], m.m[2*4 + 3],
					  m.m[3*4 + 0], m.m[3*4 + 1], m.m[3*4 + 3]);
	double d4 = ddet3x3(m.m[1*4 + 0], m.m[1*4 + 1], m.m[1*4 + 2],
					  m.m[2*4 + 0], m.m[2*4 + 1], m.m[2*4 + 2],
					  m.m[3*4 + 0], m.m[3*4 + 1], m.m[3*4 + 2]);
				
	double r = m.m[0*4 + 0] * d1 - m.m[0*4 + 1] * d2 + m.m[0*4 + 2] * d3 - m.m[0*4 + 3] * d4;
	return r;	
}

inline Mat4 transpose(Mat4_arg m)
{
	Mat4 result = m;
	for(int i = 0; i < 4; ++i) {
		for(int j = 0; j < i; ++j) {
			int firstIdx = j + 4*i;
			int secondIdx = i + 4*j;
			float t = result.m[firstIdx];
			result.m[firstIdx] = result.m[secondIdx];
			result.m[secondIdx] = t;
		}
	}
	return result;
}

inline Vec3 transform_point(Mat4_arg lhs, Vec3_arg rhs)
{
	float x = lhs.m[0] * rhs.x + lhs.m[1] * rhs.y + lhs.m[2] * rhs.z + lhs.m[3];
	float y = lhs.m[4] * rhs.x + lhs.m[5] * rhs.y + lhs.m[6] * rhs.z + lhs.m[7];
	float z = lhs.m[8] * rhs.x + lhs.m[9] * rhs.y + lhs.m[10] * rhs.z + lhs.m[11];
	return Vec3(x,y,z);
}

inline Vec3 transform_vector(Mat4_arg lhs, Vec3_arg rhs)
{
	float x = lhs.m[0] * rhs.x + lhs.m[1] * rhs.y + lhs.m[2] * rhs.z;
	float y = lhs.m[4] * rhs.x + lhs.m[5] * rhs.y + lhs.m[6] * rhs.z;
	float z = lhs.m[8] * rhs.x + lhs.m[9] * rhs.y + lhs.m[10] * rhs.z;
	return Vec3(x,y,z);
}

inline Mat4 rotation_x(float rad)
{
	float a = cos(rad);
	float b = sin(rad);
	return Mat4(1.0f, 0.0f, 0.0f, 0.0f,
				0.0f, a,    -b,    0.0f,
				0.0f, b,   a,    0.0f,
				0.0f, 0.0f, 0.0f, 1.0f);
}

inline Mat4 rotation_y(float rad)
{
	float a = cos(rad);
	float b = sin(rad);
	return Mat4(a,    0.0f, b, 0.0f,
				0.0f, 1.0f, 0.0f, 0.0f,
				-b,   0.0f, a,    0.0f,
				0.0f, 0.0f, 0.0f, 1.0f);
}

inline Mat4 rotation_z(float rad)
{
	float a = cos(rad);
	float b = sin(rad);
	return Mat4 (a,    -b,    0.0f, 0.0f,
				 b,   a,    0.0f, 0.0f,
				 0.0f, 0.0f, 1.0f, 0.0f,
				 0.0f, 0.0f, 0.0f, 1.0f);
	
}

inline Mat4 translation( Vec3_arg t)
{
	return Mat4 (1.f, 0.f, 0.f, t.x,
				 0.f, 1.f, 0.f, t.y,
				 0.f, 0.f, 1.f, t.z,
				 0.f, 0.f, 0.f, 1.f );	
}
			

inline bool is_orthonormal(Mat4_arg  m)
{
	Vec3 row0( &m.m[0]);
	Vec3 row1( &m.m[4]);
	Vec3 row2( &m.m[8]);
		
	float len ;
	len = magnitude(row0);
	if( fabs(1.0f - len) > 1e-3f) return false;
	len = magnitude(row1);
	if( fabs(1.0f - len) > 1e-3f) return false;
	len = magnitude(row2);
	if( fabs(1.0f - len) > 1e-3f) return false;
		
	if(dot(row0,row1) > 1e-3f) return false;
	if(dot(row1,row2) > 1e-3f) return false;
	if(dot(row2,row0) > 1e-3f) return false;
	return true;

}

inline void matrix_multiply(float out[4], Mat4_arg  m, float const in[4])
{
	for(int row = 0, rowIdx = 0; rowIdx < 16; rowIdx+=4, ++row)
	{
		out[row] = m.m[rowIdx + 0] * in[0] + 
			m.m[rowIdx + 1] * in[1] +
			m.m[rowIdx + 2] * in[2] +
			m.m[rowIdx + 3] * in[3];
	}
}

inline Mat4 get_rotation_part4x4(Mat4_arg m)
{
	Mat4 result = m;
	result.m[3] = 0.f;
	result.m[7] = 0.f;
	result.m[11] = 0.f;
	result.m[12] = result.m[13] = result.m[14] = 0.f; result.m[15] = 1.f;
	return result;
}

// Invert a matrix that is a rotation and translation
inline Mat4 inverse_rot_trans(Mat4_arg m)
{
	Mat4 result = transpose(get_rotation_part4x4(m));
	Vec3 inv_trans = -m.Col3(3);

	inv_trans = transform_vector(result, inv_trans);
	result.m[3] = inv_trans.x;
	result.m[7] = inv_trans.y;
	result.m[11] = inv_trans.z;
	return result;
}


////////////////////////////////////////////////////////////////////////////////
// Double 4x4 Mat

class DMat4
{
public:
	double m[16];

	inline DMat4( ) {}		
	DMat4( double v ) {
		for(int i = 0; i < 16; ++i) {
			m[i] = v;
		}
	}

	struct ident_t { };
	DMat4( ident_t ) {	// use like: DMat4( (DMat4::ident_t()) )
		for(int i = 0; i < 16; ++i) {
			m[i] = 0.0;
		}
		m[0] = 1.; m[5] = 1.; m[10] = 1.; m[15] = 1.;
	}	

	DMat4( double a, double b, double c, double d,
			double e, double f, double g, double h,
			double i, double j, double k, double l,
			double mm, double n, double o, double p)
	{
		m[0] = a; m[1] = b; m[2] = c; m[3] = d;
		m[4] = e; m[5] = f; m[6] = g; m[7] = h;
		m[8] = i; m[9] = j; m[10] = k; m[11] = l;
		m[12] = mm; m[13] = n; m[14] = o; m[15] = p;
	}

	DMat4( double v[16] ) 
	{
		memcpy(m,v,sizeof(double)*16);
	}

	DMat4(DMat4_arg  other) {
		memcpy(m,other.m,sizeof(m));
	}

	DMat4& operator=(DMat4_arg  other) {
		if(this != &other) {
			memcpy(m,other.m,sizeof(m));
		} 
		return *this;
	}

	DMat4& operator*=(DMat4_arg  other) 
	{
		*this = (*this * other);
		return *this;
	}

	DMat4& operator+=(DMat4_arg  other) 
	{
		for(int i = 0; i < 16; ++i) {
			m[i] += other.m[i];
		}
		return *this;
	}

	inline DVec3 Row3(int rownum) const
	{
		//ASSERT(rownum >= 0 && rownum < 4);
		return DVec3(&m[rownum*4]);
	}

	inline DVec3 Col3(int colnum) const
	{
		//ASSERT(colnum >= 0 && colnum < 4);
		return DVec3(m[colnum],m[colnum+4],m[colnum+8]);
	}
};

inline DMat4 operator+(DMat4_arg lhs, DMat4_arg  rhs) 
{
	DMat4 result;
	for(int i = 0; i < 16; ++i) {
		for(int j = 0; j < 16; ++j) {
			result.m[i] = lhs.m[i] + rhs.m[i];
		}
	}
	return result;
}

inline DMat4 operator*(DMat4_arg lhs, DMat4_arg rhs) 
{
	DMat4 result;
				
	int resultIdx = 0;
	for(int rowIdx = 0; rowIdx < 16; rowIdx += 4)
	{
		for(int col = 0; col < 4; ++col)
		{
			result.m[resultIdx] = lhs.m[rowIdx] * rhs.m[col]
				+ lhs.m[rowIdx + 1] * rhs.m[col + 4]
				+ lhs.m[rowIdx + 2] * rhs.m[col + 8]
				+ lhs.m[rowIdx + 3] * rhs.m[col + 12];
			++resultIdx;
		}
	}
	return result;
}

inline double det(DMat4_arg m)
{	
	double d1 = ddet3x3(m.m[1*4 + 1], m.m[1*4 + 2], m.m[1*4 + 3],
					  m.m[2*4 + 1], m.m[2*4 + 2], m.m[2*4 + 3],
					  m.m[3*4 + 1], m.m[3*4 + 2], m.m[3*4 + 3]);
	double d2 = ddet3x3(m.m[1*4 + 0], m.m[1*4 + 2], m.m[1*4 + 3],
					  m.m[2*4 + 0], m.m[2*4 + 2], m.m[2*4 + 3],
					  m.m[3*4 + 0], m.m[3*4 + 2], m.m[3*4 + 3]);
	double d3 = ddet3x3(m.m[1*4 + 0], m.m[1*4 + 1], m.m[1*4 + 3],
					  m.m[2*4 + 0], m.m[2*4 + 1], m.m[2*4 + 3],
					  m.m[3*4 + 0], m.m[3*4 + 1], m.m[3*4 + 3]);
	double d4 = ddet3x3(m.m[1*4 + 0], m.m[1*4 + 1], m.m[1*4 + 2],
					  m.m[2*4 + 0], m.m[2*4 + 1], m.m[2*4 + 2],
					  m.m[3*4 + 0], m.m[3*4 + 1], m.m[3*4 + 2]);
				
	double r = m.m[0*4 + 0] * d1 - m.m[0*4 + 1] * d2 + m.m[0*4 + 2] * d3 - m.m[0*4 + 3] * d4;
	return r;	
}

inline DMat4 transpose(DMat4_arg m)
{
	DMat4 result = m;
	for(int i = 0; i < 4; ++i) {
		for(int j = 0; j < i; ++j) {
			int firstIdx = j + 4*i;
			int secondIdx = i + 4*j;
			double t = result.m[firstIdx];
			result.m[firstIdx] = result.m[secondIdx];
			result.m[secondIdx] = t;
		}
	}
	return result;
}

inline DVec3 transform_point(DMat4_arg lhs, DVec3_arg rhs)
{
	double x = lhs.m[0] * rhs.x + lhs.m[1] * rhs.y + lhs.m[2] * rhs.z + lhs.m[3];
	double y = lhs.m[4] * rhs.x + lhs.m[5] * rhs.y + lhs.m[6] * rhs.z + lhs.m[7];
	double z = lhs.m[8] * rhs.x + lhs.m[9] * rhs.y + lhs.m[10] * rhs.z + lhs.m[11];
	return DVec3(x,y,z);
}

inline DVec3 transform_vector(DMat4_arg lhs, DVec3_arg rhs)
{
	double x = lhs.m[0] * rhs.x + lhs.m[1] * rhs.y + lhs.m[2] * rhs.z;
	double y = lhs.m[4] * rhs.x + lhs.m[5] * rhs.y + lhs.m[6] * rhs.z;
	double z = lhs.m[8] * rhs.x + lhs.m[9] * rhs.y + lhs.m[10] * rhs.z;
	return DVec3(x,y,z);
}

inline DMat4 rotation_x(double rad)
{
	double a = cos(rad);
	double b = sin(rad);
	return DMat4(1.0, 0.0, 0.0, 0.0,
				0.0, a,    -b,    0.0,
				0.0, b,   a,    0.0,
				0.0, 0.0, 0.0, 1.0);
}

inline DMat4 rotation_y(double rad)
{
	double a = cos(rad);
	double b = sin(rad);
	return DMat4(a,    0.0, b, 0.0,
				0.0, 1.0, 0.0, 0.0,
				-b,   0.0, a,    0.0,
				0.0, 0.0, 0.0, 1.0);
}

inline DMat4 rotation_z(double rad)
{
	double a = cos(rad);
	double b = sin(rad);
	return DMat4 (a,    -b,    0.0, 0.0,
				 b,   a,    0.0, 0.0,
				 0.0, 0.0, 1.0, 0.0,
				 0.0, 0.0, 0.0, 1.0);
	
}

inline DMat4 translation( DVec3_arg t)
{
	return DMat4 (1., 0., 0., t.x,
				 0., 1., 0., t.y,
				 0., 0., 1., t.z,
				 0., 0., 0., 1. );	
}
			

inline bool is_orthonormal(DMat4_arg  m)
{
	DVec3 row0( &m.m[0]);
	DVec3 row1( &m.m[4]);
	DVec3 row2( &m.m[8]);
		
	double len ;
	len = magnitude(row0);
	if( fabs(1.0 - len) > 1e-3) return false;
	len = magnitude(row1);
	if( fabs(1.0 - len) > 1e-3) return false;
	len = magnitude(row2);
	if( fabs(1.0 - len) > 1e-3) return false;
		
	if(dot(row0,row1) > 1e-3) return false;
	if(dot(row1,row2) > 1e-3) return false;
	if(dot(row2,row0) > 1e-3) return false;
	return true;

}

inline void matrix_multiply(double out[4], DMat4_arg  m, double const in[4])
{
	for(int row = 0, rowIdx = 0; rowIdx < 16; rowIdx+=4, ++row)
	{
		out[row] = m.m[rowIdx + 0] * in[0] + 
			m.m[rowIdx + 1] * in[1] +
			m.m[rowIdx + 2] * in[2] +
			m.m[rowIdx + 3] * in[3];
	}
}

inline DMat4 get_rotation_part4x4(DMat4_arg m)
{
	DMat4 result = m;
	result.m[3] = 0.;
	result.m[7] = 0.;
	result.m[11] = 0.;
	result.m[12] = result.m[13] = result.m[14] = 0.; result.m[15] = 1.;
	return result;
}

// Invert a matrix that is a rotation and translation
inline DMat4 inverse_rot_trans(DMat4_arg m)
{
	DMat4 result = transpose(get_rotation_part4x4(m));
	DVec3 inv_trans = -m.Col3(3);

	inv_trans = transform_vector(result, inv_trans);
	result.m[3] = inv_trans.x;
	result.m[7] = inv_trans.y;
	result.m[11] = inv_trans.z;
	return result;
}

