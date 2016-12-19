#pragma once

#include "vec3.hh"
#include "mat4.hh"
#include <cfloat>

////////////////////////////////////////////////////////////////////////////////
#define EPSILON 1e-3f
#define EPSILON_SQ (EPSILON * EPSILON)

////////////////////////////////////////////////////////////////////////////////
struct Plane
{
	Plane() : m_normal(), m_d() {}
	Plane(Vec3_arg normal, float d) : m_normal(normal), m_d(d) {}
	Vec3 m_normal;
	float m_d;
};

struct DPlane
{
	DPlane() : m_normal(), m_d() {}
	DPlane(DVec3_arg normal, double d) : m_normal(normal), m_d(d) {}
	DVec3 m_normal;
	double m_d;
};

struct OBB
{
	OBB() : m_center(), m_u(), m_v(), m_w() {}
	OBB(Vec3_arg center, Vec3_arg u, Vec3_arg v, Vec3_arg w)
		: m_center(center), m_u(u), m_v(v), m_w(w) {}
	Vec3 m_center;
	Vec3 m_u;
	Vec3 m_v;
	Vec3 m_w;
};

struct AABB
{
	AABB() : m_min(FLT_MAX, FLT_MAX, FLT_MAX), m_max(-FLT_MAX, -FLT_MAX, -FLT_MAX) {}
	AABB(Vec3_arg mn, Vec3_arg mx) : m_min(mn), m_max(mx) {}
	AABB(const AABB& o) : m_min(o.m_min), m_max(o.m_max) {}
	Vec3 m_min;
	Vec3 m_max;

	void Extend(Vec3_arg v) 
	{
		m_min = VecMin(v, m_min);
		m_max = VecMax(v, m_max);
	}

	void Extend(const OBB& obb)
	{	
		Extend(obb.m_center - obb.m_u - obb.m_v - obb.m_w);
		Extend(obb.m_center + obb.m_u - obb.m_v - obb.m_w);
		Extend(obb.m_center + obb.m_u + obb.m_v - obb.m_w);
		Extend(obb.m_center - obb.m_u + obb.m_v - obb.m_w);
		Extend(obb.m_center - obb.m_u - obb.m_v + obb.m_w);
		Extend(obb.m_center + obb.m_u - obb.m_v + obb.m_w);
		Extend(obb.m_center + obb.m_u + obb.m_v + obb.m_w);
		Extend(obb.m_center - obb.m_u + obb.m_v + obb.m_w);
	}
};

////////////////////////////////////////////////////////////////////////////////
enum TriPlaneIntersectType
{
	TRI_PLANE_INTERSECT_ALL_BELOW = 0x0,
	TRI_PLANE_INTERSECT_ALL_ABOVE = 0x07
};

////////////////////////////////////////////////////////////////////////////////
inline float DistToPlane(Vec3_arg p, const Plane& plane);
inline Vec3 ClosestPointOnAABBToPoint(const AABB& aabb, Vec3_arg pt);
inline Vec3 FurthestPointOnAABBToPoint(const AABB& aabb, Vec3_arg pt);
inline float DistSqAABBToPoint(const AABB& aabb, Vec3_arg pt);
inline bool AABBIntersectsSphere(const AABB& aabb, Vec3_arg sphereCenter, float sphereRadius);
inline bool AABBIntersectsShell(const AABB& aabb, Vec3_arg sphereCenter, float minRadius, float maxRadius);
inline bool AABBAbovePlane(const AABB& aabb, const Plane& plane);
inline bool AABBAbovePlane(const AABB& aabb, const DPlane& plane);
inline bool AABBContains(const AABB& aabb, Vec3_arg pt);
inline bool AABBContains(const AABB& aabb, DVec3_arg pt);
inline Vec3 MakeSplitNormal(int splitdir);
inline void MakeSplitPlane(Plane& plane, int dir, const AABB& bounds);
inline bool ComputeCircumcircle(Vec3_arg a, Vec3_arg b, Vec3_arg c, Vec3& outCenter, float &outRadiusSq);
inline bool ComputeCircumsphere(Vec3_arg a, Vec3_arg b, Vec3_arg c, Vec3_arg d, Vec3 &outCenter, float &outRadiusSq);
inline bool ComputeCircumsphere(DVec3_arg a, DVec3_arg b, DVec3_arg c, DVec3_arg d, DVec3 &outCenter, double &outRadiusSq);
inline void PlaneIntersectsTriangleList(const Plane& plane, int numTriangles, const Vec3* triangleData, int *results);
inline Vec3 ProjectN(Vec3_arg vec, Vec3_arg n);


float DistToPlane(Vec3_arg p, const Plane& plane)
{
    return fabs(dot(plane.m_normal, p) - plane.m_d);
}

Vec3 ClosestPointOnAABBToPoint(const AABB& aabb, Vec3_arg pt)
{
    Vec3 result = pt;
    result = VecMax(result, aabb.m_min);
    result = VecMin(result, aabb.m_max);
    return result;
}

Vec3 FurthestPointOnAABBToPoint(const AABB& aabb, Vec3_arg pt)
{
    Vec3 result = pt;
    result = VecMax(result, aabb.m_max);
    result = VecMin(result, aabb.m_min);
    return result;
}

float DistSqAABBToPoint(const AABB& aabb, Vec3_arg pt)
{
    float distSq = 0.f;
    for(int i = 0; i < 3; ++i)
    {
        if(pt[i] < aabb.m_min[i]) distSq += (aabb.m_min[i] - pt[i]) * (aabb.m_min[i] - pt[i]);
        if(pt[i] > aabb.m_max[i]) distSq += (pt[i] - aabb.m_max[i]) * (pt[i] - aabb.m_max[i]);
    }
    return distSq;
}

bool AABBIntersectsSphere(const AABB& aabb, Vec3_arg sphereCenter, float sphereRadius)
{
    float distSq = DistSqAABBToPoint(aabb, sphereCenter);
    return distSq < (sphereRadius * sphereRadius);
}

bool AABBIntersectsShell(const AABB& aabb, Vec3_arg sphereCenter, float minRadius, float maxRadius)
{
    float distSq = DistSqAABBToPoint(aabb, sphereCenter);
    if(distSq < (maxRadius * maxRadius))
    {
        Vec3 furthestPt = FurthestPointOnAABBToPoint(aabb, sphereCenter);
        distSq = magnitude_squared(furthestPt - sphereCenter);
        return distSq > (minRadius * minRadius);
    }
    return false;
}

bool AABBAbovePlane(const AABB& aabb, const Plane& plane)
{
    Vec3 diff = 0.5f * aabb.m_max - 0.5f * aabb.m_min;
    float radius = fabs(diff[0]) + fabs(diff[1]) + fabs(diff[2]);
    Vec3 center = aabb.m_min + diff;
    float planeDist = dot(plane.m_normal, center) - plane.m_d;
    return planeDist >= -radius;
}

bool AABBAbovePlane(const AABB& aabb, const DPlane& plane)
{
    DVec3 diff = 0.5 * DVec3(aabb.m_max) - 0.5 * DVec3(aabb.m_min);
    double radius = fabs(diff[0]) + fabs(diff[1]) + fabs(diff[2]);
    DVec3 center = DVec3(aabb.m_min) + diff;
    double planeDist = dot(plane.m_normal, center) - plane.m_d;
    return planeDist >= -radius;
}


Vec3 MakeSplitNormal(int splitdir)
{
    return Vec3(
                float(splitdir == 0),
                float(splitdir == 1),
                float(splitdir == 2));
}

void MakeSplitPlane(Plane& plane, int dir, const AABB& bounds)
{
    plane.m_normal = MakeSplitNormal(dir);
    plane.m_d = dot(0.5f*(bounds.m_max + bounds.m_min), plane.m_normal);
}

bool ComputeCircumcircle(Vec3_arg a, Vec3_arg b, Vec3_arg c, Vec3& outCenter, float &outRadiusSq)
{
    // Solving with Cramer's rule - just solving one t-parameter of two lines intersecting.
    Vec3 ab = b - a;
    Vec3 bc = c - b;
    Vec3 abHalf = 0.5f * a + 0.5f * b;
    Vec3 bcHalf = 0.5f * b + 0.5f * c;
    Vec3 midpointDiff = bcHalf - abHalf;
    
    Vec3 normal = cross(ab, bc);
    Vec3 abNormal = cross(normal, ab);
    Vec3 bcNormal = cross(normal, bc);
    
    float dot0 = dot(abNormal, abNormal);
    float dot1 = dot(abNormal, bcNormal);
    float dot3 = -dot(bcNormal, bcNormal);
    
    float denom = dot0 * dot3 + dot1 * dot1; // dot1 * dot1 = - (-dot1) * (dot1)
    if( fabs(denom) < EPSILON )
        return false;
    
    float dot4 = dot(abNormal, midpointDiff);
    float dot5 = dot(bcNormal, midpointDiff);
    
    float t = ( dot4 * dot3 + dot1 * dot5 ) / denom; // + dot 1 = - (-dot1)
    outCenter = abHalf + t * abNormal ;
    outRadiusSq = magnitude_squared(a - outCenter);
    
    return true;
    
}

bool ComputeCircumsphere(Vec3_arg a, Vec3_arg b, Vec3_arg c, Vec3_arg d, Vec3 &outCenter, float &outRadiusSq)
{
    Mat4 aMat(a.x, a.y, a.z, 1.f,
              b.x, b.y, b.z, 1.f,
              c.x, c.y, c.z, 1.f,
              d.x, d.y, d.z, 1.f);
    double aCoef = ddet(aMat);
    
    if(fabs(aCoef) < EPSILON)
        return false;
    
    Mat4 bMat = aMat;
    bMat.m[0] = a.x * a.x + a.y * a.y + a.z * a.z;
    bMat.m[4] = b.x * b.x + b.y * b.y + b.z * b.z;
    bMat.m[8] = c.x * c.x + c.y * c.y + c.z * c.z;
    bMat.m[12] = d.x * d.x + d.y * d.y + d.z * d.z;
    double bxCoef = ddet(bMat);
    
    bMat.m[1] = a.x;
    bMat.m[5] = b.x;
    bMat.m[9] = c.x;
    bMat.m[13] = d.x;
    double byCoef = -ddet(bMat);
    
    bMat.m[2] = a.y;
    bMat.m[6] = b.y;
    bMat.m[10] = c.y;
    bMat.m[14] = d.y;
    double bzCoef = ddet(bMat);
    
    bMat.m[3] = a.z;
    bMat.m[7] = b.z;
    bMat.m[11] = c.z;
    bMat.m[15] = d.z;
    double cCoef = ddet(bMat);
    
    double inv_double_a = 1.0 / (2.0 * aCoef);
    outCenter.x = bxCoef * inv_double_a;
    outCenter.y = byCoef * inv_double_a;
    outCenter.z = bzCoef * inv_double_a;
    
    outRadiusSq = float((bxCoef * bxCoef + byCoef * byCoef + bzCoef * bzCoef - 4.f * aCoef * cCoef) / (4.f * aCoef * aCoef));
    
    return true;
}

bool ComputeCircumsphere(DVec3_arg a, DVec3_arg b, DVec3_arg c, DVec3_arg d, DVec3 &outCenter, double &outRadiusSq)
{
    DMat4 aMat(a.x, a.y, a.z, 1.f,
               b.x, b.y, b.z, 1.f,
               c.x, c.y, c.z, 1.f,
               d.x, d.y, d.z, 1.f);
    double aCoef = det(aMat);
    
    if(std::abs(aCoef) < EPSILON)
        return false;
    
    DMat4 bMat = aMat;
    bMat.m[0] = a.x * a.x + a.y * a.y + a.z * a.z;
    bMat.m[4] = b.x * b.x + b.y * b.y + b.z * b.z;
    bMat.m[8] = c.x * c.x + c.y * c.y + c.z * c.z;
    bMat.m[12] = d.x * d.x + d.y * d.y + d.z * d.z;
    double bxCoef = det(bMat);
    
    bMat.m[1] = a.x;
    bMat.m[5] = b.x;
    bMat.m[9] = c.x;
    bMat.m[13] = d.x;
    double byCoef = -det(bMat);
    
    bMat.m[2] = a.y;
    bMat.m[6] = b.y;
    bMat.m[10] = c.y;
    bMat.m[14] = d.y;
    double bzCoef = det(bMat);
    
    bMat.m[3] = a.z;
    bMat.m[7] = b.z;
    bMat.m[11] = c.z;
    bMat.m[15] = d.z;
    double cCoef = det(bMat);
    
    double inv_double_a = 1.0 / (2.0 * aCoef);
    outCenter.x = bxCoef * inv_double_a;
    outCenter.y = byCoef * inv_double_a;
    outCenter.z = bzCoef * inv_double_a;
    
    outRadiusSq = (bxCoef * bxCoef + byCoef * byCoef + bzCoef * bzCoef - 4.f * aCoef * cCoef) / (4.f * aCoef * aCoef);
    
    return true;
}

void PlaneIntersectsTriangleList(const Plane& plane, int numTriangles, const Vec3* triangleData, int *results)
{
    int triangleOff = 0;
    for(int i = 0; i < numTriangles; ++i)
    {
        int v0r = (int(dot(plane.m_normal, triangleData[triangleOff]) - plane.m_d >= 0.f)) & 1;
        int v1r = (int(dot(plane.m_normal, triangleData[triangleOff+1]) - plane.m_d >= 0.f) & 1) << 1 ;
        int v2r = (int(dot(plane.m_normal, triangleData[triangleOff+2]) - plane.m_d >= 0.f) & 1) << 2 ;
        
        results[i] = v0r | v1r | v2r;
        triangleOff += 3;
    }
}

bool AABBContains(const AABB& aabb, Vec3_arg pt)
{
    return (pt.x >= aabb.m_min.x && pt.x <= aabb.m_max.x) &&
    (pt.y >= aabb.m_min.y && pt.y <= aabb.m_max.y) &&
    (pt.z >= aabb.m_min.z && pt.z <= aabb.m_max.z);
}

bool AABBContains(const AABB& aabb, DVec3_arg pt)
{
    return (pt.x >= double(aabb.m_min.x) && pt.x <= double(aabb.m_max.x)) &&
    (pt.y >= double(aabb.m_min.y) && pt.y <= double(aabb.m_max.y)) &&
    (pt.z >= double(aabb.m_min.z) && pt.z <= double(aabb.m_max.z));
}

Vec3 ProjectN(Vec3_arg vec, Vec3_arg n)
{
    return vec - dot(vec, n) * n;
}


////////////////////////////////////////////////////////////////////////////////
template<class T>
inline T Clamp(T val, T min, T max)
{
	if(val < min)
		return min;
	else if(val > max)
		return max;
	else
		return val;
}

template<class T>
inline T Max(T val0, T val1)
{
	return val0 > val1 ? val0 : val1;
}

template<class T>
inline T Min(T val0, T val1)
{
	return val0 < val1 ? val0 : val1;
}

template<class T>
inline void Swap(T& val0, T& val1)
{
	T temp = val0;
	val0 = val1;
	val1 = temp;
}

inline float AngleWrap(float angle)
{
	while(angle < 0.f)
		angle += 2.f * M_PI;
	while(angle > 2.f * M_PI)
		angle -= 2.f * M_PI;
	return angle;
}


