#pragma once

#include <vector>

class SparsePointGrid
{
	struct PointData {
		PointData() : m_pos(), m_faceCount(-1) {}
		PointData(Vec3_arg pos) : m_pos(pos), m_faceCount(-1) {}
		Vec3 m_pos;
		int m_faceCount;	// -1: not used yet, 0: deleted, > 0: valid for searching
	};

	struct Cell {
		Cell() : m_numPoints(0), m_points(0) {}
		~Cell() { delete[] m_points; }
		int m_numPoints;
		int * m_points;

	private:
		Cell(const Cell&);
		Cell& operator=(const Cell&);
	};

	PointData * m_allPoints;
	int m_numPoints;
	AABB m_pointsAABB;
	Cell** m_cells; 

	float m_cellDim;
	int m_cellsPerDim;
	float m_gridDim;		// total length of grid side
	float m_minGridDim;		
	float m_maxGridDim;		

public:
	enum SplitDir
	{
		SPLITDIR_X = 0,
		SPLITDIR_Y,
		SPLITDIR_Z
	};

	enum ConstraintFlags
	{
		CONSTRAINT_POINT_ABOVE = (1 << 0),
		CONSTRAINT_DELAUNAY = (1 << 1)
	};

	SparsePointGrid(float gridDims, int cellsPerDim);
	~SparsePointGrid();

	void InsertPoints(const Vec3* points, int numPoints);
	int ClosestPointToSplit(SplitDir dir, const AABB& boundsToSplit);
	int NearestNeighborAcrossPlane(int from, const Plane& plane);
	int PointWithMinCircumcircle(int v0, int v1);
	int PointWithMinCircumsphere(const AABB& bounds, int v0, int v1, int v2, int flags = 0);
	Vec3 GetPos(int pointIdx) const;
	int GetNumPoints() const { return m_numPoints; }
	bool IsValidPoint(int idx) const { return m_allPoints[idx].m_faceCount != 0; }
	void FindPointsInRadius(Vec3_arg from, float radius, std::vector<int> &outPoints);
	bool HasPointsInRadiusSq(Vec3_arg from, float radiusSq, const int *except, int numExcept);

	void AddRef(int pointIdx);
	void SubRef(int pointIdx);

	const AABB& GetAllPointsAABB() const { return m_pointsAABB; }
private:
	SparsePointGrid(const SparsePointGrid&);
	SparsePointGrid& operator=(const SparsePointGrid&);
	void Clear();

	void ToGrid(Vec3_arg v, int &ix, int &iy, int &iz);
	void ToGrid(float value, int& index);
	void ToGrid(const AABB& bounds, Vec3_arg v, int &ix, int &iy, int &iz);
	void ToGridNoClamp(Vec3_arg v, int &ix, int &iy, int &iz);
	void ToGridNoClamp(float value, int& index);
	void FindClosestPointInCell(const Cell* cell, const Plane& plane, 
		int& closestPointIndex, float &closestPointDist);
	void FindClosestPointInCellAbovePlane(const Cell* cell, Vec3_arg fromPos, const Plane& plane,
		float fromPlaneDist, int& closestPointIndex, float &closestPointDistSq);

	
};

