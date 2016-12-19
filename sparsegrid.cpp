#include <cstdio>
#include <climits>
//#include "assert.hh"
#include "math.hh"
#include "sparsegrid.hh"
//#include "debugdraw.hh"

SparsePointGrid::SparsePointGrid(float gridDims, int cellsPerDim)
	: m_allPoints(0)
	, m_numPoints(0)
	, m_pointsAABB()
	, m_cells(0)
	, m_cellDim(gridDims / cellsPerDim)
	, m_cellsPerDim(cellsPerDim)
	, m_gridDim(gridDims)
	, m_minGridDim(-0.5f * gridDims)
	, m_maxGridDim(0.5f * gridDims)
{
	m_cells = new Cell*[ cellsPerDim * cellsPerDim * cellsPerDim ];
	memset(m_cells, 0, sizeof(Cell*) * cellsPerDim * cellsPerDim * cellsPerDim);
}

SparsePointGrid::~SparsePointGrid()
{
	Clear();
}

void SparsePointGrid::Clear()
{
	int zOff = 0;
	for(int z = 0; z < m_cellsPerDim; ++z)
	{	
		int yOff = 0;
		for(int y = 0; y < m_cellsPerDim; ++y)
		{
			for(int x = 0; x < m_cellsPerDim; ++x)
			{
				Cell* cell = m_cells[x + yOff + zOff];
				delete cell;
			}
			yOff += m_cellsPerDim;
		}
		zOff += m_cellsPerDim * m_cellsPerDim;
	}

	delete[] m_cells;
	delete[] m_allPoints;

	m_cells = 0;
	m_allPoints = 0;
}
	
void SparsePointGrid::ToGrid(float value, int& index)
{
	index = Clamp(int(floorf((value - m_minGridDim) / m_cellDim)), 0, m_cellsPerDim - 1);
}
	
void SparsePointGrid::ToGrid(Vec3_arg v, int &ix, int &iy, int &iz)
{
	ToGrid(v.x, ix);
	ToGrid(v.y, iy);
	ToGrid(v.z, iz);
}
	
void SparsePointGrid::ToGrid(const AABB& bounds, Vec3_arg v, int &ix, int &iy, int &iz)
{
	ToGridNoClamp(v.x, ix);
	ToGridNoClamp(v.y, iy);
	ToGridNoClamp(v.z, iz);

	int iMin[3], iMax[3];
	ToGrid(bounds.m_min, iMin[0], iMin[1], iMin[2]);
	ToGrid(bounds.m_max, iMax[0], iMax[1], iMax[2]);

	ix = Clamp(ix, iMin[0], iMax[0]);
	iy = Clamp(iy, iMin[1], iMax[1]);
	iz = Clamp(iz, iMin[2], iMax[2]);
}

void SparsePointGrid::ToGridNoClamp(float value, int& index)
{
	index = int(floorf((value - m_minGridDim) / m_cellDim));
}
	
void SparsePointGrid::ToGridNoClamp(Vec3_arg v, int &ix, int &iy, int &iz)
{
	ToGridNoClamp(v.x, ix);
	ToGridNoClamp(v.y, iy);
	ToGridNoClamp(v.z, iz);
}


void SparsePointGrid::InsertPoints(const Vec3* points, int numPoints)
{
	// helper values for insertion
	int zDimMult = m_cellsPerDim * m_cellsPerDim;
	int yDimMult = m_cellsPerDim;

	int maxIndex = m_cellsPerDim;

	int numPointsInserted = 0;
	// accumulate cells.
	int ix, iy, iz;
	for(int i = 0; i < numPoints; ++i)
	{
		ToGridNoClamp(points[i], ix, iy, iz);
		
		if(ix >= 0 && ix < maxIndex &&
			iy >= 0 && iy < maxIndex &&
			iz >= 0 && iz < maxIndex)
		{
			int cellIdx = ix + yDimMult * iy + zDimMult * iz;
			//ASSERT(cellIdx >= 0 && cellIdx < (m_cellsPerDim * m_cellsPerDim * m_cellsPerDim));
			Cell* cell = m_cells[cellIdx];
			if(!cell)
			{
				m_cells[cellIdx] = cell = new Cell();
				AABB cellBounds;
				cellBounds.m_min = Vec3(ix * m_cellDim + m_minGridDim, 
						iy * m_cellDim + m_minGridDim, 
						iz * m_cellDim + m_minGridDim);
				cellBounds.m_max = cellBounds.m_min + Vec3(m_cellDim, m_cellDim, m_cellDim);
				//DebugDrawAABB(cellBounds);
			}
				
			++cell->m_numPoints;
			++numPointsInserted;
		}
	}

	m_allPoints = new PointData[numPointsInserted];
	m_numPoints = numPointsInserted;
	m_pointsAABB = AABB();

	// insert points and allocate cells
	int *curCounts = new int[m_cellsPerDim * m_cellsPerDim * m_cellsPerDim];
	memset(curCounts, 0, sizeof(int) * m_cellsPerDim * m_cellsPerDim * m_cellsPerDim);
	int nextPointIndex = 0;
	PointData * allPoints = m_allPoints;

	for(int i = 0; i < numPoints; ++i)
	{
		ToGridNoClamp(points[i], ix, iy, iz);

		if(ix >= 0 && ix < maxIndex &&
			iy >= 0 && iy < maxIndex &&
			iz >= 0 && iz < maxIndex)
		{
			int cellIdx = ix + yDimMult * iy + zDimMult * iz;
			//ASSERT(cellIdx >= 0 && cellIdx < (m_cellsPerDim * m_cellsPerDim * m_cellsPerDim));
			Cell* cell = m_cells[cellIdx];
			//ASSERT(cell);
			if(cell)
			{
				if(cell->m_points == 0)
				{
					cell->m_points = new int[cell->m_numPoints];
					curCounts[cellIdx] = 0;
				}
				//ASSERT(curCounts[cellIdx] < cell->m_numPoints);
				//ASSERT(nextPointIndex < numPointsInserted);

				allPoints[nextPointIndex] = points[i];
				m_pointsAABB.Extend(points[i]);
				cell->m_points[curCounts[cellIdx]++] = nextPointIndex;
				++nextPointIndex;
			}	
		}
	}

	delete[] curCounts;
}
	
int SparsePointGrid::ClosestPointToSplit(SplitDir dir, const AABB& boundsToSplit)
{
	int zDimMult = m_cellsPerDim * m_cellsPerDim;
	int yDimMult = m_cellsPerDim;

	Vec3 middle = boundsToSplit.m_min * 0.5f + boundsToSplit.m_max * 0.5f;
	int ix, iy, iz;
	ToGrid(middle, ix, iy, iz);
	int closestPointIndex = -1;
	float closestPointDist = FLT_MAX;
	float searchDistMax = m_cellDim;

	Plane splitPlane;
	MakeSplitPlane(splitPlane, dir, boundsToSplit);

	int iSearchedStart[3];
	int iSearchedEnd[3];

	for(int i = 0; i < 3; ++i)
	{
		iSearchedStart[i] = INT_MAX;
		iSearchedEnd[i] = INT_MIN;
	}
		
	int iStart[3];
	int iEnd[3];

	//ASSERT(dir >= 0 && dir < 3);
	for(int i = 0; i < 3; ++i)
	{
		if(i != dir)
		{
			ToGrid(boundsToSplit.m_min[i], iStart[i]);
			ToGrid(boundsToSplit.m_max[i], iEnd[i]);
		}
	}
			
	while(closestPointIndex == -1)
	{
		ToGrid(middle[dir] - searchDistMax, iStart[dir]);
		ToGrid(middle[dir] + searchDistMax, iEnd[dir]);
		
		int dim;
		for(dim = 0; dim < 3; ++dim)
		{
			if(iEnd[dim] > iSearchedEnd[dim] || iStart[dim] < iSearchedStart[dim])
				break;
		}
		if(dim == 3)
			break; // We're entirely in already searched cells, so we can quit.

		for(int iz = iStart[2]; iz <= iEnd[2]; ++iz)
		{
			for(int iy = iStart[1]; iy <= iEnd[1]; ++iy)
			{
				for(int ix = iStart[0]; ix <= iEnd[0]; ++ix)
				{
					if(ix >= iSearchedStart[0] && ix <= iSearchedEnd[0] &&
						iy >= iSearchedStart[1] && iy <= iSearchedEnd[1] &&
						iz >= iSearchedStart[2] && iz <= iSearchedEnd[2])
						continue;

					int cellIdx = ix + yDimMult * iy + zDimMult * iz;
					Cell* cell = m_cells[cellIdx];
					if(cell)
						FindClosestPointInCell(cell, splitPlane, closestPointIndex, closestPointDist);
				}
			}
		}

		for(int i = 0; i < 3; ++i)
		{
			iSearchedStart[i] = Min(iSearchedStart[i], iStart[i]);
			iSearchedEnd[i] = Max(iSearchedEnd[i], iEnd[i]);
		}

		searchDistMax += 2.f * m_cellDim ;
	}
	return closestPointIndex;
}
	
void SparsePointGrid::FindClosestPointInCell(const Cell* cell, 
	const Plane& plane, int& closestPointIndex, float &closestPointDist)
{
	int * points = cell->m_points;
	for(int i = 0, c = cell->m_numPoints; i < c; ++i)
	{	
		if(m_allPoints[points[i]].m_faceCount == 0)
			continue;
		float dist = DistToPlane(m_allPoints[points[i]].m_pos, plane);
		if(dist < closestPointDist)
		{
			closestPointDist = dist;
			closestPointIndex = points[i];
		}
	}
}
	
void SparsePointGrid::FindClosestPointInCellAbovePlane(const Cell* cell, 
	Vec3_arg fromPos, const Plane& plane, float fromPlaneDist,
	int& closestPointIndex, float &closestPointDistSq)
{
	int * points = cell->m_points;
	for(int i = 0, c = cell->m_numPoints; i < c; ++i)
	{
		if(m_allPoints[points[i]].m_faceCount == 0)
			continue;
		Vec3 pos = m_allPoints[points[i]].m_pos;
		float posDist = dot(pos, plane.m_normal) - plane.m_d ;
		if((fromPlaneDist < 0.f) != (posDist < 0.f))
		{
			float distSq = magnitude_squared(pos - fromPos);
			if(distSq < closestPointDistSq)
			{
				closestPointDistSq = distSq;
				closestPointIndex = points[i];
			}
		}
	}
}

int SparsePointGrid::NearestNeighborAcrossPlane(int from, const Plane& plane)
{
	Vec3 fromPos = GetPos(from);

	float searchDistMax = m_cellDim;

	int closestPointIndex = -1;
	float closestPointDistSq = FLT_MAX;

	float fromPlaneDist = dot(plane.m_normal, fromPos) - plane.m_d;
	
	int iSearchedStart[3];
	int iSearchedEnd[3];

	for(int i = 0; i < 3; ++i)
	{
		iSearchedStart[i] = INT_MAX;
		iSearchedEnd[i] = INT_MIN;
	}

	while(closestPointIndex == -1)
	{
		int iStart[3];
		int iEnd[3];
		ToGrid(fromPos - Vec3(searchDistMax, searchDistMax, searchDistMax), iStart[0], iStart[1], iStart[2]);
		ToGrid(fromPos + Vec3(searchDistMax, searchDistMax, searchDistMax), iEnd[0], iEnd[1], iEnd[2]);
		
		int dim;
		for(dim = 0; dim < 3; ++dim)
		{
			if(iEnd[dim] > iSearchedEnd[dim] || iStart[dim] < iSearchedStart[dim])
				break;
		}
		if(dim == 3)
			break; // We're entirely in already searched cells, so we can quit.

		int zDimOff = iStart[2] * m_cellsPerDim * m_cellsPerDim;
		for(int iz = iStart[2]; iz <= iEnd[2]; ++iz)
		{
			int yDimOff = iStart[1] * m_cellsPerDim;
			for(int iy = iStart[1]; iy <= iEnd[1]; ++iy)
			{
				for(int ix = iStart[0]; ix <= iEnd[0]; ++ix)
				{
					if(ix >= iSearchedStart[0] && ix <= iSearchedEnd[0] &&
						iy >= iSearchedStart[1] && iy <= iSearchedEnd[1] &&
						iz >= iSearchedStart[2] && iz <= iSearchedEnd[2])
						continue;

					int cellIdx = ix + yDimOff + zDimOff;
					//ASSERT(cellIdx >= 0 && cellIdx < (m_cellsPerDim * m_cellsPerDim * m_cellsPerDim));
					Cell* cell = m_cells[cellIdx];
					if(cell)
					{
						AABB cellBounds;
						cellBounds.m_min = Vec3(ix * m_cellDim + m_minGridDim, 
								iy * m_cellDim + m_minGridDim, 
								iz * m_cellDim + m_minGridDim);
						cellBounds.m_max = cellBounds.m_min + Vec3(m_cellDim, m_cellDim, m_cellDim);
						
						float distToClosestSq = DistSqAABBToPoint(cellBounds, fromPos);
						if(distToClosestSq < closestPointDistSq &&
							AABBAbovePlane(cellBounds, plane))
						{	
							//DebugDrawAABB(cellBounds);
							FindClosestPointInCellAbovePlane(cell, fromPos, plane, fromPlaneDist,
								closestPointIndex, 
								closestPointDistSq);
						}
					}
				}
				yDimOff += m_cellsPerDim;
			}
			zDimOff += m_cellsPerDim * m_cellsPerDim;
		}
		
		for(int i = 0; i < 3; ++i)
		{
			iSearchedStart[i] = Min(iSearchedStart[i], iStart[i]);
			iSearchedEnd[i] = Max(iSearchedEnd[i], iEnd[i]);
		}
		
		searchDistMax +=m_cellDim ;
	}

	return closestPointIndex;	
}

int SparsePointGrid::PointWithMinCircumcircle(int v0, int v1)
{
	Vec3 v0Pos = GetPos(v0);
	Vec3 v1Pos = GetPos(v1);
	Vec3 center = 0.5f * v0Pos + 0.5f * v1Pos;

	int bestRadiusIndex = -1;
	float bestRadius = FLT_MAX;
	float bestRadiusSq = FLT_MAX;		
	float searchDistMax = magnitude(center - v0Pos);

	int iSearchedStart[3], iSearchedEnd[3];
	for(int i = 0; i < 3; ++i)
	{
		iSearchedStart[i] = INT_MAX;
		iSearchedEnd[i] = INT_MIN;
	}

	while(true)
	{
		int iStart[3];
		int iEnd[3];

		ToGrid(center - Vec3(searchDistMax, searchDistMax, searchDistMax), iStart[0], iStart[1], iStart[2]);
		ToGrid(center + Vec3(searchDistMax, searchDistMax, searchDistMax), iEnd[0], iEnd[1], iEnd[2]);

		int dim;
		for(dim = 0; dim < 3; ++dim)
		{
			if(iEnd[dim] > iSearchedEnd[dim] || iStart[dim] < iSearchedStart[dim])
				break;
		}
		if(dim == 3)
			break; // We're entirely in already searched cells, so we can quit.

		int zDimOff = iStart[2] * m_cellsPerDim * m_cellsPerDim;
		for(int iz = iStart[2]; iz <= iEnd[2]; ++iz, zDimOff += m_cellsPerDim * m_cellsPerDim)
		{
			int yDimOff = iStart[1] * m_cellsPerDim;
			for(int iy = iStart[1]; iy <= iEnd[1]; ++iy, yDimOff += m_cellsPerDim)
			{
				for(int ix = iStart[0]; ix <= iEnd[0]; ++ix)
				{
					if(ix >= iSearchedStart[0] && ix <= iSearchedEnd[0] &&
						iy >= iSearchedStart[1] && iy <= iSearchedEnd[1] &&
						iz >= iSearchedStart[2] && iz <= iSearchedEnd[2])
						continue;

					int cellIdx = ix + yDimOff + zDimOff;
					//ASSERT(cellIdx >= 0 && cellIdx < (m_cellsPerDim * m_cellsPerDim * m_cellsPerDim));
					Cell* cell = m_cells[cellIdx];
					if(cell)
					{
						AABB cellBounds;
						cellBounds.m_min = Vec3(ix * m_cellDim + m_minGridDim, 
								iy * m_cellDim + m_minGridDim, 
								iz * m_cellDim + m_minGridDim);
						cellBounds.m_max = cellBounds.m_min + Vec3(m_cellDim, m_cellDim, m_cellDim);

						float distToClosestSq = DistSqAABBToPoint(cellBounds, center);
						if( distToClosestSq < bestRadiusSq)
						{	
							//DebugDrawAABB(cellBounds);
							int * points = cell->m_points;
							for(int i = 0, c = cell->m_numPoints; i < c; ++i)
							{
								if(m_allPoints[points[i]].m_faceCount == 0)
									continue;
								float radiusSq;
								Vec3 testCenter;
								Vec3 pos = m_allPoints[points[i]].m_pos;
								float distToPointSq = magnitude_squared(center - pos);
								if( distToPointSq < bestRadiusSq &&
									ComputeCircumcircle(v0Pos, v1Pos, pos, testCenter, radiusSq) && 
									radiusSq < bestRadiusSq)
								{
									bestRadiusIndex = points[i];
									bestRadiusSq = radiusSq;
									bestRadius = sqrtf(radiusSq);
									center = testCenter;
								}
							}
						}
					}
				}
			}
		}

		for(int i = 0; i < 3; ++i)
		{
			iSearchedStart[i] = Min(iSearchedStart[i], iStart[i]);
			iSearchedEnd[i] = Max(iSearchedEnd[i], iEnd[i]);
		}
	
		if(bestRadiusIndex == -1)
			searchDistMax += m_cellDim;
		else
			searchDistMax = bestRadius;
	}

	//DebugDrawSphere(center, bestRadius, 0.3f, 0.3f, 0.5f, 0.3f);

	return bestRadiusIndex;	
}

int SparsePointGrid::PointWithMinCircumsphere(const AABB& bounds, int v0, int v1, int v2, int flags)
{
	DVec3 v0Pos = GetPos(v0);
	DVec3 v1Pos = GetPos(v1);
	DVec3 v2Pos = GetPos(v2);
	double inv_w = 1/3.;
	DVec3 center = inv_w * v0Pos + inv_w * v1Pos + inv_w * v2Pos;

	DPlane trianglePlane;
	trianglePlane.m_normal = normalize(cross(v1Pos - v0Pos, v2Pos - v0Pos));
	trianglePlane.m_d = dot(v0Pos, trianglePlane.m_normal);

	if( !(flags & CONSTRAINT_POINT_ABOVE) )
	{
		trianglePlane.m_normal = Vec3(0,0,0);
		trianglePlane.m_d = -1.f;
	}

	int bestRadiusIndex = -1;
	double bestRadiusSq = DBL_MAX;		
	double searchDistMax = 2. * magnitude(center - v0Pos);

	int iSearchedStart[3], iSearchedEnd[3];
	for(int i = 0; i < 3; ++i)
	{
		iSearchedStart[i] = INT_MAX;
		iSearchedEnd[i] = INT_MIN;
	}

	while(true)
	{
		int iStart[3];
		int iEnd[3];

		ToGrid(bounds, center - DVec3(searchDistMax, searchDistMax, searchDistMax), iStart[0], iStart[1], iStart[2]);
		ToGrid(bounds, center + DVec3(searchDistMax, searchDistMax, searchDistMax), iEnd[0], iEnd[1], iEnd[2]);

		int dim;
		for(dim = 0; dim < 3; ++dim)
		{
			if(iEnd[dim] > iSearchedEnd[dim] || iStart[dim] < iSearchedStart[dim])
				break;
		}
		if(dim == 3)
			break; // We're entirely in already searched cells, so we can quit.

		int zDimOff = iStart[2] * m_cellsPerDim * m_cellsPerDim;
		for(int iz = iStart[2]; iz <= iEnd[2]; ++iz, zDimOff += m_cellsPerDim * m_cellsPerDim)
		{
			int yDimOff = iStart[1] * m_cellsPerDim;
			for(int iy = iStart[1]; iy <= iEnd[1]; ++iy, yDimOff += m_cellsPerDim)
			{
				for(int ix = iStart[0]; ix <= iEnd[0]; ++ix)
				{
					if(ix >= iSearchedStart[0] && ix <= iSearchedEnd[0] &&
						iy >= iSearchedStart[1] && iy <= iSearchedEnd[1] &&
						iz >= iSearchedStart[2] && iz <= iSearchedEnd[2])
						continue;

					int cellIdx = ix + yDimOff + zDimOff;
					//ASSERT(cellIdx >= 0 && cellIdx < (m_cellsPerDim * m_cellsPerDim * m_cellsPerDim));
					Cell* cell = m_cells[cellIdx];
					if(cell)
					{
						AABB cellBounds;
						cellBounds.m_min = Vec3(ix * m_cellDim + m_minGridDim, 
								iy * m_cellDim + m_minGridDim, 
								iz * m_cellDim + m_minGridDim);
						cellBounds.m_max = cellBounds.m_min + Vec3(m_cellDim, m_cellDim, m_cellDim);

						if( (flags & CONSTRAINT_POINT_ABOVE) && !AABBAbovePlane(cellBounds, trianglePlane))
							continue;

						float distToClosestSq = DistSqAABBToPoint(cellBounds, center);
						if(distToClosestSq < bestRadiusSq)
						{	
							//DebugDrawAABB(cellBounds);
							int * points = cell->m_points;
							for(int i = 0, c = cell->m_numPoints; i < c; ++i)
							{
								if(m_allPoints[points[i]].m_faceCount == 0)
									continue;
								DVec3 pos = m_allPoints[points[i]].m_pos;
								if(!AABBContains(bounds, pos))
									continue;
								double radiusSq = DBL_MAX;
								DVec3 testCenter;
								double distToPointSq = magnitude_squared(center - pos);
								double pointDist = dot(trianglePlane.m_normal, pos) - trianglePlane.m_d;
								
								if(	pointDist > 0.f &&
									distToPointSq < bestRadiusSq &&
									ComputeCircumsphere(v0Pos, v1Pos, v2Pos, pos, testCenter, radiusSq) &&
									radiusSq < bestRadiusSq)
								{
									if( flags & CONSTRAINT_DELAUNAY )
									{
										int exceptList[4] = { v0, v1, v2, points[i] };
										if(HasPointsInRadiusSq(testCenter, radiusSq, exceptList, 4))
											continue;
									}
									bestRadiusIndex = points[i];
									bestRadiusSq = radiusSq;
									center = testCenter;
								}
							}
						}
					}
				}
			}
		}

		for(int i = 0; i < 3; ++i)
		{
			iSearchedStart[i] = Min(iSearchedStart[i], iStart[i]);
			iSearchedEnd[i] = Max(iSearchedEnd[i], iEnd[i]);
		}
	
		if(bestRadiusIndex == -1)
			searchDistMax += m_cellDim;
		else 
			searchDistMax = sqrtf(bestRadiusSq);
	}

	//DebugDrawLine(inv_w * v0Pos + inv_w * v1Pos + inv_w * v2Pos, center, 1.f, 1.f, 0.f, 1.f);
	//DebugDrawVector(inv_w * v0Pos + inv_w * v1Pos + inv_w * v2Pos, 3.0 * normalize(trianglePlane.m_normal),
	//	1.f, 0.f, 1.f, 1.f);
	//DebugDrawSphere(center, sqrtf(bestRadiusSq), 0.3f, 0.3f, 0.5f, 0.3f);

	/*	if(bestRadiusSq > 4.f * magnitude_squared(0.5f * bounds.m_max - 0.5f * bounds.m_min))
	{
		return -1;
	}*/

	return bestRadiusIndex;	
}

void SparsePointGrid::FindPointsInRadius(Vec3_arg from, float radius, std::vector<int> &outPoints)
{
	outPoints.clear();
	float radiusSq = radius * radius;
	int iStart[3];
	int iEnd[3];

	ToGrid(from - Vec3(radius, radius, radius), iStart[0], iStart[1], iStart[2]);
	ToGrid(from + Vec3(radius, radius, radius), iEnd[0], iEnd[1], iEnd[2]);
	int zDimOff = iStart[2] * m_cellsPerDim * m_cellsPerDim;
	for(int iz = iStart[2]; iz <= iEnd[2]; ++iz, zDimOff += m_cellsPerDim * m_cellsPerDim)
	{
		int yDimOff = iStart[1] * m_cellsPerDim;
		for(int iy = iStart[1]; iy <= iEnd[1]; ++iy, yDimOff += m_cellsPerDim)
		{
			for(int ix = iStart[0]; ix <= iEnd[0]; ++ix)
			{
				int cellIdx = ix + yDimOff + zDimOff;
				//ASSERT(cellIdx >= 0 && cellIdx < (m_cellsPerDim * m_cellsPerDim * m_cellsPerDim));
				Cell* cell = m_cells[cellIdx];
				if(cell)
				{
					AABB cellBounds;
					cellBounds.m_min = Vec3(ix * m_cellDim + m_minGridDim, 
							iy * m_cellDim + m_minGridDim, 
							iz * m_cellDim + m_minGridDim);
					cellBounds.m_max = cellBounds.m_min + Vec3(m_cellDim, m_cellDim, m_cellDim);
					float distToClosestSq = DistSqAABBToPoint(cellBounds, from);
					if( distToClosestSq < radiusSq )
					{
						int * points = cell->m_points;
						for(int i = 0, c = cell->m_numPoints; i < c; ++i)
						{
							if(m_allPoints[points[i]].m_faceCount == 0)
								continue;
							Vec3 pos = m_allPoints[points[i]].m_pos;

							if(magnitude_squared(pos - from) < radiusSq)
								outPoints.push_back(points[i]);
						}
					}
				}
			}
		}
	}
}

bool SparsePointGrid::HasPointsInRadiusSq(Vec3_arg from, float radiusSq, const int *except, int numExcept)
{
	int iStart[3];
	int iEnd[3];

	float radius = sqrtf(radiusSq);
	ToGrid(from - Vec3(radius, radius, radius), iStart[0], iStart[1], iStart[2]);
	ToGrid(from + Vec3(radius, radius, radius), iEnd[0], iEnd[1], iEnd[2]);
	int zDimOff = iStart[2] * m_cellsPerDim * m_cellsPerDim;
	for(int iz = iStart[2]; iz <= iEnd[2]; ++iz, zDimOff += m_cellsPerDim * m_cellsPerDim)
	{
		int yDimOff = iStart[1] * m_cellsPerDim;
		for(int iy = iStart[1]; iy <= iEnd[1]; ++iy, yDimOff += m_cellsPerDim)
		{
			for(int ix = iStart[0]; ix <= iEnd[0]; ++ix)
			{
				int cellIdx = ix + yDimOff + zDimOff;
				//ASSERT(cellIdx >= 0 && cellIdx < (m_cellsPerDim * m_cellsPerDim * m_cellsPerDim));
				Cell* cell = m_cells[cellIdx];
				if(cell)
				{
					AABB cellBounds;
					cellBounds.m_min = Vec3(ix * m_cellDim + m_minGridDim, 
							iy * m_cellDim + m_minGridDim, 
							iz * m_cellDim + m_minGridDim);
					cellBounds.m_max = cellBounds.m_min + Vec3(m_cellDim, m_cellDim, m_cellDim);
					float distToClosestSq = DistSqAABBToPoint(cellBounds, from);
					if( distToClosestSq < radiusSq )
					{
						int * points = cell->m_points;
						for(int i = 0, c = cell->m_numPoints; i < c; ++i)
						{
							if(m_allPoints[points[i]].m_faceCount == 0)
								continue;
							int j = 0;
							for(j = 0; j < numExcept; ++j)
							{
								if(points[i] == except[j])
									break;
							}
							if(j != numExcept) 
								continue;

							Vec3 pos = m_allPoints[points[i]].m_pos;
							if(magnitude_squared(pos - from) < radiusSq)
							{
								return true;
							}
						}
					}
				}
			}
		}
	}
	return false;
}

Vec3 SparsePointGrid::GetPos(int pointIdx) const
{
	//ASSERT(pointIdx >= 0 && pointIdx < m_numPoints);
	return m_allPoints[pointIdx].m_pos;
}

void SparsePointGrid::AddRef(int pointIdx)
{
	//ASSERT(pointIdx >= 0 && pointIdx < m_numPoints);
	PointData& data = m_allPoints[pointIdx];
	//ASSERT(data.m_faceCount != 0);
	if(data.m_faceCount < 0) {
		data.m_faceCount = 1;
	} else {
		++data.m_faceCount;
	}
	printf("idx %d ++ref = %d\n", pointIdx, data.m_faceCount);
}

void SparsePointGrid::SubRef(int pointIdx)
{
	//ASSERT(pointIdx >= 0 && pointIdx < m_numPoints);
	PointData& data = m_allPoints[pointIdx];
	if(data.m_faceCount > 0)
	{
		--data.m_faceCount;
	}
	printf("idx %d --ref = %d\n", pointIdx, data.m_faceCount);
}

