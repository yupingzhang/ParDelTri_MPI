#include <cstdio>
#include "math.hh"
#include "triangulator.hh"
//#include "debugdraw.hh"

Triangulator::Triangulator(SparsePointGrid *grid)
	: m_grid(grid)
	, m_nodesTodo()
	, m_bStarted(false)
	, m_results()
{
}

Triangulator::~Triangulator()
{
	for(std::list<SplitNode*>::iterator iter = m_nodesTodo.begin(),
		end = m_nodesTodo.end(); iter != end; ++iter)
	{
		delete *iter;
	}
}

bool Triangulator::IsDone() const
{
	return m_bStarted && m_nodesTodo.empty();
}

void Triangulator::Step(int np, int pe_id, int& rec_call)
{
	//ASSERT(m_grid);
    if (!m_grid) {
        return;
    }
    
	if(m_bStarted == false)
	{
		m_bStarted = true;
		SplitNode* node = new SplitNode(0, m_grid->GetAllPointsAABB());
		node->rec_index = 1;
		m_nodesTodo.push_front(node);
		StartWall(node);
	}
	else
	{
		SplitNode* curNode = 0;
		while(!m_nodesTodo.empty())
		{
			curNode = m_nodesTodo.front();
			rec_call = curNode->rec_index;

			if(curNode->m_activeFaces.Empty())
			{
 				int id = pe_id + 1;

                //make sure only one of the processes returns its i-level tetrahedra wall
                if (!must_save(np, id, rec_call)) {
                    m_results.clear();
                }

                //Debug
		//                if (curNode->m_children[0] == NULL)
                //{
                //	printf("No child on rec_call: %s\n", rec_call*2);
                //}
                //if (curNode->m_children[1] == NULL)
                //{
                //	printf("No child on rec_call: %s\n", rec_call*2 + 1);
	       	//}
                
                // parse PosHalf
                if (curNode->m_children[0] && on_path(np, id, rec_call*2))
                {
                    curNode->m_children[0]->rec_index = rec_call*2;
                    //printf("Add child 0, go to PosHalf, recursive call: %d\n", rec_call*2);
                    m_nodesTodo.push_back(curNode->m_children[0]);
                    //InitActiveFaces(curNode->m_children[0]);
                }
                // parse NegHalf
                if(curNode->m_children[1] && on_path(np, id, rec_call*2 + 1))
                {
                    curNode->m_children[1]->rec_index = rec_call*2 + 1;
                    //printf("Add child 1, go to NegHalf, recursive call: %d\n", rec_call*2+1);
                    m_nodesTodo.push_back(curNode->m_children[1]);
                    //InitActiveFaces(curNode->m_children[1]);
                }
                
                // delete current node
                m_nodesTodo.pop_front();
                delete curNode;
                curNode = 0;
                
                // start a new node
                if (!m_nodesTodo.empty()) {
                    InitActiveFaces(m_nodesTodo.front());
                }
            }
			else
				break;
		}

		if(curNode)
		{
			ContinueWall(curNode);
		}
	}
}
	
void Triangulator::StartWall(SplitNode* node)
{
	//DebugDrawAABB(node->m_bounds);
	// Make the first tetrahedron
	int v0 = m_grid->ClosestPointToSplit(SparsePointGrid::SplitDir(node->m_splitDir), node->m_bounds);
	if(v0 == -1)
	{
		printf("StartWall: Failed to find first point.\n");
		return;
	}
	//DebugDrawPoint(m_grid->GetPos(v0), 1.f, 0.f, 0.f);

	Plane splitPlane;
	MakeSplitPlane(splitPlane, node->m_splitDir, node->m_bounds);

	//DebugDrawPlane(node->m_bounds, splitPlane);
	int v1 = m_grid->NearestNeighborAcrossPlane(v0, splitPlane);
	if(v1 == -1)
	{
		printf("StartWall: Failed to find second point.\n");
		return;
	}
	//DebugDrawPoint(m_grid->GetPos(v1), 1.f, 0.f, 0.f);

//#ifdef DEBUG
//	float d0 = dot(splitPlane.m_normal, m_grid->GetPos(v0)) - splitPlane.m_d;
//	float d1 = dot(splitPlane.m_normal, m_grid->GetPos(v1)) - splitPlane.m_d;
//	ASSERT( (d0 < 0.f) != (d1 < 0.f) );
//#endif

	int v2 = m_grid->PointWithMinCircumcircle(v0, v1);
	if(v2 == -1)
	{
		printf("StartWall: Failed to find third point.\n");
		return;
	}
	//DebugDrawPoint(m_grid->GetPos(v2), 0.f, 1.f, 0.f);

	int v3 = m_grid->PointWithMinCircumsphere(node->m_bounds, v0, v1, v2, SparsePointGrid::CONSTRAINT_DELAUNAY);
	if(v3 == -1)
	{
		printf("StartWall: Failed to find last point of first simplex.\n");
		return;
	}
	//DebugDrawPoint(m_grid->GetPos(v3), 1.f, 0.f, 1.f);
	
	AddSimplex(node, splitPlane, v0, v1, v2, v3 );
	
	node->m_activeFaces.DebugDrawFaces(m_grid);
}
	
void Triangulator::AddSimplex(SplitNode* node, const Plane& plane, int v0, int v1, int v2, int v3)
{
	Vec3 v0Pos = m_grid->GetPos(v0);
	Vec3 v1Pos = m_grid->GetPos(v1);
	Vec3 v2Pos = m_grid->GetPos(v2);
	Vec3 v3Pos = m_grid->GetPos(v3);
	// Sort the first 3 vertices so they are counter clockwise from the outside of the tetrahedron
	Vec3 middle = 0.25f * v0Pos + 0.25f * v1Pos + 0.25f * v2Pos + 0.25f * v3Pos;

	Vec3 v0v1 = v1Pos - v0Pos;
	Vec3 v0v2 = v2Pos - v0Pos;
	Vec3 normal = cross(v0v1, v0v2);
	// if the middle is "above" the plane, then we need to switch the order of v1 and v2.
	if(dot(normal, middle) > dot(normal, v0Pos))
	{
		int temp = v1;
		v1 = v2;
		v2 = temp;

		Vec3 tempPos = v1Pos;
		v1Pos = v2Pos;
		v2Pos = tempPos;

		normal = -normal;
	}

	m_results.push_back(Tetrahedron(v0,v1,v2,v3));
	Vec3 triangleData[3*4] =
	{
		v0Pos, v1Pos, v2Pos,
		v1Pos, v0Pos, v3Pos,
		v0Pos, v2Pos, v3Pos,
		v1Pos, v3Pos, v2Pos,
	};

	int triangleIndices[3*4] =
	{
		v0, v1, v2,
		v1, v0, v3,
		v0, v2, v3,
		v1, v3, v2
	};

	int results[4];
	PlaneIntersectsTriangleList(plane, 4, triangleData, results);

	// at least one triangle should be intersecting the middle plane.
	//ASSERT( results[0] != results[1] || results[0] != results[2] || results[0] != results[3] ||
	//	results[1] != results[2] || results[1] != results[3] ||
	//	results[2] != results[3]);
    if ( !(results[0] != results[1] || results[0] != results[2] || results[0] != results[3] ||
         results[1] != results[2] || results[1] != results[3] ||
         results[2] != results[3]))
    {
        printf("Skip this simplex, at least one triangle should be intersecting the middle plane.\n");
        return;
    }

	for(int i = 0; i < 4; ++i)
	{
		int *indices = &triangleIndices[3*i];
		if(results[i] == TRI_PLANE_INTERSECT_ALL_ABOVE)
		{
			if(node->m_children[0] == 0)
			{
				AABB bounds = node->m_bounds;
				bounds.m_min[node->m_splitDir] = 0.5f * bounds.m_min[node->m_splitDir]
					+ 0.5f * bounds.m_max[node->m_splitDir];
				node->m_children[0] = new SplitNode((node->m_splitDir + 1) % 3, bounds);
				//m_nodesTodo.push_back(node->m_children[0]);
			}
			node->m_children[0]->m_activeFaces.AddOrRemove(indices[0], indices[1], indices[2]);
				
		}
		else if(results[i] == TRI_PLANE_INTERSECT_ALL_BELOW)
		{
			if(node->m_children[1] == 0)
			{
				AABB bounds = node->m_bounds;
				bounds.m_max[node->m_splitDir] = 0.5f * bounds.m_min[node->m_splitDir]
					+ 0.5f * bounds.m_max[node->m_splitDir];
				node->m_children[1] = new SplitNode((node->m_splitDir + 1) % 3, bounds);
				//m_nodesTodo.push_back(node->m_children[1]);
			}
			node->m_children[1]->m_activeFaces.AddOrRemove(indices[0], indices[1], indices[2]);
		} 
		else
		{
			node->m_activeFaces.AddOrRemove(indices[0], indices[1], indices[2]);
		}
	}
}

void Triangulator::InitActiveFaces(SplitNode* node)
{
	Plane splitPlane;
	MakeSplitPlane(splitPlane, node->m_splitDir, node->m_bounds);

	std::vector<int> allFaceIndices;
	node->m_activeFaces.GetAllFaces(allFaceIndices);
	
	std::vector<Vec3> allFaceData(allFaceIndices.size());
	for(int i = 0, c = allFaceIndices.size(); i < c; ++i)
		allFaceData[i] = m_grid->GetPos(allFaceIndices[i]);

	std::vector<int> allResults(allFaceIndices.size() / 3);
	PlaneIntersectsTriangleList(splitPlane, allResults.size(), &allFaceData[0], &allResults[0]);

	for(int iIndex = 0, iResult = 0, c = allFaceIndices.size(); iIndex < c; iIndex += 3, iResult += 1)
	{
		int *indices = &allFaceIndices[iIndex];	
		if(allResults[iResult] == TRI_PLANE_INTERSECT_ALL_ABOVE)
		{
			//VERIFY(node->m_activeFaces.Remove(indices[0], indices[1], indices[2]));
            //debug
            if (!node->m_activeFaces.Remove(indices[0], indices[1], indices[2])) {
                printf("Cannot remove this face from the active face list...\n");
                break;
            }
            
            if(node->m_children[0] == 0)
			{
				AABB bounds = node->m_bounds;
				bounds.m_min[node->m_splitDir] = 0.5f * bounds.m_min[node->m_splitDir]
					+ 0.5f * bounds.m_max[node->m_splitDir];
				node->m_children[0] = new SplitNode((node->m_splitDir + 1) % 3, bounds);
				//m_nodesTodo.push_back(node->m_children[0]);
			}
			node->m_children[0]->m_activeFaces.AddOrRemove(indices[0], indices[1], indices[2]);
				
		}
		else if(allResults[iResult] == TRI_PLANE_INTERSECT_ALL_BELOW)
		{
			//VERIFY(node->m_activeFaces.Remove(indices[0], indices[1], indices[2]));
            //debug
            if (!node->m_activeFaces.Remove(indices[0], indices[1], indices[2])) {
                printf("Cannot remove this face from the active face list...\n");
                break;
            }
            
			if(node->m_children[1] == 0)
			{
				AABB bounds = node->m_bounds;
				bounds.m_max[node->m_splitDir] = 0.5f * bounds.m_min[node->m_splitDir]
					+ 0.5f * bounds.m_max[node->m_splitDir];
				node->m_children[1] = new SplitNode((node->m_splitDir + 1) % 3, bounds);
				//m_nodesTodo.push_back(node->m_children[1]);
			}
			node->m_children[1]->m_activeFaces.AddOrRemove(indices[0], indices[1], indices[2]);
		} 
	}

}

void Triangulator::ContinueWall(SplitNode* node)
{
	//ASSERT(!node->m_activeFaces.Empty());
    //debug
    if (!node->m_activeFaces.Empty()) {
        
        // Build a simplex from the next active face.
        Plane splitPlane;
        MakeSplitPlane(splitPlane, node->m_splitDir, node->m_bounds);

        //DebugDrawAABB(node->m_bounds);
        //DebugDrawPlane(node->m_bounds, splitPlane);

        int verts[3];
        if(!node->m_activeFaces.GetFront(verts))
        {
            printf("Could not find next active face.\n");
            return;
        }

    //	for(int i = 0; i < 3; ++i)
    //	{
    //		DebugDrawPoint(m_grid->GetPos(verts[i]), 0.f, 0.f, 1.f);
    //	}
        
        int v3 = m_grid->PointWithMinCircumsphere(node->m_bounds, verts[0], verts[1], verts[2],
            SparsePointGrid::CONSTRAINT_POINT_ABOVE | SparsePointGrid::CONSTRAINT_DELAUNAY);
        if(v3 == -1)
        {
            node->m_activeFaces.Remove(verts[0], verts[1], verts[2]);
        }
        else
        {
            //DebugDrawPoint(m_grid->GetPos(v3), 1.f, 0.f, 1.f);
            AddSimplex(node, splitPlane, verts[0], verts[1], verts[2], v3 );
        }

    //	node->m_activeFaces.DebugDrawFaces(m_grid);
    }
}

////////////////////////////////////////////////////////////////////////////////
// TriangleHashList

TriangleHashList::TriangleHashList(int numBuckets)
	: m_buckets(numBuckets, static_cast<HashValue*>(0))
	, m_lastInserted(0)
	, m_stats()
{
	memset(&m_stats, 0, sizeof(m_stats));
}

TriangleHashList::~TriangleHashList()
{
	for(int i = 0, c = m_buckets.size(); i < c; ++i)
	{	
		HashValue* value = m_buckets[i];
		if(value)
		{
			HashValue* next = value->next;
			delete value;
			value = next;
		}
	}
}

static int MakeHash(int v0, int v1, int v2)
{
	// Totally arbitrary primes. Should probably verify that this is half-decent.
	int hash = (31 * v0) ^ (3023 * v1) ^ (1783 * v2);
	return hash;
}

bool TriangleHashList::AddOrRemove(int v0, int v1, int v2)
{
	int sorted[3];
	sorted[0] = v0;
	sorted[1] = v1;
	sorted[2] = v2;
	bool frontFacing = true;

	if(sorted[1] < sorted[0])
	{
		Swap(sorted[1], sorted[0]);
		frontFacing = !frontFacing;
	}
	if(sorted[2] < sorted[0])
	{
		Swap(sorted[2], sorted[0]);
		frontFacing = !frontFacing;
	}
	if(sorted[2] < sorted[1])
	{
		Swap(sorted[1], sorted[2]);
		frontFacing = !frontFacing;
	}

	//ASSERT(m_lastInserted == 0 || m_lastInserted->prevInserted == 0);

    if (m_lastInserted == 0 || m_lastInserted->prevInserted == 0) {
        

	int hash = MakeHash(sorted[0], sorted[1], sorted[2]);
	int idx = hash % m_buckets.size();

	if(RemoveFromBuckets(idx, sorted[0], sorted[1], sorted[2])) 
	{
		return false;
	}

	HashValue *value = new HashValue;
	value->v0 = sorted[0];
	value->v1 = sorted[1];
	value->v2 = sorted[2];
	value->frontFacing = frontFacing;
	value->nextInserted = m_lastInserted;
	if(m_lastInserted)
		m_lastInserted->prevInserted = value;
	m_lastInserted = value;
	++m_stats.numInsertions;
	if(m_buckets[idx])
	{
		++m_stats.numCollisions;
		value->next = m_buckets[idx];
	}
	m_buckets[idx] = value;
    
    }
    
	return true;
}

bool TriangleHashList::GetFront(int *verts)
{
	if(m_lastInserted)
	{
		HashValue *result = m_lastInserted;
		verts[0] = result->v0;
		verts[1] = result->v1;
		verts[2] = result->v2;
		if(!result->frontFacing)
			Swap(verts[1], verts[2]);
		return true;	
	}
	else return false;
}
	
void TriangleHashList::DebugDrawFaces(const SparsePointGrid * grid)
{
	HashValue * cur = m_lastInserted;
	while(cur)
	{	
		int v1 = cur->v1;
		int v2 = cur->v2;
		if(!cur->frontFacing)
			Swap(v1,v2);
		//Vec3 v0Pos = grid->GetPos(cur->v0);
		//Vec3 v1Pos = grid->GetPos(v1);
		//Vec3 v2Pos = grid->GetPos(v2);
		//DebugDrawTriangle(v0Pos, v1Pos, v2Pos, 0.3f, 0.3f, 0.5f, 0.30f);
		//DebugDrawVector( (v0Pos + v1Pos + v2Pos) / 3.f, 3.f * normalize(cross(v1Pos - v0Pos, v2Pos - v0Pos)),
		//	0.f, 1.f, 0.f, 1.f);

		cur = cur->nextInserted;
	}
}

bool TriangleHashList::Remove(int v0, int v1, int v2)
{
	int sorted[3];
	sorted[0] = v0;
	sorted[1] = v1;
	sorted[2] = v2;

	if(sorted[1] < sorted[0])
	{
		Swap(sorted[1], sorted[0]);
	}
	if(sorted[2] < sorted[0])
	{
		Swap(sorted[2], sorted[0]);
	}
	if(sorted[2] < sorted[1])
	{
		Swap(sorted[1], sorted[2]);
	}

	return RemoveFromBuckets(sorted[0], sorted[1], sorted[2]);
}

bool TriangleHashList::RemoveFromBuckets(int v0, int v1, int v2)
{
	int hash = MakeHash(v0, v1, v2);
	int idx = hash % m_buckets.size();
	return RemoveFromBuckets(idx, v0, v1, v2);
}

bool TriangleHashList::RemoveFromBuckets(int idx, int v0, int v1, int v2)
{
	HashValue* value = m_buckets[idx];
	HashValue* last = 0;
	
	while(value)
	{
		if(value->v0 == v0 &&
			value->v1 == v1 &&
			value->v2 == v2)
		{
			if(last) {
				last->next = value->next;
			} else {
				m_buckets[idx] = value->next;
			}

			if(value->prevInserted) value->prevInserted->nextInserted = value->nextInserted;
			if(value->nextInserted) value->nextInserted->prevInserted = value->prevInserted;

			if(value == m_lastInserted)
			{
				m_lastInserted = value->nextInserted;
				//ASSERT(m_lastInserted == 0 || m_lastInserted->prevInserted == 0);
			}

			delete value;
			return true;
		}

		last = value;
		value = value->next;
	}
	return false;
}

void TriangleHashList::GetAllFaces(std::vector<int>& outIndices) const
{
	HashValue* cur = m_lastInserted;
	while(cur)
	{
		int v0 = cur->v0;
		int v1 = cur->v1;
		int v2 = cur->v2;
		if(!cur->frontFacing)
			Swap(v1,v2);
		outIndices.push_back(v0);
		outIndices.push_back(v1);
		outIndices.push_back(v2);
		cur = cur->nextInserted;
	}
}

