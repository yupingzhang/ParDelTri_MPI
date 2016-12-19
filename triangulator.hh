#pragma once

#include <list>
#include <vector>
#include "sparsegrid.hh"

class TriangleHashList
{
	struct HashValue
	{
		HashValue() : v0(-1), v1(-1), v2(-1), frontFacing(false),
			next(0), nextInserted(0), prevInserted(0) {}

		int v0,v1,v2;
		bool frontFacing;
		HashValue *next;

		HashValue *nextInserted;
		HashValue *prevInserted;
	};

	std::vector< HashValue * > m_buckets;
	HashValue* m_lastInserted;

	struct StatsType
	{
		int numInsertions;
		int numCollisions;	
	} m_stats;
public:
	TriangleHashList(int numBuckets);
	~TriangleHashList();

	bool AddOrRemove(int v0, int v1, int v2);
	bool Remove(int v0, int v1, int v2);
	bool GetFront(int *verts);
	void GetAllFaces(std::vector<int>& outFaceIndices) const;	
	bool Empty() const { return m_lastInserted == 0; }

	void DebugDrawFaces(const SparsePointGrid * grid);
private:

	TriangleHashList(const TriangleHashList&);
	TriangleHashList& operator=(const TriangleHashList&);
	bool RemoveFromBuckets(int v0, int v1, int v2);
	bool RemoveFromBuckets(int idx, int v0, int v1, int v2);
};

class Triangulator
{
public:
	struct Tetrahedron
	{
		Tetrahedron(int v0, int v1, int v2, int v3)
			: v0(v0), v1(v1), v2(v2), v3(v3) {}
		int v0, v1, v2, v3;
	};
private:

	struct SplitNode
	{
		SplitNode(int dir, const AABB& bounds) 
			: m_splitDir(dir)
			, m_bounds(bounds)
			, m_activeFaces( Max(4,4*int(bounds.m_max[dir] - bounds.m_min[dir])))
            , rec_index(0)
		{
			m_children[0] = m_children[1] = 0;
		}

		int m_splitDir;
		AABB m_bounds;
		TriangleHashList m_activeFaces;
		SplitNode* m_children[2];
        
        // keep track of the recursive index
        int rec_index;
        
	private:
		SplitNode(const SplitNode&);
		SplitNode& operator=(const SplitNode&);
	};

	SparsePointGrid* m_grid;
	std::list< SplitNode* > m_nodesTodo;
	bool m_bStarted;
	std::vector< Tetrahedron > m_results;

public:
	Triangulator(SparsePointGrid* grid);
	~Triangulator();
	bool IsDone() const;
	void Step(int np, int pe_id, int& rec_call);

	int GetNumTetrahedrons() const { return m_results.size(); }
	const Tetrahedron& GetTetrahedron(int idx) const { return m_results[idx]; }

private:
	Triangulator(const Triangulator&);
	Triangulator& operator=(const Triangulator&);
	void StartWall(SplitNode* node);
	void ContinueWall(SplitNode* node);
	void AddSimplex(SplitNode* node, const Plane& plane, int v0, int v1, int v2, int v3);
	void InitActiveFaces(SplitNode* node);
    
    int count; // dividing times
};



// function on_path
// return true if the recursive call of index i is in the pass assigned to id
// recursive call index is bfs order
inline bool on_path(int np, int pe_id, int rec_call)
{
    // find the recursive call location in the tree
    int level = 0, index = 0;
    level = log2(rec_call);                  // level start at 0
    index = rec_call - pow(2.0, level);    // index start at 0
    
    // compute how many processes at level i execute on the same node
    float np_node = np / pow(2.0, level);
    
    int min = 0, max = 0;         // process id range on this node
    // min = int(np_node*float(index));
    //max = int(np_node*float(index + 1) + 0.5);
    min = int(floor(np_node*float(index)));
    max = int(ceil(np_node*float(index + 1)));

    if ((pe_id > min) && (pe_id <= max))
    {
      //printf("On_path yes: np: %d pe_id: %d rec_call: %d \n", np, pe_id, rec_call);
        return true;
    }
    else
    {
      //printf("On_path no: np: %d pe_id: %d rec_call: %d \n", np, pe_id, rec_call);
        return false;
    }
}


// function must_save
// input: np: number of processes, pe_id: process id start at 1, rec_call: recursive call index
// prevent replication in the output data
// only one of the m/2i pe_id actually returns its i-level tetrahedra wall
inline bool must_save(int np, int pe_id, int rec_call)
{
    // find the recursive call location in the tree
    int level = 0, index = 0;
    level = log2(rec_call);
    index = rec_call - pow(2.0, level);
    
    // compute how many processes at level i execute the same wall
    float np_node = np / pow(2.0, level);
    int min = int(np_node*float(index) + 1.0);
    
    if (pe_id == min)
    {
      //printf("Must_save YES: np: %d pe_id: %d rec_call: %d \n", np, pe_id, rec_call);
        return true;
    }
    else
    {
      //printf("Must_save No: np: %d pe_id: %d rec_call: %d \n", np, pe_id, rec_call);
        return false;
    }
}


