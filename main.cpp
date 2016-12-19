//
//  main.cpp
//  DelTri
//
//  Created by Yuping Zhang on 12/4/16.
//  Copyright © 2016 yzhang. All rights reserved.
//

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include "math.hh"
#include "vec3.hh"
#include "triangulator.hh"

#include <mpi.h>

namespace patch
{
  template < typename T > std::string to_string( const T& n )
  {
    std::ostringstream stm ;
    stm << n ;
    return stm.str() ;
  }
}

static int g_numPoints = 10000;					// default number of points
static unsigned int g_seed = 12345U;			// default random number seed
static int g_gridDims = 20;						// Default number of cells per dimension in the point grid
static Triangulator *s_triangulator;			// global triangulator
static SparsePointGrid *s_grid;					// global point grid
static const char* g_szPointFile;				// file to read points from
static const char* g_szMeshFile;				// file to write resulting tetrahedrons to

void GeneratePoints(std::vector<Vec3>& points)
{
    points.clear();
    points.resize(g_numPoints);
    
    unsigned int seed = g_seed;
    const float inv_randmax = 1.f/RAND_MAX;
    for(int i = 0, c = g_numPoints; i < c; ++i)
    {
        int ix = rand_r(&seed);
        int iy = rand_r(&seed);
        int iz = rand_r(&seed);
        
        Vec3 &pt = points[i];
        pt.x = ((ix * inv_randmax) * 2.f - 1.f) * 100.f;
        pt.y = ((iy * inv_randmax) * 2.f - 1.f) * 100.f;
        pt.z = ((iz * inv_randmax) * 2.f - 1.f) * 100.f;
    }
}

void LoadPoints(const char* szPointFile, std::vector<Vec3>& points)
{
    points.clear();
    //printf("Loading points from %s\n", szPointFile);
    
    FILE *fp = fopen(szPointFile, "rb");
    if(!fp)
    {
        printf("Failed to load point file %s\n", szPointFile);
        return;
    }
    
    Vec3 pt;
    char line[256];
    char *read;
    do
    {
        read = fgets(line, sizeof(line), fp);
        if(read)
        {
            int numRead = sscanf(read, "%f %f %f\n", &pt.x, &pt.y, &pt.z);
            if(numRead == 3)
            {
                points.push_back(pt);
            }
        }
    } while(read);
    fclose(fp);
    
    //printf("Read %d points.\n", int(points.size()));
}


void WriteVolumeMesh(const SparsePointGrid *grid, const Triangulator* triangulator, const char *szMeshFile)
{
    FILE* fp = fopen(szMeshFile, "wb");
    if(!fp)
    {
        printf("Failed to open %s for writing.\n", szMeshFile);
        return;
    }
    
    Vec3 tetPos[4];
    for(int i = 0, c = triangulator->GetNumTetrahedrons(); i < c; ++i)
    {
        const Triangulator::Tetrahedron& tet = triangulator->GetTetrahedron(i);
        tetPos[0] = grid->GetPos(tet.v0);
        tetPos[1] = grid->GetPos(tet.v1);
        tetPos[2] = grid->GetPos(tet.v2);
        tetPos[3] = grid->GetPos(tet.v3);
        
        for(int j = 0; j < 4; ++j)
            fprintf(fp, "%f %f %f ", tetPos[j].x, tetPos[j].y, tetPos[j].z);
        fprintf(fp, "\n");
    }
    fclose(fp);
    printf("Wrote %d tetrahedrons to %s.\n", triangulator->GetNumTetrahedrons(), szMeshFile);
}


void WriteMeshIndex(int id, int np, const SparsePointGrid *grid, const Triangulator* triangulator, const char *szMeshFile)
{
    std::string outfile = std::string(g_szMeshFile) + "_np" + patch::to_string(np) + "_id" + patch::to_string(id) \
    + ".txt";
    FILE* fp = fopen(outfile.c_str(), "wb");
    if(!fp)
    {
      printf("Failed to open %s for writing.\n", outfile.c_str());
        return;
    }

    fprintf(fp, "Process %d: Wrote %d tetrahedrons to %s.\n", id, triangulator->GetNumTetrahedrons(), outfile.c_str());    
    for(int i = 0, c = triangulator->GetNumTetrahedrons(); i < c; ++i)
    {
        const Triangulator::Tetrahedron& tet = triangulator->GetTetrahedron(i);
        fprintf(fp, "%d %d %d %d\n", tet.v0, tet.v1, tet.v2, tet.v3);
    }

    fclose(fp);
    
    int total = 0;
    if(id == 0)
    {
      printf("Process %d: Wrote %d tetrahedrons to %s.\n", id, triangulator->GetNumTetrahedrons(), outfile.c_str());  
      total += triangulator->GetNumTetrahedrons();
    }

    //printf("Wrote %d tetrahedrons to %s.\n", triangulator->GetNumTetrahedrons(), szMeshFile);
    char message[100];
    int dest = 0, tag = 0; 
    MPI_Status  status; 

    if (id != 0) {
      /* Create message */
      //sprintf(message, "Wrote %d tetrahedrons to %s.\n", triangulator->GetNumTetrahedrons(), szMeshFile);
      sprintf(message, "%d\n", triangulator->GetNumTetrahedrons());

      /* Use strlen+1 so that '\0' gets transmitted */                                   
      MPI_Send(message, strlen(message)+1, MPI_CHAR,
	       dest, tag, MPI_COMM_WORLD);
    } else { /* my_rank == 0 */

      for (int source = 1; source < np; source++) {
	      MPI_Recv(message, 100, MPI_CHAR, source, tag,
		  MPI_COMM_WORLD, &status);
	  if (message)
          {
	         total += atoi(message);
          }
	  std::string out = std::string(g_szMeshFile) + "_np" + patch::to_string(np) + "_id" + patch::to_string(source) + ".txt";
                                            
          printf("Process %d: Wrote %d tetrahedrons to %s.\n", source, atoi(message), out.c_str());
      }

      printf("Total tetrahedrons: %d\n", total);
    }
}


void AddOuterPoints(std::vector<Vec3>& points, const AABB& aabb)
{
    Vec3 bounds[] = { aabb.m_min + Vec3(EPSILON, EPSILON, EPSILON), 
        aabb.m_max - Vec3(EPSILON, EPSILON, EPSILON)};
    for(int i = 0; i < 8; ++i)
    {
        int zSel = (i >> 2) & 1;
        int ySel = (i >> 1) & 1;
        int xSel = i & 1;
        
        Vec3 pt(bounds[xSel][0], bounds[ySel][1], bounds[zSel][2]);
        points.push_back(pt);
    }
}


void initPointset()
{
    std::vector<Vec3> points;
    
    std::string pf = "";
    if (g_szPointFile) {
        pf = g_szPointFile;
    }
    if(pf != "")
    {
        LoadPoints(g_szPointFile, points);
    }
    else
    {
        //printf("Initializing points...\n");
        GeneratePoints(points);
    }
    
    if(points.empty())
    {
        exit(1);
    }
    //AddOuterPoints(points, AABB(Vec3(-128, -128, -128), Vec3(128, 128, 128)));
    
    //printf("building grid...\n");
   
    SparsePointGrid * grid = new SparsePointGrid(256.f, g_gridDims);
    grid->InsertPoints(&points[0], points.size());
    
    s_grid = grid;
    s_triangulator = new Triangulator(grid);
}


void triangulate(int np, int pe_id)
{
  //printf("Process %d: Triangulating...\n", pe_id);
    int rec_call = 1;               //recursive call index has to start with 1
    
    //int stepCount = 0;
    //char progress[] = { '/','-','\\', '|' };
    
    //int tmp = 0;
    while(!s_triangulator->IsDone())
    {
      //        if (tmp != rec_call)
      //  {
      //      printf("Process %d: Recursive call index: %d\n", pe_id, rec_call);
      //      tmp = rec_call;
      //  }

        s_triangulator->Step(np, pe_id, rec_call);
	    // ++stepCount;
          
        // if (stepCount % 100 == 0) {
        //     int idx = (stepCount / 100) % sizeof(progress);
        //     printf("%c (%d)\r", progress[idx], stepCount);
        // }
    }
    
    //printf("Process %d: Total step: %d\n", pe_id, stepCount);
}


int main(int argc, char *argv[])
{
    int id, np;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    char message[100];

    MPI_Status status;

    MPI_Init(&argc, &argv);

    /* Check processes: */
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    /* process command-line inputs: */
    if (argc < 3)
    {
        if (id == 0)
        {
            printf ("Usage: \"./triangulator input.txt output.txt\" \n");
        }
        MPI_Abort (MPI_COMM_WORLD, 1);
    }

    
    g_szPointFile = argv[1];
    g_szMeshFile = argv[2];

    // initialize pointset
    initPointset();    

    //take timing
    double start, finish, elapse, t_time = 0.0;
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    //triangulate();
    triangulate(np, id);

    MPI_Barrier(MPI_COMM_WORLD);
    finish = MPI_Wtime();
    
    elapse = finish - start;
    MPI_Reduce (&elapse, &t_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (id == 0)
    {
        printf("Elapsed time = %9.2f seconds\n", t_time);
    }

    // write the result to file
    //WriteVolumeMesh(s_grid, s_triangulator, g_szMeshFile);
    WriteMeshIndex(id, np, s_grid, s_triangulator, g_szMeshFile);
    
    MPI_Finalize();

    return 0;
}



