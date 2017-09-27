//Computational Fabrication Assignment #1
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <chrono> 
#include "../include/CompFab.h"
#include "../include/Mesh.h"

/*
//Ray-Triangle Intersection
//Returns 1 if triangle and ray intersect, 0 otherwise
int rayTriangleIntersection(CompFab::Ray &ray, CompFab::Triangle &triangle)
{
    //Moved Implementation to Mesh.cpp to be used in KD tree
}
*/

//Triangle list (global)
typedef std::vector<CompFab::Triangle> TriangleList;

TriangleList g_triangleList;
Mesh *g_mesh;
KdNode *g_kdTree;
CompFab::VoxelGrid *g_voxelGrid;

//Number of intersections with surface made by a ray originating at voxel and cast in direction.
int numSurfaceIntersections(CompFab::Vec3 &voxelPos, CompFab::Vec3 &dir)
{
    
    unsigned int numHits = 0;
    
    /********* ASSIGNMENT *********/
    /* Check and return the number of times a ray cast in direction dir, 
     * from voxel center voxelPos intersects the surface */
    
    CompFab::Ray ray = CompFab::Ray(voxelPos, dir);
    // for(auto& tri: g_triangleList) {
    //     //numHits += rayTriangleIntersection(ray, tri);   
    // }

    numHits += g_kdTree->rayIntersection(ray, 0);    
    return numHits;
}

bool loadMesh(char *filename, unsigned int dim)
{
    g_triangleList.clear();
    
    g_mesh = new Mesh(filename, true);
    
    CompFab::Vec3 v1, v2, v3;

    //copy triangles to global list
    for(unsigned int tri =0; tri<g_mesh->t.size(); ++tri)
    {
        v1 = g_mesh->v[g_mesh->t[tri][0]];
        v2 = g_mesh->v[g_mesh->t[tri][1]];
        v3 = g_mesh->v[g_mesh->t[tri][2]];
        g_triangleList.push_back(CompFab::Triangle(v1,v2,v3));
    }

    //Create Voxel Grid
    CompFab::Vec3 bbMax, bbMin;
    BBox(*g_mesh, bbMin, bbMax);
    
    //Build Voxel Grid
    double bbX = bbMax[0] - bbMin[0];
    double bbY = bbMax[1] - bbMin[1];
    double bbZ = bbMax[2] - bbMin[2];
    double spacing;
    
    if(bbX > bbY && bbX > bbZ)
    {
        spacing = bbX/(double)(dim-2);
    } else if(bbY > bbX && bbY > bbZ) {
        spacing = bbY/(double)(dim-2);
    } else {
        spacing = bbZ/(double)(dim-2);
    }
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    
    g_voxelGrid = new CompFab::VoxelGrid(bbMin-hspacing, dim, dim, dim, spacing);
    
    return true;
   
}

void saveVoxelsToObj(const char * outfile)
{
 
    Mesh box;
    Mesh mout;
    int nx = g_voxelGrid->m_dimX;
    int ny = g_voxelGrid->m_dimY;
    int nz = g_voxelGrid->m_dimZ;
    double spacing = g_voxelGrid->m_spacing;
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    
    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {
                if(!g_voxelGrid->isInside(ii,jj,kk)){
                    continue;
                }
                CompFab::Vec3 coord(((double)ii)*spacing, ((double)jj)*spacing, ((double)kk)*spacing);
                CompFab::Vec3 box0 = coord - hspacing;
                CompFab::Vec3 box1 = coord + hspacing;
                makeCube(box, box0, box1);
                mout.append(box);
            }
        }
    }

    mout.save_obj(outfile);
}

int main(int argc, char **argv)
{
    unsigned int dim = 128; //dimension of voxel grid (e.g. 32x32x32)
    //Load OBJ
    if(argc < 3)
    {
        std::cout<<"Usage: Voxelizer InputMeshFilename OutputMeshFilename \n";
        return 0;
    }
    
    std::cout<<"Load Mesh : "<<argv[1]<<"\n";
    loadMesh(argv[1], dim);
    
    /********* ASSIGNMENT *********/
    /* Iterate over all voxels in g_voxelGrid and test whether they are inside our outside of the
     * surface defined by the triangles in g_triangleList */
    unsigned int num_samples = 3;
    double cutoff = 0.5; // Fraction of rays that need to be inside the object in order to fill the voxel


    CompFab::Vec3 voxelPos(0.0, 0.0, 0.0);
    double spacing = g_voxelGrid->m_spacing;


    auto start = std::chrono::high_resolution_clock::now();

    // Build the KD-tree
    std::cout << "Building KD-Tree..." << std::endl;
    g_kdTree = new KdNode();
    g_kdTree = g_kdTree->build(g_triangleList, 0);
    std::cout << "Done" << std::endl;

    // Fill the voxels
    std::cout << "Generating voxels..." << std::endl;
    for(int xi = 0; xi < g_voxelGrid->m_dimX; ++xi) {
        for(int yi = 0; yi < g_voxelGrid->m_dimY; ++yi) {
            for(int zi = 0; zi < g_voxelGrid->m_dimZ; ++zi) {

                CompFab::Vec3 offset(((double)xi)*spacing,
                                    ((double)yi)*spacing,
                                    ((double)zi)*spacing);
                voxelPos = g_voxelGrid->m_lowerLeft + offset;

                // Sample multiple directions and use a cutoff fraction to decide if it is filled
                int num_inside = 0;
                for(int si = 0; si < num_samples; si++) {
                    double theta = (double)rand()/(double)RAND_MAX * 2.0 * M_PI; // [0, 2pi]
                    double z = (double)rand()/(double)RAND_MAX * 2.0 - 1.0; // [-1, 1]
                    double r = sqrt(1.0 - z*z);
                    CompFab::Vec3 direction(r * cos(theta), r*sin(theta), z);

                    int n_intersections = numSurfaceIntersections(voxelPos, direction);
                    if(n_intersections % 2 != 0) {
                        num_inside++;
                    }
                }

                if((double)num_inside / (double)num_samples > cutoff) {
                    bool& insideArrayRef = g_voxelGrid->isInside(xi, yi, zi);
                    insideArrayRef = true;
                }
            }
        }   
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Done. Took: " << elapsed.count() << "s" << std::endl;
    
    //Write out voxel data as obj
    std::cout << "Saving obj..." << std::endl;
    saveVoxelsToObj(argv[2]);
    std::cout << "Done" << std::endl;
    
    delete g_voxelGrid;
    delete g_mesh;
    delete g_kdTree;
}