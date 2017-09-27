//
//  CompFab.cpp
//  voxelizer
//
//
//

#include "../include/CompFab.h"


using namespace CompFab;





CompFab::Vec3Struct::Vec3Struct()
{
    m_x = m_y = m_z = 0.0;
}

CompFab::Vec3Struct::Vec3Struct(double x, double y, double z)
{
    m_x = x;
    m_y = y;
    m_z = z;
}

void CompFab::Vec3Struct::normalize() {
    
    double magnitude = sqrt(m_x*m_x+m_y*m_y+m_z*m_z);
    
    if(magnitude > EPSILON)
    {
        m_x /= magnitude;
        m_y /= magnitude;
        m_z /= magnitude;
    }
}

//Data Types
CompFab::Vec3iStruct::Vec3iStruct()
{
    m_x = m_y = m_z = 0.0;
}

CompFab::Vec3iStruct::Vec3iStruct(double x, double y, double z)
{
    m_x = x;
    m_y = y;
    m_z = z;
}

CompFab::Vec2fStruct::Vec2fStruct()
{
    m_x = m_y = 0.0;
}

CompFab::Vec2fStruct::Vec2fStruct(double x, double y)
{
    m_x = x;
    m_y = y;
}

CompFab::RayStruct::RayStruct()
{
    m_origin[0] = m_origin[1] = m_origin[2] = 0.0;
    m_direction[0] = 1.0;
    m_direction[1] = m_direction[2] = 0.0;
}

CompFab::RayStruct::RayStruct(Vec3 &origin, Vec3 &direction)
{
    m_origin = origin;
    m_direction = direction;
}

CompFab::TriangleStruct::TriangleStruct(Vec3 &v1, Vec3 &v2,Vec3 &v3)
{
    m_v1 = v1;
    m_v2 = v2;
    m_v3 = v3;
}

CompFab::TriangleStruct::TriangleStruct()
{
}

bool CompFab::operator==(const Vec3& lhs, const Vec3& rhs)
{
    return (lhs.m_x == rhs.m_x) &&
            (lhs.m_y == rhs.m_y) &&
            (lhs.m_z == rhs.m_z);
}

bool CompFab::operator==(const Triangle& lhs, const Triangle& rhs)
{
    return (lhs.m_v1 == rhs.m_v1) &&
            (lhs.m_v2 == rhs.m_v2) &&
            (lhs.m_v3 == rhs.m_v3);
}

std::ostream& CompFab::operator <<(std::ostream& os, const Vec3& v)
{
    return os << "<" << v.m_x << "," << v.m_y << "," << v.m_z << ">";
}

CompFab::Vec3 CompFab::operator-(const Vec3 &v1, const Vec3 &v2)
{
    Vec3 v3;
    v3[0] = v1[0] - v2[0];
    v3[1] = v1[1] - v2[1];
    v3[2] = v1[2] - v2[2];

    return v3;
}

CompFab::Vec3 CompFab::operator+(const Vec3 &v1, const Vec3 &v2)
{
    Vec3 v3;
    v3[0] = v1[0] + v2[0];
    v3[1] = v1[1] + v2[1];
    v3[2] = v1[2] + v2[2];
    
    return v3;
}


//Cross Product
Vec3 CompFab::operator%(const Vec3 &v1, const Vec3 &v2)
{
    Vec3 v3;
    v3[0] = v1[1]*v2[2] - v1[2]*v2[1];
    v3[1] = v1[2]*v2[0] - v1[0]*v2[2];
    v3[2] = v1[0]*v2[1] - v1[1]*v2[0];

    return v3;
}

//Dot Product
double CompFab::operator*(const Vec3 &v1, const Vec3 &v2)
{
    return v1.m_x*v2.m_x + v1.m_y*v2.m_y+v1.m_z*v2.m_z;
}

Vec3 CompFab::operator*(const Vec3 &v1, const double s)
{
    return Vec3(v1.m_x * s, v1.m_y * s, v1.m_z * s);
}


//Grid structure for Voxels
CompFab::VoxelGridStruct::VoxelGridStruct(Vec3 lowerLeft, unsigned int dimX, unsigned int dimY, unsigned int dimZ, double spacing)
{
    m_lowerLeft = lowerLeft;
    m_dimX = dimX;
    m_dimY = dimY;
    m_dimZ = dimZ;
    m_size = dimX*dimY*dimZ;
    m_spacing = spacing;
    
    //Allocate Memory
    m_insideArray = new bool[m_size];
    
    for(unsigned int ii=0; ii<m_size; ++ii)
    {
        m_insideArray[ii] = false;
    }
    
}

CompFab::VoxelGridStruct::~VoxelGridStruct()
{
    delete[] m_insideArray;
}





