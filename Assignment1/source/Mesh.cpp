#include "../include/Mesh.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#ifdef _WIN32
#define NOMINMAX //Stop errors with std::max
#include <windows.h>
#endif
#include <stdio.h>
#include <cstdlib>
#include <utility>
#include <map>
#include <sstream>
#include <string.h>
//#include "util.h"

typedef double real_t;

///@brief is a point inside a box
bool ptInBox(const CompFab::Vec3 & mn,
    const CompFab::Vec3 mx, const CompFab::Vec3 & x)
{
  for(int dim = 0 ;dim<3;dim++){
    if(x[dim]<mn[dim] || x[dim] > mx[dim]){
        return false;
    }
  }
  return true;
}

void makeCube(Mesh & m, const CompFab::Vec3 & mn,
    const CompFab::Vec3 mx)
{
  CompFab::Vec3 ss = mx -mn;
  m=UNIT_CUBE;
  for(unsigned int ii = 0;ii<m.v.size();ii++){
    m.v[ii][0] = mn[0] + ss[0]*m.v[ii][0];
    m.v[ii][1] = mn[1] + ss[1]*m.v[ii][1];
    m.v[ii][2] = mn[2] + ss[2]*m.v[ii][2];

  }
}

void Mesh::append(const Mesh & m)
{
  unsigned int offset = v.size();
  unsigned int ot = t.size();
  v.insert(v.end(),m.v.begin(),m.v.end());
  t.insert(t.end(),m.t.begin(), m.t.end());
  for(unsigned int ii = ot;ii<t.size();ii++){
    for(int jj = 0 ;jj<3;jj++){
      t[ii][jj] += offset;
    }
  }
}

Mesh & Mesh::operator= (const Mesh& m)
{
  v = m.v;
  t = m.t;
  n = m.n;
  return *this;
}

///@brief cube [0,1]^3
CompFab::Vec3 CUBE_VERT[8]={
    CompFab::Vec3 (0, 0, 0),
    CompFab::Vec3 (1, 0, 0),
    CompFab::Vec3 (1, 1, 0),
    CompFab::Vec3 (0, 1, 0),
    CompFab::Vec3 (0, 0, 1),
    CompFab::Vec3 (1, 0, 1),
    CompFab::Vec3 (1, 1, 1),
    CompFab::Vec3 (0, 1, 1)
};

CompFab::Vec3i CUBE_TRIG[12]={CompFab::Vec3i(0,3,1),
CompFab::Vec3i(1, 3, 2),
CompFab::Vec3i(5, 4, 0),
CompFab::Vec3i(5, 0, 1),
CompFab::Vec3i(6, 5, 1),
CompFab:: Vec3i(1, 2, 6),
CompFab:: Vec3i(3, 6, 2),
CompFab:: Vec3i(3, 7, 6),
CompFab:: Vec3i(4, 3, 0),
CompFab:: Vec3i(4, 7, 3),
CompFab:: Vec3i(7, 4, 5),
CompFab:: Vec3i(7, 5, 6)};
Mesh UNIT_CUBE(CUBE_VERT,CUBE_TRIG);

Mesh::Mesh():v(0),t(0){}

Mesh::Mesh(const std::vector<CompFab::Vec3>&_v,
    const std::vector<CompFab::Vec3i>&_t):v(_v),t(_t)
{
  compute_norm();
}

Mesh::Mesh(const CompFab::Vec3 * _v,
  const CompFab::Vec3i * _t)
{
  v.assign(_v,_v+8);
  t.assign(_t,_t+12);
  
  compute_norm();
}

void Mesh::save(std::ostream & out, std::vector<CompFab::Vec3> * vert)
{
  std::string vTok("v");
  std::string fTok("f");
  std::string texTok("vt");
  char bslash='/';
  std::string tok;
  if(vert==0){
    vert = &v;
  }
  for(size_t ii=0;ii<vert->size();ii++){
    out<<vTok<<" "<<(*vert)[ii][0]<<" "<<(*vert)[ii][1]<<" "<<(*vert)[ii][2]<<"\n";
  }
  if(tex.size()>0){
    for(size_t ii=0;ii<tex.size();ii++){
      out<<texTok<<" "<<tex[ii][0]<<" "<<tex[ii][1]<<"\n";
    }
    for(size_t ii=0;ii<t.size();ii++){
      out<<fTok<<" "<<t[ii][0]+1<<bslash<<texId[ii][0]+1<<" "
      <<t[ii][1]+1<<bslash<<texId[ii][1]+1<<" "
      <<t[ii][2]+1<<bslash<<texId[ii][2]+1<<"\n";
    }
  }else{
    for(size_t ii=0;ii<t.size();ii++){
      out<<fTok<<" "<<t[ii][0]+1<<" "<<
          t[ii][1]+1<<" "<<t[ii][2]+1<<"\n";
    }
  }
  out<<"#end\n";
}

void Mesh::save(const char * filename)
{
  std::ofstream out;
  out.open(filename);
  save(out);
  out.close();
}


void Mesh::load(std::istream &in)
{
  read_obj(in);
}

void Mesh::read_obj(std::istream & f)
{
  std::string line;
  std::string vTok("v");
  std::string fTok("f");
  std::string texTok("vt");
  char bslash='/',space=' ';
  std::string tok;
  while(1) {
    std::getline(f,line);
    if(f.eof()) {
      break;
    }
    if(line == "#end"){
      break;
    }
    if(line.size()<3) {
      continue;
    }
    if(line.at(0)=='#') {
      continue;
    }
    std::stringstream ss(line);
    ss>>tok;
    if(tok==vTok) {
      CompFab::Vec3 vec;
      ss>>vec[0]>>vec[1]>>vec[2];
      v.push_back(vec);
    } else if(tok==fTok) {
      bool hasTexture = false;
      if (line.find(bslash) != std::string::npos) {
        std::replace(line.begin(), line.end(), bslash, space);
        hasTexture = true;
      }
      std::stringstream facess(line);
      facess>>tok;
      std::vector<int> vidx;
      std::vector<int> texIdx;
      int x;
      while(facess>>x){
        vidx.push_back(x);
        if(hasTexture){
          facess>>x;
          texIdx.push_back(x);
        }
      }
      texIdx.resize(vidx.size());
      for(int ii = 0;ii<vidx.size()-2;ii++){
        CompFab::Vec3i trig, textureId;
        trig[0] = vidx[0]-1;
        textureId[0] = texIdx[0]-1;
        for (int jj = 1; jj < 3; jj++) {
          trig[jj] = vidx[ii+jj]-1;
          textureId[jj] = texIdx[ii+jj]-1;
        }
        t.push_back(trig);
        texId.push_back(textureId);
      }
    } else if(tok==texTok) {
        CompFab::Vec2f texcoord;
        ss>>texcoord[0];
        ss>>texcoord[1];
        tex.push_back(texcoord);
    }
  }
  std::cout<<"Num Triangles: "<< t.size()<<"\n";
}

void Mesh::read_ply(std::istream & f)
{
  std::string line;
  std::string vertLine("element vertex");
  std::string faceLine("element face");
  std::string texLine("property float s");
  std::string endHeaderLine("end_header");
  while(true) {
    std::getline(f,line);
    if(std::string::npos!=line.find(vertLine)) {
      break;
    }
  }
  std::string token;
  std::stringstream ss(line);
  ss>>token>>token;
  int nvert;
  ss>>nvert;
  bool hasTex=false;
  while(true) {
    std::getline(f,line);
    if(std::string::npos!=line.find(faceLine)) {
      break;
    }
    if(std::string::npos!=line.find(texLine)) {
      hasTex=true;
    }
  }
  std::stringstream ss1(line);
  ss1>>token>>token;
  int nface;
  ss1>>nface;
  while(true) {
    std::getline(f,line);
    if(std::string::npos!=line.find(endHeaderLine)) {
      break;
    }
  }

  v.resize(nvert);
  t.resize(nface);
  if(hasTex) {
    tex.resize(nvert);
  }
  for (int ii =0; ii<nvert; ii++) {
    for (int jj=0; jj<3; jj++) {
      f>>v[ii][jj];
    }
    if(hasTex) {
      for (int jj=0; jj<2; jj++) {
        f>>tex[ii][jj];
      }
      tex[ii][1]=1-tex[ii][1];;
    }
  }
  for (int ii =0; ii<nface; ii++) {
    int nidx;
    f>>nidx;
    for (int jj=0; jj<3; jj++) {
      f>>t[ii][jj];
    }
  }
}

void Mesh::save_obj(const char * filename)
{
  std::ofstream out(filename);
  if(!out.good()){
    std::cout<<"cannot open output file"<<filename<<"\n";
    return;
  }
  save(out);
  out.close();
}

void Mesh::update()
{}

Mesh::Mesh(const char * filename,bool normalize)
{
  load_mesh(filename,normalize);
}


void Mesh::load_mesh(const char * filename, bool normalize)
{
  std::ifstream f ;
  f.open(filename);
  if(!f.good()) {
    std::cout<<"Error: cannot open mesh "<<filename<<"\n";
    return;
  }
  switch(filename[strlen(filename)-1]) {
  case 'y':
    read_ply(f);
    break;
  case 'j':
    read_obj(f);
    break;
  default:
    break;
  }
  if(normalize){
    rescale();
  }
  compute_norm();

  f.close();
}

void Mesh::rescale()
{
  if(v.size()==0){
    std::cout<<"empty mesh\n";
    return;
  }
  CompFab::Vec3 mn=v[0],mx=v[0];

  //scale and translate to [0 , 1]
  for (unsigned int dim = 0; dim<3; dim++) {
    for( size_t ii=0; ii<v.size(); ii++) {
      mn[dim]= std::min(v[ii][dim],mn[dim]);
      mx[dim] = std::max(v[ii][dim],mx[dim]);
    }
    real_t translate = -mn[dim];
    for(size_t ii=0; ii<v.size(); ii++) {
      v[ii][dim]=(v[ii][dim]+translate);
    }
  }

  real_t scale = 1/(mx[0]-mn[0]);
  for(unsigned int dim=1; dim<3; dim++) {
    scale=std::min(1/(mx[dim]-mn[dim]),scale);
  }

  for(size_t ii=0; ii<v.size(); ii++) {
    for (unsigned int dim = 0; dim<3; dim++) {
      v[ii][dim]=v[ii][dim]*scale;
    }
  }
}

void Mesh::compute_norm()
{
    CompFab::Vec3 ZERO;
    
    n.resize(v.size(), ZERO);
  for(unsigned int ii=0; ii<t.size(); ii++) {
    CompFab::Vec3 a = v[t[ii][1]] - v[t[ii][0]];
    CompFab::Vec3 b = v[t[ii][2]] - v[t[ii][0]];
    b=a%b;
    b.normalize();
    for(int jj=0; jj<3; jj++) {
      n[t[ii][jj]]+=b;
      if(t[ii][jj]>=(int)n.size() || t[ii][jj]<0){
        std::cout<<ii<<" "<<jj<<" "<<t[ii][jj]<<" normal computation error\n";
      }
    }
  }
  for(unsigned int ii=0; ii<v.size(); ii++) {
    n[ii].normalize();
  }
}

void BBox(const Mesh & m,
    CompFab::Vec3 & mn, CompFab::Vec3 & mx)
{
  BBox(m.v, mn, mx);
}

bool is_nbr(const CompFab::Vec3i & a, const CompFab::Vec3i&b, int vert)
{
  for (int ii=0; ii<3; ii++) {

    int va=a[ii];
    if(va<=vert) {
      continue;
    }

    for (unsigned int jj=0; jj<3; jj++) {
      int vb=b[jj];
      if(vb<=vert) {
        continue;
      }
      if(va==vb) {
        return true;
      }
    }
  }
  return false;
}


void adjlist(const Mesh & m, std::vector<std::vector<int> > & adjMat)
{
  if(adjMat.size()==m.t.size()) {
    return;
  }
  std::vector<std::vector<int> >trigList;
  trigList.resize(m.v.size());
  for (unsigned int ii=0; ii<m.t.size(); ii++) {
    for (unsigned int jj=0; jj<3; jj++) {
      int vidx=m.t[ii][jj];
      trigList[vidx].push_back(ii);
    }
  }
  adjMat.resize(m.t.size());
  for (unsigned int ii=0; ii<m.v.size(); ii++) {
    int n_nbr=trigList[ii].size();
    for (int jj=0; jj<n_nbr; jj++) {
      int tj=trigList[ii][jj];
      for (int kk=(jj+1); kk<n_nbr; kk++) {
        int tk=trigList[ii][kk];
        if(is_nbr(m.t[tj],m.t[tk],ii)) {
          adjMat[tj].push_back(tk);
          adjMat[tk].push_back(tj);
        }

      }
    }
  }
}


void BBox(const std::vector<CompFab::Vec3 >& v,
    CompFab::Vec3 & mn, CompFab::Vec3 & mx)
{
  mn = v[0];
  mx = v[0];
  for(unsigned int ii = 1 ;ii<v.size();ii++){
    for(int dim = 0 ; dim<3;dim++){
      if(v[ii][dim]<mn[dim]){
        mn[dim] = v[ii][dim];
      }
      if(v[ii][dim]>mx[dim]){
        mx[dim] = v[ii][dim];
      }
    }
  }
}

CompFab::Vec3 vectorMax(const CompFab::Vec3 &a, const CompFab::Vec3 &b) {
  return CompFab::Vec3(std::max(a.m_x,b.m_x), std::max(a.m_y,b.m_y), std::max(a.m_z,b.m_z));
}

CompFab::Vec3 vectorMin(const CompFab::Vec3 &a, const CompFab::Vec3 &b) {
  return CompFab::Vec3(std::min(a.m_x,b.m_x), std::min(a.m_y,b.m_y), std::min(a.m_z,b.m_z));
}

BoundingBox::BoundingBox(const CompFab::Vec3 &minP, const CompFab::Vec3 &maxP) {
  mMinP = minP - EPSILON_VEC;
  mMaxP = maxP + EPSILON_VEC;

  mExtents = (maxP - minP) * 0.5;
}

BoundingBox::BoundingBox(std::vector<CompFab::Triangle> &faces) {
    if(faces.size() == 0) return;

    mMaxP = CompFab::Vec3(-DBL_INF, -DBL_INF, -DBL_INF);
    mMinP = CompFab::Vec3(DBL_INF, DBL_INF, DBL_INF);

    for(int i = 0; i < (int)faces.size(); i++) {
        for(int v = 0; v < 3; v++) {
            mMaxP = vectorMax(mMaxP, faces[i][v]);
            mMinP = vectorMin(mMinP, faces[i][v]);
        }
    }

    mExtents = (mMaxP - mMinP) * 0.5;
}

BoundingBox::BoundingBox() {
}


int BoundingBox::longestAxis() {
    int axis;
    double longest = -DBL_INF;

    for(int i = 0; i < 3; i++) {
      double len = mMaxP[i] - mMinP[i];
      if(len > longest) {
        longest = len;
        axis = i;
      }
    }
    return axis;
}

int BoundingBox::hit(const CompFab::Ray &ray) {
  //Alg from http://tavianator.com/2011/05/fast-branchless-raybounding-box-intersections/
  CompFab::Vec3 rayOrg = ray.m_origin;
  CompFab::Vec3 rayDir = ray.m_direction;

  double tx1 = (mMinP.m_x - rayOrg.m_x)/rayDir.m_x;
  double tx2 = (mMaxP.m_x - rayOrg.m_x)/rayDir.m_x;

  double tmin = std::min(tx1, tx2);
  double tmax = std::max(tx1, tx2);

  double ty1 = (mMinP.m_y - rayOrg.m_y)/rayDir.m_y;
  double ty2 = (mMaxP.m_y - rayOrg.m_y)/rayDir.m_y;

  tmin = std::max(tmin, std::min(ty1, ty2));
  tmax = std::min(tmax, std::max(ty1, ty2));

  double tz1 = (mMinP.m_z - rayOrg.m_z)/rayDir.m_z;
  double tz2 = (mMaxP.m_z - rayOrg.m_z)/rayDir.m_z;

  tmin = std::max(tmin, std::min(tz1, tz2));
  tmax = std::min(tmax, std::max(tz1, tz2));

  if(tmax >= 0.0 && tmin >= 0.0 && tmax >= tmin) {
    return 2;
  } else if (tmax >= 0.0 and tmin <= 0.0) {
    return 1;
  }
  else {
    return 0;
  }
}

KdNode::KdNode() {
}

KdNode::~KdNode() {
    delete mLeft;
    delete mRight;
}


KdNode* KdNode::build(std::vector<CompFab::Triangle> &faces, int depth) const {
    KdNode *node = new KdNode();
    node->mFaces = faces;
    node->mLeft = NULL;
    node->mRight = NULL;
    node->mBbox = BoundingBox(faces);

    //cout << depth << endl;
    if(faces.size() == 0 || faces.size() == 1) {
        //cout << "hello 0" << endl;
        return node;
    } 

    CompFab::Vec3 midpt;
    for(int i = 0; i < faces.size(); i++) {
      for (int j = 0; j < 3; ++j) {
        midpt = midpt + faces[i][j];
      }
    }
    midpt = midpt * (1.0 / ((double)faces.size() * 3.0));

    std::vector<CompFab::Triangle> leftFaces;
    std::vector<CompFab::Triangle> rightFaces;
    int axis = node->mBbox.longestAxis();

    for(int i = 0; i < (int)faces.size(); i++) {
        CompFab::Vec3 triMidpt;
        for (int j = 0; j < 3; ++j) {
          triMidpt = triMidpt + faces[i][j];
        }
        triMidpt = triMidpt * (1.0 / (3.0));

        if(midpt[axis] >= triMidpt[axis]) { //Cheat and always take first vert of tri?
            rightFaces.push_back(faces[i]);
        } else {
            leftFaces.push_back(faces[i]);
        }
    }


    if(leftFaces.size() == 0 || rightFaces.size() == 0) {
        return node;
    } else {
        node->mLeft = build(leftFaces, depth + 1);
        node->mRight = build(rightFaces, depth + 1);
    }

    return node;
}


int KdNode::rayIntersection(const CompFab::Ray &ray, int depth) {

    if(mBbox.hit(ray)) {
        // if(depth >= 7) {
        //   return mBbox.hit(ray);
        // }

        if(mLeft) {
            int hitLeft = mLeft->rayIntersection(ray, depth + 1);
            int hitRight = mRight->rayIntersection(ray, depth + 1);
            return hitLeft + hitRight;
        }
        else {
            int hitTri = 0;
            for(int i = 0; i < (int)mFaces.size(); i++) {
                hitTri += rayTriangleIntersection(ray, mFaces[i]);
            }

            return hitTri;
        }
    } 
    else {
        return 0;
    }

    return 0;
}

int rayTriangleIntersection(const CompFab::Ray &ray, CompFab::Triangle &triangle)
{
    /********* ASSIGNMENT *********/
    /* Ray-Triangle intersection test: Return 1 if ray intersects triangle, 
     * 0 otherwise */

    // Implementation based on Möller–Trumbore intersection algorithm
    // http://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm

    const CompFab::Vec3 &V1 = triangle.m_v1;
    const CompFab::Vec3 &V2 = triangle.m_v2;
    const CompFab::Vec3 &V3 = triangle.m_v3;

    CompFab::Vec3 e1, e2;  //Edge1, Edge2
    CompFab::Vec3 P, Q, T;
    double det, inv_det, u, v;
    double t; // Ray intersection distance

    //Find vectors for two edges sharing V1
    e1 = V2 - V1;
    e2 = V3 - V1;

    //Begin calculating determinant - also used to calculate u parameter
    P = ray.m_direction % e2;
    //if determinant is near zero, ray lies in plane of triangle
    det = e1 * P;
    //NOT CULLING
    if(det > -TRI_EPSILON && det < TRI_EPSILON) return 0;
    inv_det = 1.0 / det;

    //calculate distance from V1 to ray origin
    T = ray.m_origin - V1;

    //Calculate u parameter and test bound
    u = (T * P) * inv_det;

    //The intersection lies outside of the triangle
    if(u < 0.0 || u > 1.0) return 0;

    //Prepare to test v parameter
    Q = T % e1;

    //Calculate V parameter and test bound
    v = (ray.m_direction * Q) * inv_det;
    //The intersection lies outside of the triangle
    if(v < 0.0 || u + v  > 1.0) return 0;

    t = (e2 * Q) * inv_det;

    if(t > TRI_EPSILON) { // ray intersection
        return 1;
    }

    return 0;
}
