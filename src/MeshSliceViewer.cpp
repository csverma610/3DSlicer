#include "MeshSliceViewer.hpp"
#include <iostream>
#include <fstream>
#include <GL/glut.h>

using namespace std;


void MeshSliceViewer::init() {
    glEnable(GL_DEPTH_TEST);
    glEnable ( GL_STENCIL_TEST );


    sliceDisplayID = glGenLists(1);
//   primitives.push_back(new OpenCSG::DisplayListPrimitive(sliceDisplayID, OpenCSG::Intersection, 1));

    OpenCSG::setOption(OpenCSG::AlgorithmSetting, OpenCSG::Automatic);
    OpenCSG::setOption(OpenCSG::OffscreenSetting, OpenCSG::AutomaticOffscreenType);

    GLuint id1 = glGenLists(1);
    glNewList(id1, GL_COMPILE);
    glPushMatrix();
    glTranslatef(-0.25f, 0.0f, 0.0f);
    glutSolidSphere(1.0, 20, 20);
    glPopMatrix();
    glEndList();

    GLuint id2 = glGenLists(1);
    glNewList(id2, GL_COMPILE);
    glPushMatrix();
    glTranslatef(0.25f, 0.0f, 0.0f);
    glutSolidSphere(1.0, 20, 20);
    glPopMatrix();
    glEndList();

    GLuint id3 = glGenLists(1);
    glNewList(id3, GL_COMPILE);
    glPushMatrix();
    glTranslatef(0.0f, 0.0f, 0.5f);
    glScalef(0.5f, 0.5f, 2.0f);
    glutSolidSphere(1.0, 20, 20);
    glPopMatrix();
    glEndList();

    primitives.push_back(new OpenCSG::DisplayListPrimitive(id1, OpenCSG::Intersection, 1));
    primitives.push_back(new OpenCSG::DisplayListPrimitive(id2, OpenCSG::Intersection, 1));
    primitives.push_back(new OpenCSG::DisplayListPrimitive(id3, OpenCSG::Subtraction, 1));

}
//////////////////////////////////////////////////////////////

void MeshSliceViewer::readMesh( const string &filename)
{
    mesh.nodes.clear();
    mesh.edges.clear();
    mesh.faces.clear();

    ifstream ifile(filename.c_str(), ios::in);
    if( ifile.fail() ) {
        cout << "Error: Can't open input file " << endl;
        return;
    }

    string str;
    ifile >> str;
    if( str != "OFF") {
        cout << "Error: input file not in the off format" << endl;
        return;
    }

    int numPoints, numFaces, numEdges;
    ifile >> numPoints >> numFaces >> numEdges;

    if( numPoints < 1) {
        cout << "Warning: Input file has no points " << endl;
        return;
    }

    if( numPoints ) {
        float x, y, z;
        mesh.nodes.resize(numPoints);
        for( size_t i = 0; i < numPoints; i++) {
            ifile >> x >> y >> z;
            mesh.nodes[i] = {x,y,z};
        }
    }

    int n0, n1;
    if( numFaces) {
        set<pair<size_t,size_t>> eSet;

        mesh.faces.resize(numFaces);

        int nn;
        for( size_t i = 0; i < numFaces; i++) {
            ifile >> nn; assert(nn == 3);
            for( int j = 0; j < nn; j++)
                ifile >> mesh.faces[i][j];
            for( int j = 0; j < nn; j++) {
                int v0 = mesh.faces[i][(j+1)%nn];
                int v1 = mesh.faces[i][(j+2)%nn];
                int vmin = min(v0,v1);
                int vmax = max(v0,v1);
                eSet.insert( make_pair(vmin,vmax));
            }
        }
        mesh.edges.resize(eSet.size());
        int index = 0;
        for(auto edge: eSet) {
            n0 = edge.first;
            n1 = edge.second;
            assert( n0 != n1);
            mesh.edges[index++]  = {n0,n1};
        }
    }

    if( numEdges ) {
        mesh.edges.resize(numEdges);
        for(size_t i = 0; i < numEdges; i++) {
            ifile >> n0 >> n1;
            mesh.edges[i]  = {n0,n1};
        }
    }
    mesh_initialized = 0;

}
//////////////////////////////////////////////////////////////
void MeshSliceViewer:: initMesh()
{
    meshDisplayID = glGenLists(1);

//    primitives.push_back(new OpenCSG::DisplayListPrimitive(meshDisplayID, OpenCSG::Intersection, 1));

    size_t numNodes = mesh.nodes.size();

    assert( numNodes );

    float xmin = std::numeric_limits<float>::max();
    float ymin = std::numeric_limits<float>::max();
    float zmin = std::numeric_limits<float>::max();

    float xmax = -0.99*xmin;
    float ymax = -0.99*ymin;
    float zmax = -0.99*zmin;

    for( size_t i = 0; i < numNodes; i++) {
        float x = mesh.nodes[i][0];
        float y = mesh.nodes[i][1];
        float z = mesh.nodes[i][2];
        xmin = min( xmin, x);
        xmax = max( xmax, x);

        ymin = min( ymin, y);
        ymax = max( ymax, y);

        zmin = min( zmin, z);
        zmax = max( zmax, z);
    }

    float xc = 0.5*(xmax + xmin);
    float yc = 0.5*(ymax + ymin);
    float zc = 0.5*(zmax + zmin);

    for( size_t i = 0; i < numNodes; i++) {
        mesh.nodes[i][0] -= xc;
        mesh.nodes[i][1] -= yc;
        mesh.nodes[i][2] -= zc;
    }

    float xlen = fabs(xmax - xmin);
    float ylen = fabs(ymax - ymin);
    float zlen = fabs(zmax - zmin);
    float maxlen = max(zlen, max(xlen,ylen));

    for( size_t i = 0; i < numNodes; i++) {
        mesh.nodes[i][0] /= maxlen;
        mesh.nodes[i][1] /= maxlen;
        mesh.nodes[i][2] /= maxlen;
    }

    xmin = std::numeric_limits<float>::max();
    ymin = std::numeric_limits<float>::max();
    zmin = std::numeric_limits<float>::max();
    xmax = -0.99*xmin;
    ymax = -0.99*ymin;
    zmax = -0.99*zmin;
    for( size_t i = 0; i < numNodes; i++) {
        float x = mesh.nodes[i][0];
        float y = mesh.nodes[i][1];
        float z = mesh.nodes[i][2];
        xmin = min( xmin, x);
        xmax = max( xmax, x);

        ymin = min( ymin, y);
        ymax = max( ymax, y);

        zmin = min( zmin, z);
        zmax = max( zmax, z);
    }

    minCorner = {xmin, ymin, zmin};
    maxCorner = {xmax, ymax, zmax};

    sliceID   = 0;
}
//////////////////////////////////////////////////////////////

void MeshSliceViewer:: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        return;
    }

    /*
        if( e->key() == Qt::Key_1) {
            displayPlane[0] = !displayPlane[0];
        }
        if( e->key() == Qt::Key_2) {
            displayPlane[1] = !displayPlane[1];
        }

        if( e->key() == Qt::Key_3) {
            displayPlane[2] = !displayPlane[2];
        }

        if( e->key() == Qt::Key_B) {
            displayBox = !displayBox;
        }

        if( e->key() == Qt::Key_A) {
            displayAxis = !displayAxis;
        }
    */
    if( e->key() == Qt::Key_Plus)  {
        sliceID++;
        slice_initialized = 0;
    }

    if( e->key() == Qt::Key_Minus) {
        sliceID--;
        slice_initialized = 0;
    }

    /*
        if( e->key() == Qt::Key_R) {
            currSliceID[0] = minCorner[0];
            currSliceID[1] = minCorner[1];
            currSliceID[2] = minCorner[2];
        }

        if( e->key() == Qt::Key_C) {
          xSlice.clear();
          ySlice.clear();
          zSlice.clear();
        }

        if( e->key() == Qt::Key_X) {
            slicer.getSlice(currSliceID[0], 0, xSlice);
            cout <<"XSliceIndex " << currSliceID[0] << " #Voxels: " << xSlice.size() << endl;
            currSliceID[0] = currSliceID[0] + increment;
            if( currSliceID[0] >= maxCorner[0] )
                currSliceID[0] = minCorner[0];
        }

        if( e->key() == Qt::Key_Y) {
            slicer.getSlice(currSliceID[1], 1, ySlice);
            cout <<"YSliceIndex " << currSliceID[1] << " #Voxels: " << ySlice.size() << endl;
            currSliceID[1] = currSliceID[1] + increment;
            if( currSliceID[1] >= maxCorner[1] )
                currSliceID[1] = minCorner[1];

        }

        if( e->key() == Qt::Key_Z) {
            slicer.getSlice(currSliceID[2], 2, zSlice);
            currSliceID[2] = currSliceID[2] + increment;
            if( currSliceID[2] >= maxCorner[2] )
                currSliceID[2] = minCorner[2];
        }

        if( e->key() == Qt::Key_G)
            displayGrid = !displayGrid;

        if( e->key() == Qt::Key_P) {
            displayPointCloud = !displayPointCloud;
            if( displayPointCloud ) getActiveVoxelCloud();
        }

        if( e->key() == Qt::Key_H) {
            this->showEntireScene();
        }
    */
    QGLViewer::keyPressEvent(e);
    this->update();
}
//////////////////////////////////////////////////////////////

void MeshSliceViewer::drawCube(const Array3F &origin, const Array3F &length)
{
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL);
    glColor3f( 1.0, 0.0, 0.0);
    float x = origin[0];
    float y = origin[1];
    float z = origin[2];

    float dx = length[0];
    float dy = length[1];
    float dz = length[2];

    glNormal3f( -1.0, 0.0, 0.0);
    glBegin(GL_QUADS);
    glVertex3f(x,y,z);
    glVertex3f(x,y,z+dz);
    glVertex3f(x,y+dy, z+dz);
    glVertex3f(x,y+dy,z);
    glEnd();

    glNormal3f( 1.0, 0.0, 0.0);
    glBegin(GL_QUADS);
    glVertex3f(x+dx,y,z);
    glVertex3f(x+dx,y+dy,z);
    glVertex3f(x+dx,y+dy, z+dz);
    glVertex3f(x+dx,y,z+dz);
    glEnd();

    glNormal3f(0.0, -1.0, 0.0);
    glBegin(GL_QUADS);
    glVertex3f(x,y,z);
    glVertex3f(x+dx,y,z);
    glVertex3f(x+dx,y, z+dz);
    glVertex3f(x,y,z+dz);
    glEnd();

    glNormal3f(0.0, 1.0, 0.0);
    glBegin(GL_QUADS);
    glVertex3f(x,y+dy,z);
    glVertex3f(x,y+dy,z+dz);
    glVertex3f(x+dx,y+dy, z+dz);
    glVertex3f(x+dx,y+dy,z);
    glEnd();

    glNormal3f(0.0, 0.0, -1.0);
    glBegin(GL_QUADS);
    glVertex3f(x,y,z);
    glVertex3f(x,y+dy,z);
    glVertex3f(x+dx,y+dy, z);
    glVertex3f(x+dx,y, z);
    glEnd();

    glNormal3f(0.0, 0.0, +1.0);
    glBegin(GL_QUADS);
    glVertex3f(x,y,z+dz);
    glVertex3f(x+dx,y,z+dz);
    glVertex3f(x+dx,y+dy, z+dz);
    glVertex3f(x,y+dy, z+dz);
    glEnd();
}

//////////////////////////////////////////////////////////////

void MeshSliceViewer::mouseMoveEvent(QMouseEvent *e)
{
    /*
      if( this->selectedName() >= 0) {
      QPoint qPoint = e->pos();
      qglviewer::Vec pixelPoint;
      pixelPoint[0] = qPoint.x();
      pixelPoint[1] = qPoint.y();
      qglviewer::Vec  point = camera()->unprojectedCoordinatesOf(pixelPoint);
      }
    */
    QGLViewer::mouseMoveEvent(e);
}
//////////////////////////////////////////////////////////////

void MeshSliceViewer::selectPlane()
{
//   int pid = this->selectedName();
}
//////////////////////////////////////////////////////////////
void MeshSliceViewer::draw()
{
    if( !mesh_initialized ) {
        initMesh();
        glNewList(meshDisplayID, GL_COMPILE);
        glPushMatrix();
        glColor3f( 1.0, 1.0, 1.0);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        for( size_t i = 0; i < mesh.faces.size(); i++) {
            int v0 = mesh.faces[i][0];
            int v1 = mesh.faces[i][1];
            int v2 = mesh.faces[i][2];
            glBegin(GL_TRIANGLES);
            glVertex3f( mesh.nodes[v0][0], mesh.nodes[v0][1], mesh.nodes[v0][2] );
            glVertex3f( mesh.nodes[v1][0], mesh.nodes[v1][1], mesh.nodes[v1][2] );
            glVertex3f( mesh.nodes[v2][0], mesh.nodes[v2][1], mesh.nodes[v2][2] );
            glEnd();
        }

        /*
                glColor3f( 1.0, 0.0, 0.0);
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                for( size_t i = 0; i < mesh.faces.size(); i++) {
                    int v0 = mesh.faces[i][0];
                    int v1 = mesh.faces[i][1];
                    int v2 = mesh.faces[i][2];
                    glBegin(GL_TRIANGLES);
                    glVertex3f( mesh.nodes[v0][0], mesh.nodes[v0][1], mesh.nodes[v0][2] );
                    glVertex3f( mesh.nodes[v1][0], mesh.nodes[v1][1], mesh.nodes[v1][2] );
                    glVertex3f( mesh.nodes[v2][0], mesh.nodes[v2][1], mesh.nodes[v2][2] );
                    glEnd();
                }
        */
        glPopMatrix();
        glEndList();
        mesh_initialized = 1;
    }

    if( !slice_initialized ) {
        if( mesh_initialized) {
            glNewList(sliceDisplayID, GL_COMPILE);
            glPushMatrix();
            float dz   = (maxCorner[2] - minCorner[2])/100.0;
            float xlen = maxCorner[0] - minCorner[0];
            float ylen = maxCorner[1] - minCorner[1];
            Array3F origin = {minCorner[0], minCorner[1], minCorner[2] + sliceID*dz};
            Array3F length = {xlen, ylen,dz};
            drawCube(origin, length);
            glPopMatrix();
            glEndList();
            slice_initialized = 1;
            this->update();
        }
    }

    /*
        assert( primitives.size() == 2);
        primitives[0]->render();
        primitives[1]->render();
    */

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0.0, 2.0, 5.0,  /* eye is at (0,2,5) */
              0.0, 0.0, 0.0,  /* center is at (0,0,0) */
              0.0, 1.0, 0.0); /* up is in positive Y direction */
    glRotatef(0.0, 0.0f, 1.0f, 0.0f);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    OpenCSG::render(primitives);

    glDepthFunc(GL_EQUAL);
    for (auto i = primitives.begin(); i != primitives.end(); ++i) {
        (*i)->render();
    }
    glDepthFunc(GL_LESS);

    /*
        if(displayAxis) {
            glEnable(GL_LIGHTING);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glPushMatrix();
            float maxlen = max(boxLength[0], boxLength[1] );
            maxlen = max(maxlen, boxLength[1]);
            drawAxis(maxlen);
            glPopMatrix();
        }

        if(displayPointCloud) {
            glDisable(GL_LIGHTING);
            glEnable(GL_FOG);
            glColor3f(1.0, 0.0, 0.0);
            glPointSize(1);
            glPushMatrix();
            glBegin(GL_POINTS);
            for( auto xyz: pointCloud)
                glVertex3f(xyz[0], xyz[1], xyz[2] );
            glEnd();
            glPopMatrix();
            glDisable(GL_FOG);
        }


        if( displayBox ) {
            glDisable(GL_LIGHTING);
            glColor3f(1.0, 1.0, 1.0);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glPushMatrix();
            drawCube( origin[0], boxLength[0],
                      origin[1], boxLength[1],
                      origin[2], boxLength[2] );
            glPopMatrix();
        }

        if( displayGrid ) {
            glEnable(GL_POLYGON_OFFSET_FILL);
            glPolygonOffset(1.0, 2.0);
            glDisable(GL_LIGHTING);
            glColor3f(0.1, 0.1, 0.1);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glPushMatrix();
            for( const auto &c : xSlice) drawCube( c[0], c[1], c[2] );
            for( const auto &c : ySlice) drawCube( c[0], c[1], c[2] );
            for( const auto &c : zSlice) drawCube( c[0], c[1], c[2] );
            glPopMatrix();
        }


        glColor3f(1.0, 1.0, 1.0);
        glEnable(GL_LIGHTING);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glPushMatrix();
        for( const auto &c : xSlice) drawCube( c[0], c[1], c[2] );
        for( const auto &c : ySlice) drawCube( c[0], c[1], c[2] );
        for( const auto &c : zSlice) drawCube( c[0], c[1], c[2] );
        glPopMatrix();

        glDisable(GL_POLYGON_OFFSET_FILL);

        glDisable(GL_DEPTH_TEST);
        glDisable(GL_LIGHTING);
        if( displayPlane[0] ) drawXPlane(currSliceID[0]-1);
        if( displayPlane[1] ) drawYPlane(currSliceID[1]-1);
        if( displayPlane[2] ) drawZPlane(currSliceID[2]-1);
        glEnable(GL_DEPTH_TEST);
    */
}

int main(int argc, char **argv)
{
    if( argc != 2)  {
        cout << "Usage: " << argv[0] << " meshfile " << endl;
        return 1;
    }
    QApplication application(argc, argv);
    glutInit(&argc, argv);


    MeshSliceViewer viewer;

    viewer.readMesh( argv[1] );
    viewer.setWindowTitle("Mesh Slicer");
    viewer.show();
    return application.exec();
}
///////////////////////////

void MeshSliceViewer::drawXPlane(int i)
{
    /*
        glColor4f(1.0, 0.0, 0.0, 0.1);

        glEnable(GL_BLEND);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        Point3F xyz;
        glBegin(GL_QUADS);
        xyz = getWorldCoord(i, minCorner[1], minCorner[2]);
        xyz[0] += 0.5*voxelSize[0];
        glVertex3fv(&xyz[0]);

        xyz = getWorldCoord(i, maxCorner[1]+1, minCorner[2]);
        xyz[0] += 0.5*voxelSize[0];
        glVertex3fv(&xyz[0]);

        xyz = getWorldCoord(i, maxCorner[1]+1, maxCorner[2]+1);
        xyz[0] += 0.5*voxelSize[0];
        glVertex3fv(&xyz[0]);

        xyz = getWorldCoord(i, minCorner[1], maxCorner[2]+1);
        xyz[0] += 0.5*voxelSize[0];
        glVertex3fv(&xyz[0]);

        glEnd();
        glDisable(GL_BLEND);

        glColor3f(1.0, 0.0, 0.0);
        glLineWidth(2);

        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

        glBegin(GL_QUADS);
        xyz = getWorldCoord(i, minCorner[1], minCorner[2]);
        xyz[0] += 0.5*voxelSize[0];
        glVertex3fv(&xyz[0]);

        xyz = getWorldCoord(i, maxCorner[1]+1, minCorner[2]);
        xyz[0] += 0.5*voxelSize[0];
        glVertex3fv(&xyz[0]);

        xyz = getWorldCoord(i, maxCorner[1]+1, maxCorner[2]+1);
        xyz[0] += 0.5*voxelSize[0];
        glVertex3fv(&xyz[0]);

        xyz = getWorldCoord(i, minCorner[1], maxCorner[2]+1);
        xyz[0] += 0.5*voxelSize[0];
        glVertex3fv(&xyz[0]);

        glEnd();
        glLineWidth(1);
    */

}

///////////////////////////////////////////////////////////////////////////////

void MeshSliceViewer::drawYPlane(int j)
{
    /*
        Point3F xyz;
        glColor4f(0.0, 1.0, 0.0, 0.1);

        glEnable(GL_BLEND);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        glBegin(GL_QUADS);
        xyz = getWorldCoord(minCorner[0], j, minCorner[2]);
        xyz[1] += 0.5*voxelSize[1];
        glVertex3fv(&xyz[0]);

        xyz = getWorldCoord(maxCorner[0]+1, j, minCorner[2]);
        xyz[1] += 0.5*voxelSize[1];
        glVertex3fv(&xyz[0]);

        xyz = getWorldCoord(maxCorner[0]+1, j, maxCorner[2]+1);
        xyz[1] += 0.5*voxelSize[1];
        glVertex3fv(&xyz[0]);

        xyz = getWorldCoord(minCorner[0], j, maxCorner[2]+1);
        xyz[1] += 0.5*voxelSize[1];
        glVertex3fv(&xyz[0]);

        glEnd();
        glDisable(GL_BLEND);

        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

        glLineWidth(2);
        glColor3f(0.0, 1.0, 0.0);

        glBegin(GL_QUADS);
        xyz = getWorldCoord(minCorner[0], j, minCorner[2]);
        xyz[1] += 0.5*voxelSize[1];
        glVertex3fv(&xyz[0]);

        xyz = getWorldCoord(maxCorner[0]+1, j, minCorner[2]);
        xyz[1] += 0.5*voxelSize[1];
        glVertex3fv(&xyz[0]);

        xyz = getWorldCoord(maxCorner[0]+1, j, maxCorner[2]+1);
        xyz[1] += 0.5*voxelSize[1];
        glVertex3fv(&xyz[0]);

        xyz = getWorldCoord(minCorner[0], j, maxCorner[2]+1);
        xyz[1] += 0.5*voxelSize[1];
        glVertex3fv(&xyz[0]);

        glEnd();
        glLineWidth(1);
    */
}

///////////////////////////////////////////////////////////////////////////////

void MeshSliceViewer::drawZPlane(int k)
{
    /*
        Point3F xyz;
        glColor4f(0.0, 0.0, 1.0, 0.1);

        glEnable(GL_BLEND);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        glBegin(GL_QUADS);
        xyz = getWorldCoord(minCorner[0], minCorner[1], k);
        xyz[2] += 0.5*voxelSize[2];
        glVertex3fv(&xyz[0]);

        xyz = getWorldCoord(maxCorner[0]+1, minCorner[1], k);
        xyz[2] += 0.5*voxelSize[2];
        glVertex3fv(&xyz[0]);

        xyz = getWorldCoord(maxCorner[0]+1, maxCorner[1]+1, k);
        xyz[2] += 0.5*voxelSize[2];
        glVertex3fv(&xyz[0]);

        xyz = getWorldCoord(minCorner[0], maxCorner[1]+1, k);
        xyz[2] += 0.5*voxelSize[2];
        glVertex3fv(&xyz[0]);

        glEnd();
        glDisable(GL_BLEND);

        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

        glLineWidth(2);
        glColor3f(0.0, 0.0, 1.0);

        glBegin(GL_QUADS);
        xyz = getWorldCoord(minCorner[0], minCorner[1], k);
        xyz[2] += 0.5*voxelSize[2];
        glVertex3fv(&xyz[0]);

        xyz = getWorldCoord(maxCorner[0]+1, minCorner[1], k);
        xyz[2] += 0.5*voxelSize[2];
        glVertex3fv(&xyz[0]);

        xyz = getWorldCoord(maxCorner[0]+1, maxCorner[1]+1, k);
        xyz[2] += 0.5*voxelSize[2];
        glVertex3fv(&xyz[0]);

        xyz = getWorldCoord(minCorner[0], maxCorner[1]+1, k);
        xyz[2] += 0.5*voxelSize[2];
        glVertex3fv(&xyz[0]);

        glEnd();
        glLineWidth(1);
    */
}

