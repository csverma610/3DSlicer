///////////////////////////////////////////////////////////////////////////////
// Written by : Chaman Singh Verma
// Date  :  11th May, 2018
// Code Adapted from OpenCSG example. Offscreen part was taken from Internet.
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <array>
#include <set>
#include <cassert>
#include <limits>
#include <math.h>
#include <iomanip>

#include <GL/glew.h>
#include <GL/glut.h>
#include "displaylistPrimitive.h"
#include <opencsg.h>

using namespace std;

enum {
    GF_STANDARD, GF_DC, GF_OQ, SCS_STANDARD, SCS_DC, SCS_OQ, ALGO_AUTOMATIC,
    OFFSCREEN_AUTOMATIC, OFFSCREEN_FBO, OFFSCREEN_PBUFFER
};

std::vector<OpenCSG::Primitive*> primitives;

bool    spin = 0;
float   rot  = 0.0f;
int     numSlices = 100;
int     sliceID   = 0;
int     sliceDir  = 1;
float   slicePos  = 0.0;
float   sliceStep = 0.0;
float   eyePos    = 2.0;
bool    animate   = 0;

bool    displayPlane = 1;
bool    displayMesh  = 1;
bool    displayIntersection  = 1;
bool    displayBox   = 1;
bool    displayWireframe = 1;
bool    displayAxis = 1;

bool    cameraType   = 0;
float   zoom = 0.0;
bool    full_screen = 0;

int     winWidth  = 1000;
int     winHeight = 1000;
float   frameRatio = 1.0;

std::string filename;

typedef std::array<float,3> Array3F;
typedef std::array<int,2>   Array2I;
typedef std::array<int,3>   Array3I;

struct Mesh
{
    std::vector<Array3F> nodes;
    std::vector<Array2I> edges;
    std::vector<Array3I> faces;
};

Mesh mesh;
Array3F minCorner, maxCorner;
Array3F center;

////////////////////////////////////////////////////////////////////////////////

void clearPrimitives() {
    for (auto i = primitives.begin(); i != primitives.end(); ++i) {
        OpenCSG::DisplayListPrimitive* p =
            static_cast<OpenCSG::DisplayListPrimitive*>(*i);
        glDeleteLists(1, p->getDisplayListId());
        delete p;
    }

    primitives.clear();
}
////////////////////////////////////////////////////////////////////////////////

void setScreen(int w, int h) {
    glutReshapeWindow(winWidth, winHeight);
}

////////////////////////////////////////////////////////////////////////////////

void setCamera()
{
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    if( cameraType == 0)
        gluPerspective(40.0, 1.0, 0.001, 10000.0);
    else {
        float aspect = (maxCorner[0]-minCorner[0] + 2*zoom)/(maxCorner[1]-minCorner[1] + 2*zoom);
        float xmin = minCorner[0] - zoom;
        float xmax = maxCorner[0] + zoom;
        float ymin = minCorner[1] - zoom;
        float ymax = maxCorner[1] + zoom;
        if( aspect > 1.0)
            glOrtho( xmin, xmax, aspect*ymin,  aspect*ymax, 0.00001, 100000.0);
        else
            glOrtho(xmin/aspect, xmax/aspect, ymin, ymax, 0.00001, 100000.0);
    }
}

////////////////////////////////////////////////////////////////////////////////

void readMesh()
{
    cout << "Reading input mesh ..." << endl;
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
}

////////////////////////////////////////////////////////////////////////////////

void initMesh()
{
    cout << "Initializing mesh ..." << endl;
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

    xlen = fabs(xmax - xmin);
    ylen = fabs(ymax - ymin);
    zlen = fabs(zmax - zmin);

    minCorner = {xmin, ymin, zmin};
    maxCorner = {xmax, ymax, zmax};
    slicePos  =  minCorner[2];
    sliceStep = (maxCorner[2] - minCorner[2])/(double)numSlices;
    center[0] = 0.5*(maxCorner[0] + minCorner[0] );
    center[1] = 0.5*(maxCorner[1] + minCorner[1] );
    center[2] = 0.5*(maxCorner[2] + minCorner[2] );
}

////////////////////////////////////////////////////////////////////////////////

void setMesh()
{
    GLuint id = glGenLists(1);

    glNewList(id, GL_COMPILE);
    glPushMatrix();

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
    glPopMatrix();
    glEndList();
    primitives[0] = new OpenCSG::DisplayListPrimitive(id, OpenCSG::Intersection, 1);
}

////////////////////////////////////////////////////////////////////////////////

void xRotation(double angle)
{
    double cost = cos(M_PI*angle/180.0);
    double sint = sin(M_PI*angle/180.0);

    size_t numNodes = mesh.nodes.size();
    for( size_t i = 0; i < numNodes; i++) {
        float y = mesh.nodes[i][1];
        float z = mesh.nodes[i][2];
        mesh.nodes[i][1] = cost*y - sint*z;
        mesh.nodes[i][2] = sint*y + cost*z;
    }
}

////////////////////////////////////////////////////////////////////////////////

void yRotation(double angle)
{
    double cost = cos(M_PI*angle/180.0);
    double sint = sin(M_PI*angle/180.0);

    size_t numNodes = mesh.nodes.size();
    for( size_t i = 0; i < numNodes; i++) {
        float x = mesh.nodes[i][0];
        float z = mesh.nodes[i][2];
        mesh.nodes[i][0] =  cost*x + sint*z;
        mesh.nodes[i][2] = -sint*x + cost*z;
    }
}

////////////////////////////////////////////////////////////////////////////////

void zRotation(double angle)
{
    double cost = cos(M_PI*angle/180.0);
    double sint = sin(M_PI*angle/180.0);

    size_t numNodes = mesh.nodes.size();
    for( size_t i = 0; i < numNodes; i++) {
        float x = mesh.nodes[i][0];
        float y = mesh.nodes[i][1];
        mesh.nodes[i][0] =  cost*x - sint*y;
        mesh.nodes[i][1] =  sint*x + cost*y;
    }
}

////////////////////////////////////////////////////////////////////////////////

void drawCube(const Array3F &origin, const Array3F &length)
{
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

////////////////////////////////////////////////////////////////////////////////

void setSlice()
{
    if( slicePos > maxCorner[2] ) {
        slicePos = maxCorner[2];
        sliceID  = 0;
        sliceDir = -1;
    }

    if( slicePos < minCorner[2] ) {
        slicePos = minCorner[2];
        sliceID  = 0;
        sliceDir = +1;
    }

    cout << sliceID << "   SlicePos:  " << slicePos << endl;

    GLuint id = glGenLists(1);

    glNewList(id, GL_COMPILE);
    glPushMatrix();
    float xlen = maxCorner[0] - minCorner[0];
    float ylen = maxCorner[1] - minCorner[1];
    float zlen = (maxCorner[0] - minCorner[0])/(double)winHeight;
    Array3F origin = {minCorner[0], minCorner[1], slicePos};
    Array3F length = {xlen, ylen, zlen};
    drawCube(origin, length);
    glPopMatrix();
    glEndList();
    primitives[1] = new OpenCSG::DisplayListPrimitive(id, OpenCSG::Intersection, 1);

}

////////////////////////////////////////////////////////////////////////////////

void setModel()
{
    clearPrimitives();
    primitives.resize(2);
    setMesh();
    setSlice();
}

////////////////////////////////////////////////////////////////////////////////
void drawBox()
{
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    Array3F len;
    len[0] = maxCorner[0] - minCorner[0];
    len[1] = maxCorner[1] - minCorner[1];
    len[2] = maxCorner[2] - minCorner[2];
    glColor3f( 1.0, 1.0, 1.0);
    drawCube( minCorner, len);
}

////////////////////////////////////////////////////////////////////////////////

void drawAxis()
{
    glDisable (GL_BLEND);
    glDisable(GL_CULL_FACE);

    glLineWidth(5.0);
    glColor3f( 1.0, 0.0, 0.0);
    glBegin(GL_LINES);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(1000.0, 0.0, 0.0);
    glEnd();

    glColor3f( 0.0, 1.0, 0.0);
    glBegin(GL_LINES);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 1000.0, 0.0);
    glEnd();

    glColor3f( 0.0, 0.0, 1.0);
    glBegin(GL_LINES);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 1000.0);
    glEnd();
    glLineWidth(1.0);
}

////////////////////////////////////////////////////////////////////////////////

void display()
{
    eyePos = max( maxCorner[2] + 0.001F, eyePos);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0.0, 0.0, eyePos,    /* eye is at (0,2,5) */
              0.0, 0.0, 0.0,   /* center is at (0,0,0) */
              0.0, 1.0, 0.0);  /* up is in positive Y direction */
    glRotatef(rot, 0.0f, 1.0f, 0.0f);

    glDisable(GL_LIGHTING);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    OpenCSG::render(primitives);

    if( displayIntersection) {
        glColor3f( 1.0, 0.0, 0.0);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glDepthFunc(GL_EQUAL);
        if(primitives.size()  == 2) {
            primitives[0]->render();
            primitives[1]->render();
        }
        glDepthFunc(GL_LESS);
    }

    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glEnable(GL_CULL_FACE);

    if(primitives.size()  == 2) {
        if( displayMesh) {
            if( displayWireframe)
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            else
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glColor4f( 0.9, 0.9, 0.9, 0.5);
            primitives[0]->render();
        }

        if( displayPlane) {
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glColor3f( 0.0, 1.0, 1.0);
            primitives[1]->render();

            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glLineWidth(2.0);
            glColor3f( 0.0, 0.0, 1.0);
            primitives[1]->render();
            glLineWidth(1.0);
        }
    }

    if( displayBox) {
        glDisable(GL_CULL_FACE);
        drawBox();
    }

    if( displayAxis) drawAxis();

    glutSwapBuffers();
}

////////////////////////////////////////////////////////////////////////////////

int saveAs(const string &filename, const std::vector<unsigned char> &image, const std::array<int,2> &dim)
{
    FILE *f = fopen ( filename.c_str(), "wb" );
    if ( !f ) {
        cout << "Warning: Can not  open file " << endl;
        return 1;
    }
    cout << "Saving Image : " << filename << endl;
    int cols = dim[0];
    int rows = dim[1];

    fprintf(f, "P5\n" );
    fprintf(f, "%d %d\n", cols, rows );
    fprintf(f, "255\n" );

    vector<unsigned char> data(cols*rows);

    for( int j = 0; j < rows; j++) {
        for( int i = 0; i < cols; i++) {
            data[(rows-j-1)*cols + i] = image[j*cols+i];
        }
    }

    fwrite( &data[0], rows*cols, 1, f);

    fclose(f);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int saveScreenShots()
{
    int numWidth;
    if( numSlices < 1000000) numWidth  = 5;
    if( numSlices < 100000) numWidth  = 5;
    if( numSlices < 10000)  numWidth  = 4;
    if( numSlices < 1000)   numWidth  = 3;
    if( numSlices < 100)    numWidth  = 2;
    if( numSlices < 10)     numWidth  = 1;

    int texWidth  = frameRatio*winWidth;
    int texHeight = frameRatio*winHeight;

    int    fboStatus = 0;
    GLuint frameBuffer = 0;
    GLuint colorRenderBuffer;
    GLuint depthRenderBuffer;
    GLuint texture;

    glGenFramebuffers(1, &frameBuffer);
    glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer);

    glGenRenderbuffers(1, &colorRenderBuffer);
    glBindRenderbuffer(GL_RENDERBUFFER, colorRenderBuffer);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_RGB, texWidth, texHeight);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, colorRenderBuffer);

    glGenRenderbuffers(1, &depthRenderBuffer);
    glBindRenderbuffer(GL_RENDERBUFFER, depthRenderBuffer);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT16, texWidth, texHeight);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthRenderBuffer);

    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
//  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB,  texWidth, texHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texture, 0);
//  glEnable(GL_TEXTURE_2D);

    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, texture, 0);
    GLenum drawBuffers[2] = {GL_COLOR_ATTACHMENT0, GL_DEPTH_ATTACHMENT};
    glDrawBuffers(1, drawBuffers); // "1" is the size of DrawBuffers
    fboStatus = 1;
    if( glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
        fboStatus = 0;
    else
        cout << "FrameBuffer is ready to use" << endl;

    glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer);

    glViewport(0,0, texWidth, texHeight); // Render on the whole framebuffer, complete from the lower left corner to the upper right

    zoom  = 0.01;
    cameraType = 1;
    setCamera();
    displayPlane = 0;
    displayMesh  = 0;
    displayBox   = 0;
    displayAxis  = 0;

    std::array<int,2> dim = {texWidth, texHeight};

    size_t numPixels = texWidth*texHeight;
    vector<unsigned char> rgb(3*numPixels);
    vector<unsigned char> data(numPixels);

    slicePos  = minCorner[2];
    sliceDir  = 1;
    size_t index = 0;
    while(1) {
        if( sliceDir == -1) break;
        setModel();
        display();
        glReadBuffer(GL_FRONT);
        glReadPixels(0, 0, texWidth, texHeight, GL_RGB, GL_UNSIGNED_BYTE, &rgb[0] );
        for( size_t i = 0; i < numPixels; i++) data[i] = rgb[3*i];
        ostringstream oss;
        oss << "Slice" << setw(numWidth) << setfill('0') << index++ << ".pgm";
        filename = oss.str();
        saveAs(filename, data, dim);
        slicePos += sliceStep;
    }

    if( fboStatus) glDeleteFramebuffers( 1, &frameBuffer);

    glDisable(GL_TEXTURE_2D);

    glViewport(0,0, winWidth, winHeight); // Render on the whole framebuffer, complete from the lower left corner to the upper right

    cameraType = 0;
    setCamera();
    displayPlane = 1;
    displayMesh = 1;
    displayBox  = 1;
    slicePos  =  minCorner[2];
    setModel();
    display();
}

////////////////////////////////////////////////////////////////////////////////

void idle() {

    static int ancient = 0;
    static int last = 0;
    static int msec = 0;
    last = msec;
    msec = glutGet(GLUT_ELAPSED_TIME);

    if (spin) {
        rot += (msec-last)/100.0f;
        while (rot >= 360.0f)
            rot -= 360.0f;
    }

    if( animate ) {
        slicePos = slicePos + sliceDir*sliceStep;
        setModel();
    }

    display();
}

////////////////////////////////////////////////////////////////////////////////

void key(unsigned char k, int, int) {
    switch (k) {
    case '+':
        sliceDir = 1;
        slicePos += sliceStep;
        spin     = 0;
        zoom    += 0.01;
        sliceID++;
        setModel();
        break;
    case '-':
        sliceDir = -1;
        slicePos -= sliceStep;
        spin     = 0;
        zoom    -= 0.01;
        sliceID--;
        setModel();
        break;
    case 'o':
        spin = !spin;
        break;
    case 'b':
        displayBox = !displayBox;
        break;
    case 'm':
        displayMesh = !displayMesh;
        break;
    case 'p':
        displayPlane = !displayPlane;
        break;
    case 'i':
        displayIntersection = !displayIntersection;
        break;
    case 'n':
        slicePos =  slicePos + sliceDir*sliceStep;
        sliceID  =  sliceID  + sliceDir;
        spin     = 0;
        setModel();
        break;
    case 'r':
        rot  = 0.0;
        zoom = 0.0;
        eyePos = 2.0;
        setCamera();
        break;
    case 'w':
        displayWireframe = !displayWireframe;
        break;
    case 'c':
        cameraType = !cameraType;
        setCamera();
        break;
    case 'a':
        animate = !animate;
        break;
    case 'x':
        xRotation(90.0);
        initMesh();
        setModel();
        break;
    case 'y':
        yRotation(90.0);
        initMesh();
        setModel();
        break;
    case 'z':
        zRotation(90.0);
        initMesh();
        setModel();
        break;
    case 's':
        saveScreenShots();
        break;
    case 't':
        displayAxis = !displayAxis;
        break;
    default:
        break;
    }
    display();
}

////////////////////////////////////////////////////////////////////////////////

void SpecialKeys(int k, int, int) {

    switch (k) {
    case GLUT_KEY_UP:
        eyePos += 0.1;
        if( cameraType ) {
            zoom += 0.01;
            setCamera();
        }
        break;
    case GLUT_KEY_DOWN:
        eyePos -= 0.1;
        if( cameraType ) {
            zoom -= 0.01;
            setCamera();
        }
        break;
    case GLUT_KEY_LEFT:
        rot -= 1.0;
        break;
    case GLUT_KEY_RIGHT:
        rot += 1.0;
        break;
    case GLUT_KEY_PAGE_UP:
        slicePos =  maxCorner[2];
        sliceDir =  -1;
        setModel();
        break;
    case GLUT_KEY_PAGE_DOWN:
        slicePos =  minCorner[2];
        sliceDir =  +1;
        setModel();
        break;
    }

    display();
}

////////////////////////////////////////////////////////////////////////////////

void menu(int value) {
    switch (value) {
    case ALGO_AUTOMATIC:
        OpenCSG::setOption(OpenCSG::AlgorithmSetting, OpenCSG::Automatic);
        break;
    case GF_STANDARD:
        OpenCSG::setOption(OpenCSG::AlgorithmSetting, OpenCSG::Goldfeather);
        OpenCSG::setOption(OpenCSG::DepthComplexitySetting, OpenCSG::NoDepthComplexitySampling);
        break;
    case GF_DC:
        OpenCSG::setOption(OpenCSG::AlgorithmSetting, OpenCSG::Goldfeather);
        OpenCSG::setOption(OpenCSG::DepthComplexitySetting, OpenCSG::DepthComplexitySampling);
        break;
    case GF_OQ:
        OpenCSG::setOption(OpenCSG::AlgorithmSetting, OpenCSG::Goldfeather);
        OpenCSG::setOption(OpenCSG::DepthComplexitySetting, OpenCSG::OcclusionQuery);
        break;
    case SCS_STANDARD:
        OpenCSG::setOption(OpenCSG::AlgorithmSetting, OpenCSG::SCS);
        OpenCSG::setOption(OpenCSG::DepthComplexitySetting, OpenCSG::NoDepthComplexitySampling);
        break;
    case SCS_DC:
        OpenCSG::setOption(OpenCSG::AlgorithmSetting, OpenCSG::SCS);
        OpenCSG::setOption(OpenCSG::DepthComplexitySetting, OpenCSG::DepthComplexitySampling);
        break;
    case SCS_OQ:
        OpenCSG::setOption(OpenCSG::AlgorithmSetting, OpenCSG::SCS);
        OpenCSG::setOption(OpenCSG::DepthComplexitySetting, OpenCSG::OcclusionQuery);
        break;
    case OFFSCREEN_AUTOMATIC:
        OpenCSG::setOption(OpenCSG::OffscreenSetting, OpenCSG::AutomaticOffscreenType);
        break;
    case OFFSCREEN_FBO:
        OpenCSG::setOption(OpenCSG::OffscreenSetting, OpenCSG::FrameBufferObject);
        break;
    case OFFSCREEN_PBUFFER:
        OpenCSG::setOption(OpenCSG::OffscreenSetting, OpenCSG::PBuffer);
        break;
    default:
        break;
    }
    display();
}

////////////////////////////////////////////////////////////////////////////////
void initLights()
{
    GLfloat light_diffuse[]   = { 1.0f,  0.0f,  0.0f,  1.0f};  // Red diffuse light
    GLfloat light_position0[] = {-1.0f, -1.0f, -1.0f,  0.0f};  // Infinite light location
    GLfloat light_position1[] = { 1.0f,  1.0f,  1.0f,  0.0f};  // Infinite light location

    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position0);
    glEnable(GL_LIGHT0);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT1, GL_POSITION, light_position1);
    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHTING);
}

////////////////////////////////////////////////////////////////////////////////

void init()
{
    // gray background
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

//  initLights();

    glEnable(GL_NORMALIZE);

    // Use depth buffering for hidden surface elimination
    glEnable(GL_DEPTH_TEST);

    setCamera();

    readMesh();
    initMesh();

    setModel();

    OpenCSG::setOption(OpenCSG::AlgorithmSetting, OpenCSG::Goldfeather);
    OpenCSG::setOption(OpenCSG::DepthComplexitySetting, OpenCSG::NoDepthComplexitySampling);

    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
//  glDepthFunc(GL_LEQUAL);

}

////////////////////////////////////////////////////////////////////////////////

void Usage( const string &exe)
{
    cout << exe << " meshfile numslices winWidth frameRatio" << endl;

}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    if( argc != 5) {
        Usage(argv[0] );
        return 1;
    }

    filename  = argv[1];
    numSlices = atoi(argv[2]);
    winHeight = atoi(argv[3] );
    winWidth  = winHeight;
    frameRatio = atof( argv[4] );

    glutInit(&argc, argv);

    glutInitWindowSize(winWidth, winHeight);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL);
    glutCreateWindow("OpenCSG example application");
    glutSpecialFunc(SpecialKeys);

    int err = glewInit();
    if (GLEW_OK != err) {
        // problem: glewInit failed, something is seriously wrong
        std::cerr << "GLEW Error: " << glewGetErrorString(err) << std::endl;
        return 1;
    }

    int menuAlgorithm = glutCreateMenu(menu);
    glutAddMenuEntry("Goldfeather standard",GF_STANDARD);
    glutAddMenuEntry("Goldfeather depth complexity sampling", GF_DC);
    glutAddMenuEntry("Goldfeather occlusion query", GF_OQ);
    glutAddMenuEntry("SCS standard", SCS_STANDARD);
    glutAddMenuEntry("SCS depth complexity sampling", SCS_DC);
    glutAddMenuEntry("SCS occlusion query", SCS_OQ);
    glutAddMenuEntry("Automatic", ALGO_AUTOMATIC);

    int menuSettings = glutCreateMenu(menu);
    glutAddMenuEntry("Automatic", OFFSCREEN_AUTOMATIC);
    glutAddMenuEntry("Frame buffer object", OFFSCREEN_FBO);
    glutAddMenuEntry("PBuffer", OFFSCREEN_PBUFFER);

    glutCreateMenu(menu);
    glutAddSubMenu("CSG Algorithms", menuAlgorithm);
    glutAddSubMenu("Settings", menuSettings);

    glutAttachMenu(GLUT_RIGHT_BUTTON);

    glutDisplayFunc(display);
    glutKeyboardFunc(key);
    glutReshapeFunc(setScreen);

    menu(OFFSCREEN_AUTOMATIC);

    glutIdleFunc(idle);
    init();
    glutMainLoop();

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

