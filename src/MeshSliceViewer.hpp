///////////////////////////////////////////////////////////////////
//Written by  :  Chaman Singh Verma
//Date        :  1st Oct, 2017
//Description :  The following program visualizes slices extracted
//               from VDBSlicer.
//License     :  This program depends on "libQGLViewer" and "QT"
//               which are non-free for commerical usages. Use this
//               program only if you accept the above licenses.
///////////////////////////////////////////////////////////////////

#ifdef foreach
#undef foreach
#endif

#include <QGLViewer/qglviewer.h>
#include <qapplication.h>

#include <string>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QDialog>
#include <vector>
#include <set>
#include <cassert>

#include <opencsg.h>
#include "displaylistPrimitive.h"

enum {
    ALGO_AUTOMATIC, GF_STANDARD, GF_DC, GF_OQ, SCS_STANDARD, SCS_DC, SCS_OQ,
    OFFSCREEN_AUTOMATIC, OFFSCREEN_FBO, OFFSCREEN_PBUFFER
};

class MeshSliceViewer : public QGLViewer
{
    typedef std::array<float,3> Array3F;
    typedef std::array<int,2>   Array2I;
    typedef std::array<int,3>   Array3I;

    struct Mesh 
    {
        std::vector<Array3F> nodes;
        std::vector<Array2I> edges;
        std::vector<Array3I> faces;
    };



public:
    void readMesh(const std::string &file);

protected:
    virtual void draw();
    virtual void init();
    virtual void keyPressEvent( QKeyEvent *e);
    virtual void mouseMoveEvent(QMouseEvent *e);

private:
    Mesh  mesh;

    GLuint  meshDisplayID, sliceDisplayID;
    bool    mesh_initialized  = 0;
    bool    slice_initialized = 0;
    Array3F minCorner, maxCorner;

    int sliceDir = 2;
    int sliceID  = 0;

    std::vector<OpenCSG::Primitive*> primitives;

    void initMesh();

    void drawXPlane(int i);
    void drawYPlane(int i);
    void drawZPlane(int i);
    void selectPlane();
    void drawCube(const Array3F &origin, const Array3F &len);

    void clearPrimitives() {
        for (std::vector<OpenCSG::Primitive*>::const_iterator i = primitives.begin(); i != primitives.end(); ++i) {
            OpenCSG::DisplayListPrimitive* p =
                static_cast<OpenCSG::DisplayListPrimitive*>(*i);
            glDeleteLists(1, p->getDisplayListId());
            delete p;
        }

        primitives.clear();
    }




};
