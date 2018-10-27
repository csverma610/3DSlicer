OBJS = area.o batch.o channelManager.o context.o frameBufferObject.o frameBufferObjectExt.o occlusionQuery.o opencsgRender.o openglHelper.o pBufferTexture.o primitive.o primitiveHelper.o renderGoldfeather.o renderSCS.o scissorMemo.o settings.o stencilManager.o ../RenderTexture/RenderTexture.o displaylistPrimitive.o

OBJS += main.o

CPPFLAGS = -std=c++14 -O3 -fPIC
CPPFLAGS += -I$(GLEW_ROOT)/include -I../include  

#CPPFLAGS += -I$(QGLVIEWER_ROOT)/include
#CPPFLAGS += -I$(QTDIR)/include -I$(QTDIR)/include/QtCore -I$(QTDIR)/include/QtXml -I$(QTDIR)/include/QtWidgets -I$(QTDIR)/include/QtGui -I$(QTDIR)/include/QtOpenGL


LIBS  = -lGL -lGLU -lglut -L$(GLEW_ROOT)/lib -lGLEW -lX11
#LIBS += -L$(QGLVIEWER_ROOT)/lib -lQGLViewer
#LIBS += -L$(QTDIR)/lib -lQt5Widgets -lQt5Core -lQt5Gui -lQt5OpenGL


csg:$(OBJS)
	g++ -o csg $(OBJS) $(LIBS)


.o:.cpp
	g++ $(CPPFLAGS) $<

clean:
	\rm -rf *.o *.pgm csg
