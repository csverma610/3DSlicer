OBJS = area.o batch.o channelManager.o context.o frameBufferObject.o frameBufferObjectExt.o occlusionQuery.o opencsgRender.o openglHelper.o pBufferTexture.o primitive.o primitiveHelper.o renderGoldfeather.o renderSCS.o scissorMemo.o settings.o stencilManager.o ../RenderTexture/RenderTexture.o displaylistPrimitive.o

OBJS += main.o

CPPFLAGS = -std=c++14 -O3 -fPIC
CPPFLAGS += -I$(GLEW_DIR)/include -I../include  


LIBS  = -lGL -lGLU -lglut -L$(GLEW_DIR)/lib -lGLEW -lX11


csg:$(OBJS)
	g++ -o csg $(OBJS) $(LIBS)


.o:.cpp
	g++ $(CPPFLAGS) $<

clean:
	\rm -rf *.o *.pgm csg
