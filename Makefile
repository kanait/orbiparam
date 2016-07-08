CC = g++
TARGET = orbiparam

CXXFLAGS := -Wall -I../../megapoly/include/vecmath-c++-1.2-1.4 -I../../megapoly/render -I/usr/local/include -I/usr/local/Cellar/eigen/3.2.8/include/eigen3

RELEASE := 1
ifeq ($(RELEASE),1)
CXXFLAGS += -O3 -g
else
CXXFLAGS += -g
endif

SRC_DIR := .
SRC := orbiparam.cxx
INC := MyMesh.hxx Orbifold.hxx ShortestPathDijkstra.hxx ToMeshR.hxx VertexPQ.hxx
#INC += $(shell find ../render -name "*.hxx")
#INC += $(shell find ../render -name "*.h")

LIBS_PATH = 

LIBS = -L/usr/local/lib -lOpenMeshCore -lGLEW -lglfw3 -framework OpenGL
LIBS += -lpng -lz

#STATIC_LIBS = /usr/local/lib/libboost_filesystem.a
#STATIC_LIBS += /usr/local/lib/libboost_system.a
#STATIC_LIBS += /usr/local/lib/libboost_thread.a

all:
	@make $(TARGET)

depend: .depend

.depend: $(SRC)
	@echo "Making dependencies ..."
	@rm -f ./.depend
	@$(CC) $(CXXFLAGS) -MM $^ >> ./.depend

include .depend

#%.o : %.cxx
orbiparam.o : orbiparam.cxx
	@echo "Compiling $< ..."
	@$(CC) $(CXXFLAGS) -o $@ -c $<

$(TARGET): $(SRC:.cxx=.o)
	@echo "Linking $@..."
	@$(CC) -o $@ $(SRC:.cxx=.o) $(LIBS_PATH) $(LIBS) $(STATIC_LIBS)

.PHONY: clean
clean:
	\rm -f $(SRC_DIR)/*.o
	\rm -f $(TARGET)
	\rm -f .depend

