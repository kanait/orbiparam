#
# Copyright (c) Mark J. Kilgard, 1996, 1997.
#

LN = ln -s
MV = mv
RM = -rm -rf

TARGETS = param

CFLAGS = -fpermissive -I ../../Eigen/eigen-eigen-10219c95fe65 -I ../../OpenMesh/OpenMesh-3.3/src
LDLIBS = -L../../OpenMesh/OpenMesh-3.3/build_linux/Build/lib -lOpenMeshCore

SRCS = param.cxx
HEADERS =
OBJS = $(SRCS:.c=.o)

default : $(TARGETS)

param : $(OBJS)
	$(CC) -o $@ $(CFLAGS) $(OBJS) $(LDFLAGS) $(LDLIBS)

clean:
	$(RM) $(OBJS) $(TARGETS).exe *.stackdump
