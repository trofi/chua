# Dependencies are:
# - pkg-config (pkg-config provider)
# - freeglut (glut provider)
# - mesa (GLU provider)
# - libglvnd (GL provider)

CFLAGS := $(shell pkg-config --cflags glut) $(shell pkg-config --cflags glu) $(shell pkg-config --cflags gl) -Wall -Wextra -O2
LDLIBS := $(shell pkg-config --libs   glut) $(shell pkg-config --libs   glu) $(shell pkg-config --libs   gl) -lm

all: chua
