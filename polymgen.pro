TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    src/data-structs/polymesh.cpp \
    src/helpers/spatial/spatial-algs.cpp \
    src/helpers/logger.cpp \
    src/polysgen/polysgen.cpp \
    src/core/relations.cpp \
    src/core/edge.cpp \
    src/core/polyhedral-set.cpp \
    src/core/surface/surface-edge.cpp \
    src/core/surface/surface-face.cpp \
    src/core/surface/front/surface-front-edge.cpp \
    src/core/surface/front/surface-front-vert.cpp \
    src/core/surface/surface.cpp \
    src/core/tetr.cpp \
    example.cpp \
    src/core/polyhedron.cpp \
    src/core/face.cpp \
    src/core/front/polyhedron-front-edge.cpp \
    src/core/front/polyhedron-front-face.cpp \
    src/data-structs/polyshell.cpp


HEADERS += \
    src/core/front/edge.h \
    src/core/front/face.h \
    src/core/shell/edge.h \
    src/core/shell/face.h \
    src/core/shell/front/edge.h \
    src/core/shell/front/vert.h \
    src/core/shell/vert.h \
    src/core/surface/edge.h \
    src/core/surface/face.h \
    src/core/surface/front/edge.h \
    src/core/surface/front/vert.h \
    src/core/surface/vert.h \
    src/data-structs/polymesh.h \
    src/helpers/mathconsts.h \
    src/helpers/spatial/algs.h \
    src/helpers/spatial/vec.h \
    src/helpers/logger.h \
    src/polysgen/polysgen.h \
    src/core/filetype.h \
    src/core/genparams.h \
    src/core/relations.h \
    src/core/edge.h \
    src/core/polyhedral-set.h \
    src/core/shell/shell.h \
    src/core/surface/surface.h \
    src/core/tetr.h \
    src/definitions.h \
    src/real-type.h \
    src/core/polyhedron.h \
    src/core/face.h \
    src/core/vert.h \
    src/data-structs/polyshell.h

