TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    data-structs/polymesh.cpp \
    helpers/spatial/spatial-algs.cpp \
    helpers/logger.cpp \
    polysgen/polysgen.cpp \
    core/relations.cpp \
    core/edge.cpp \
    core/polyhedral-set.cpp \
    core/surface/surface-edge.cpp \
    core/surface/surface-face.cpp \
    core/surface/front/surface-front-edge.cpp \
    core/surface/front/surface-front-vert.cpp \
    core/surface/surface.cpp \
    core/tetr.cpp \
    main.cpp \
    core/polyhedron.cpp \
    core/face.cpp \
    core/front/polyhedron-front-edge.cpp \
    core/front/polyhedron-front-face.cpp \
    data-structs/polyshell.cpp


HEADERS += \
    core/front/edge.h \
    core/front/face.h \
    core/shell/edge.h \
    core/shell/face.h \
    core/shell/front/edge.h \
    core/shell/front/vert.h \
    core/shell/vert.h \
    core/surface/edge.h \
    core/surface/face.h \
    core/surface/front/edge.h \
    core/surface/front/vert.h \
    core/surface/vert.h \
    data-structs/polymesh.h \
    helpers/spatial/algs.h \
    helpers/spatial/vec.h \
    helpers/logger.h \
    polysgen/polysgen.h \
    core/filetype.h \
    core/genparams.h \
    core/relations.h \
    core/edge.h \
    core/polyhedral-set.h \
    core/shell/shell.h \
    core/surface/surface.h \
    core/tetr.h \
    definitions.h \
    real-type.h \
    core/polyhedron.h \
    core/face.h \
    core/vert.h \
    data-structs/polyshell.h

