TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    data-structures/polymesh.cpp \
    data-structures/polystruct.cpp \
    helpers/spatial-algs/spatial-algs.cpp \
    helpers/spatial-algs/vec.cpp \
    helpers/logger.cpp \
    helpers/timer.cpp \
    polygen/polygen.cpp \
    spatial-objs/front/plane/front-plane-edge.cpp \
    spatial-objs/front/plane/front-plane-vertex.cpp \
    spatial-objs/front/surface/front-surface-edge.cpp \
    spatial-objs/shell/shell-edge.cpp \
    spatial-objs/shell/shell-vertex.cpp \
    spatial-objs/edge.cpp \
    spatial-objs/polyhedral-set.cpp \
    spatial-objs/tetr.cpp \
    spatial-objs/vertex.cpp \
    main.cpp \
    real-type.cpp \
    spatial-objs/polyhedron.cpp \
    spatial-objs/face.cpp \
    spatial-objs/shell/shell-face.cpp \
    spatial-objs/front/surface/front-surface-face.cpp


HEADERS += \
    data-structures/polymesh.h \
    data-structures/polystruct.h \
    helpers/spatial-algs/spatial-algs.h \
    helpers/spatial-algs/vec.h \
    helpers/logger.h \
    helpers/timer.h \
    polygen/polygen.h \
    spatial-objs/front/plane/front-plane-edge.h \
    spatial-objs/front/plane/front-plane-vertex.h \
    spatial-objs/front/surface/front-surface-edge.h \
    spatial-objs/shell/shell-edge.h \
    spatial-objs/shell/shell-vertex.h \
    spatial-objs/edge.h \
    spatial-objs/polyhedral-set.h \
    spatial-objs/tetr.h \
    spatial-objs/vertex.h \
    definitions.h \
    helpers/cosd-values.h \
    real-type.h \
    spatial-objs/polyhedron.h \
    spatial-objs/face.h \
    spatial-objs/shell/shell-face.h \
    spatial-objs/front/surface/front-surface-face.h

