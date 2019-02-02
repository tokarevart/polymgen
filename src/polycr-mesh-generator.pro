TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    Crystallite3.cpp \
    Edge3.cpp \
    Facet3.cpp \
    main.cpp \
    Polycrystal3.cpp \
    ShellEdge3.cpp \
    ShellFacet3.cpp \
    ShellVertex3.cpp \
    Simplex3.cpp \
    Vertex3.cpp \
    data-structures/PolyMesh.cpp \
    data-structures/PolyStruct.cpp \
    helpers/spatialalgs/SpatialAlgs.cpp \
    helpers/spatialalgs/Vec3.cpp \
    helpers/Logger.cpp \
    helpers/Timer.cpp \
    polygen/PolyGen.cpp

HEADERS += \
    Crystallite3.h \
    Edge3.h \
    Facet3.h \
    Polycrystal3.h \
    ShellEdge3.h \
    ShellFacet3.h \
    ShellVertex3.h \
    Simplex3.h \
    Vertex3.h \
    data-structures/PolyMesh.h \
    data-structures/PolyStruct.h \
    helpers/spatialalgs/SpatialAlgs.h \
    helpers/spatialalgs/Vec3.h \
    helpers/Logger.h \
    helpers/Timer.h \
    polygen/PolyGen.h
