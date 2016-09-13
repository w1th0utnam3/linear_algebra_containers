TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

INCLUDEPATH += ../src/

HEADERS += \
        ../src/matrix.h \
        ../src/column_vector.h \
        ../src/vector3.h \
        ../src/quaternion.h \
    ../src/matrixbase.h

