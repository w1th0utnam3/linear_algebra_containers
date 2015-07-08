TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

HEADERS += \
    ../source/matrix.h \
    ../source/column_vector.h \
    ../source/vector3.h \
    ../source/quaternion.h
    
