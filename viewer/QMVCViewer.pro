EXT_DIR = /home/boubek/lib/libQGLViewer-2.6.3 # YOU NEED TO SET YOUR EXT_DIR CORRECTLY, AND PUT libQGlViewer (v2.6.1 or more recent)

VIEWERSHARED_DIR = shared
GEOMSHARED_DIR = ../coordinates/shared

TEMPLATE = app
TARGET = QMVCViewer
DEPENDPATH += .
INCLUDEPATH += .
INCLUDEPATH += $${VIEWERSHARED_DIR}
INCLUDEPATH += $${GEOMSHARED_DIR}
QT += opengl xml
CONFIG += qt \
    release \
    c++11
MOC_DIR = ./tmp/moc
OBJECTS_DIR = ./tmp/obj



# DEPENDENCIES :

#QGLViewer
{
 INCLUDEPATH += $${EXT_DIR}
 LIBS +=    -L$${EXT_DIR}/QGLViewer \
            -lQGLViewer
}

LIBS += -lglut \
    -lGLU
LIBS += -lgsl \
    -lgomp
LIBS += -lblas \
    -lgomp



release:QMAKE_CXXFLAGS_RELEASE += -O3 \
    -fopenmp
release:QMAKE_CFLAGS_RELEASE += -O3 \
    -fopenmp

# Input
HEADERS += $${GEOMSHARED_DIR}/point3.h \
    CageManip.h \
    $${VIEWERSHARED_DIR}/qt/Manipulator.h \
    $${VIEWERSHARED_DIR}/qt/RectangleSelection.h \
    $${VIEWERSHARED_DIR}/gl/GLUtilityMethods.h \
    $${VIEWERSHARED_DIR}/gl/BasicColors.h \
    CageManipInterface.h \
    ../coordinates/qmvc/quadutilities.h \
    ../coordinates/mvc/mvc.h \
    ../coordinates/smvc/smvc.h \
    ../coordinates/qmvc/qmvc.h \
    ../coordinates/gc/gc.h \
    ../coordinates/mec/mec.h

SOURCES += main.cpp \
    $${VIEWERSHARED_DIR}/gl/GLUtilityMethods.cpp\
    $${VIEWERSHARED_DIR}/gl/BasicColors.cpp

