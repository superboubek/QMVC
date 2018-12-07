TEMPLATE = app
TARGET = QMVCViewer

VIEWERSHARED_DIR = shared
GEOMSHARED_DIR = ../coordinates/shared
DEPENDPATH += .
INCLUDEPATH += .
INCLUDEPATH += $${VIEWERSHARED_DIR}
INCLUDEPATH += $${GEOMSHARED_DIR}
QT += opengl xml
CONFIG += qt release c++11
MOC_DIR = ./tmp/moc
OBJECTS_DIR = ./tmp/obj

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

# DEPENDENCIES :

unix {
 QGLV_DIR = /home/boubek/lib/libQGLViewer-2.6.3 # YOU NEED TO SET YOUR EXT_DIR CORRECTLY, AND PUT libQGlViewer (v2.6.1 or more recent)
 INCLUDEPATH += $${QGLV_DIR}
 LIBS +=    -L$${QGLV_DIR}/QGLViewer -lQGLViewer
 LIBS += -lglut -lGLU -lgsl -lgomp -lblas -lgomp
}

win32 {
    GLEW_DIR = 'c:/Dev/glew-2.0.0'
    INCLUDEPATH += '$$GLEW_DIR\include'
    LIBS += -L"$$GLEW_DIR/lib/Release\x64" -L"$$GLEW_DIR/bin/Release/x64" -lopengl32 -lglu32 -lglew32
    QGLV_DIR = 'C:\Dev\libQGLViewer-2.6.4'
    INCLUDEPATH += '$$QGLV_DIR'
    LIBS +=-L"$$QGLV_DIR\QGLViewer" -lQGLViewer2
    GNU_DIR =  'C:\Dev\GnuWin32'
    INCLUDEPATH += '$$GNU_DIR\include'
    LIBS += -L"$$GNU_DIR\msvc2015_64\lib\gsl" -lgsl -lcblas
    QMAKE_CXXFLAGS_WARN_ON -= -wd4251
}


release:QMAKE_CXXFLAGS_RELEASE += -O3 \
    -fopenmp
release:QMAKE_CFLAGS_RELEASE += -O3 \
    -fopenmp

win32:RC_ICONS = icons/qmvc.ico


