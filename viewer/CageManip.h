/*
    QMVC - Reference Implementation of paper: 
    
    "Mean Value Coordinates for Quad Cages in 3D", 
    jean-Marc Thiery, Pooran Memari and Tamy Boubekeur
    SIGGRAPH Asia 2018
    
    This program allows to compute QMVC for a set of 3D points contained 
    in a cage made of quad and triangles, as well as other flavors of 
    space coordinates for cages (MVC, SMVC, GC, MEC). It comes also with 
    a 3D viewer which helps deforming a mesh with a cage. 
    
    Copyright (C) 2018  jean-Marc Thiery, Pooran Memari and Tamy Boubekeur

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
#ifndef CAGEMANIP_H
#define CAGEMANIP_H

#include <vector>
#include <queue>
#include <map>
#include <set>
#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cfloat>
#include <cmath>
#include <chrono>
#include <cassert>
#include <QCoreApplication>

#include "point3.h"
#include "qt/RectangleSelection.h"
#include "qt/Manipulator.h"
#include "gl/TextureHandler.h"

#define GL_GLEXT_PROTOTYPES
#include <GL/glew.h>
#include <gl/openglincludeQtComp.h>
#include <QOpenGLFunctions_4_3_Core>
#include <QOpenGLFunctions>
#include <QGLViewer/qglviewer.h>

#include "CageManipInterface.h"

#include <QTimer>
#include <qt/QSmartAction.h>
#include <QComboBox>
#include <QFormLayout>
#include <QToolBar>
#include <QColorDialog>
#include <QInputDialog>
#include <QFileDialog>
#include <QLabel>

#include <ctime>
#include <ratio>
#include <chrono>




struct AnimationTimer{
    bool playAnimation;
    double animationPeriodInMicroSeconds;
    double uAnimation;
    std::chrono::time_point<std::chrono::high_resolution_clock> timerPrevious;
    bool loopAnimation;

    AnimationTimer() {
        playAnimation = false;
        animationPeriodInMicroSeconds = 100000000.0;
        uAnimation = 0.0;
    }

    void setLoopAnimation(bool l) {
        loopAnimation = l;
    }

    void setAnimationLengthInMilliSeconds( double totalTimeInMilliSeconds ) {
        animationPeriodInMicroSeconds = 1000 * totalTimeInMilliSeconds;
    }

    bool togglePause() {
        playAnimation = !playAnimation;
        if(playAnimation){
            timerPrevious = std::chrono::high_resolution_clock::now();
        }
        return playAnimation;
    }

    void restart( double u = 0.0 ) {
        timerPrevious = std::chrono::high_resolution_clock::now();
        uAnimation = u;
    }

    void recomputeU(  ) {
        if( playAnimation ) {
            auto timerCurrent = std::chrono::high_resolution_clock::now();
            auto elapsed = timerCurrent - timerPrevious;
            long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
            double uAnim = uAnimation + ((double)(microseconds)/animationPeriodInMicroSeconds);
            if(loopAnimation)
                uAnim = uAnim - floor(uAnim);
            timerPrevious = timerCurrent;
            uAnimation = uAnim;
        }
    }

    bool isPlaying() {
        return playAnimation;
    }

    void multiplySpeedBy( double ffactor ) {
        animationPeriodInMicroSeconds /= ffactor;
    }

    double getU() {
        return uAnimation;
    }
};


static constexpr GLenum faceTarget[6] = {
    GL_TEXTURE_CUBE_MAP_POSITIVE_X_EXT,
    GL_TEXTURE_CUBE_MAP_NEGATIVE_X_EXT,
    GL_TEXTURE_CUBE_MAP_POSITIVE_Z_EXT,
    GL_TEXTURE_CUBE_MAP_NEGATIVE_Z_EXT,
    GL_TEXTURE_CUBE_MAP_POSITIVE_Y_EXT,
    GL_TEXTURE_CUBE_MAP_NEGATIVE_Y_EXT
};

class CMViewer : public QGLViewer , public QOpenGLFunctions_4_3_Core
{
    Q_OBJECT

    CMInterface< point3d > cm_interface;

    QWidget * controls;
    int nbSelectionTools;
    QComboBox *selectionToolCombo;


    // Selection :
    int selectionTool;
    RectangleSelection * _selection_rectangle;

    // Manipulator :
    SimpleManipulator * manipulator;

    // texture:
    ScalarTextureHandler * scalarTexture;

    int coordinate_id_to_show;


    AnimationTimer animationTimer;
    int m_video_lowAngleRevolutionDegrees;
    int m_video_degre_Start;
    int m_video_degre_End;

    // we revolve arount the point recorded by qglviewer (the user can change it easily):
    qglviewer::Vec revolvePoint;
    qglviewer::Vec position;
    double distanceToRevolvePoint;

    // record the camera orientation, so that we can restore it afterwards:
    qglviewer::Vec upVector;
    qglviewer::Vec rightVector;
    qglviewer::Vec viewDirection;

    bool envMappingMode;

public :

    inline QWidget * getControlWidget () { return controls; }

    //--------------------------------
    // Environnement Map Settings.
    //--------------------------------

    /* In case your <GL/gl.h> does not advertise EXT_texture_lod_bias... */
#ifndef GL_EXT_texture_lod_bias
# define GL_MAX_TEXTURE_LOD_BIAS_EXT         0x84fd
# define GL_TEXTURE_FILTER_CONTROL_EXT       0x8500
# define GL_TEXTURE_LOD_BIAS_EXT             0x8501
#endif

    /* In case your <GL/gl.h> does not advertise EXT_texture_cube_map... */
#ifndef GL_EXT_texture_cube_map
# define GL_NORMAL_MAP_EXT                   0x8511
# define GL_REFLECTION_MAP_EXT               0x8512
# define GL_TEXTURE_CUBE_MAP_EXT             0x8513
# define GL_TEXTURE_BINDING_CUBE_MAP_EXT     0x8514
# define GL_TEXTURE_CUBE_MAP_POSITIVE_X_EXT  0x8515
# define GL_TEXTURE_CUBE_MAP_NEGATIVE_X_EXT  0x8516
# define GL_TEXTURE_CUBE_MAP_POSITIVE_Y_EXT  0x8517
# define GL_TEXTURE_CUBE_MAP_NEGATIVE_Y_EXT  0x8518
# define GL_TEXTURE_CUBE_MAP_POSITIVE_Z_EXT  0x8519
# define GL_TEXTURE_CUBE_MAP_NEGATIVE_Z_EXT  0x851A
# define GL_PROXY_TEXTURE_CUBE_MAP_EXT       0x851B
# define GL_MAX_CUBE_MAP_TEXTURE_SIZE_EXT    0x851C
#endif


    static void loadCubeMapFace (const std::string & sideFilename,
                                 const std::string & topDownFilename) {
        QImage * sideImage = new QImage (QString (sideFilename.c_str ()));
        QImage * topDownImage = new QImage (QString (topDownFilename.c_str ()));
        for (unsigned int i = 0; i < 4; i++)
            gluBuild2DMipmaps (faceTarget[i], 4,
                               sideImage->width (), sideImage->height (),
                               GL_RGBA, GL_UNSIGNED_BYTE, sideImage->bits ());
        for (unsigned int i = 4; i < 6; i++)
            gluBuild2DMipmaps (faceTarget[i], 4,
                               topDownImage->width (), topDownImage->height (),
                               GL_RGBA, GL_UNSIGNED_BYTE, topDownImage->bits ());
    }

    void initEnvMap () {
        int cubeMapMode = GL_REFLECTION_MAP_EXT;
        int cubeMapWrap = GL_REPEAT;
        std::string dirPath = std::string (qPrintable (QCoreApplication::applicationDirPath ()));
        loadCubeMapFace (dirPath + "/images/cm_lines.png", dirPath + "/images/cm_grid.png");
        glTexParameteri (GL_TEXTURE_CUBE_MAP_EXT, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri (GL_TEXTURE_CUBE_MAP_EXT, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glEnable(GL_TEXTURE_CUBE_MAP_EXT);
        assert (cubeMapMode == GL_NORMAL_MAP_EXT || cubeMapMode == GL_REFLECTION_MAP_EXT);
        glTexGeni (GL_S, GL_TEXTURE_GEN_MODE, cubeMapMode);
        glTexGeni (GL_T, GL_TEXTURE_GEN_MODE, cubeMapMode);
        glTexGeni (GL_R, GL_TEXTURE_GEN_MODE, cubeMapMode);
        glTexParameteri (GL_TEXTURE_CUBE_MAP_EXT, GL_TEXTURE_WRAP_S, cubeMapWrap);
        glTexParameteri (GL_TEXTURE_CUBE_MAP_EXT, GL_TEXTURE_WRAP_T, cubeMapWrap);
        glDisable (GL_TEXTURE_CUBE_MAP_EXT);
    }

    CMViewer(QGLWidget * parent = NULL) : QGLViewer(parent) , QOpenGLFunctions_4_3_Core()
    {
        _selection_rectangle = new RectangleSelection;

        connect(_selection_rectangle,SIGNAL(add(QRectF const &)),this,SLOT(selectSlot(QRectF const &)));
        connect(_selection_rectangle,SIGNAL(remove(QRectF const &)),this,SLOT(unselectSlot(QRectF const &)));
        connect(_selection_rectangle,SIGNAL(apply()),this,SLOT(computeManipulatorForSelectionSlot()));

        manipulator = new SimpleManipulator;

        connect( manipulator , SIGNAL(moved()) , this , SLOT(cageChanged()) );
        connect( manipulator , SIGNAL(mouseReleased()) , this , SLOT(manipulatorReleased()) );

        selectionTool = 0;

        coordinate_id_to_show = -1;

        setWindowTitle(QString("QMVC Viewer. Press 'H' for help!"));

        m_video_lowAngleRevolutionDegrees = 10;
        m_video_degre_Start = 0;
        m_video_degre_End = 480;
        animationTimer.setAnimationLengthInMilliSeconds(10000);

        envMappingMode = false;

        //
        construct_controls();
    }



    void construct_controls()
    {
        controls = new QWidget( this );
        QFormLayout * layout = new QFormLayout;
        controls->setLayout( layout );
        controls->setWindowFlags(Qt::WindowStaysOnTopHint | Qt::Tool);
        controls->setWindowTitle( "QMVC Controls" );
        controls->setFocusPolicy( Qt::ClickFocus );

        DetailedAction * open_mesh = new DetailedAction( QIcon("./icons/open.png") , "Open Mesh" , "Open Mesh" , this , this , SLOT(open_mesh()) );
        DetailedAction * open_cage = new DetailedAction( QIcon("./icons/open-cage.png") , "Open Binding Cage" , "Open Binding Cage" , this , this , SLOT(open_cage()) );
        // DetailedAction * computeMEC = new DetailedAction( QIcon("./icons/work.png") , "Clamp and compute MEC" , "Clamp and compute MEC" , this , this , SLOT(computeMEC()) );
        DetailedAction * open_deformed_cage = new DetailedAction( QIcon("./icons/open-deformed-cage.png") , "Open Deformed Cage" , "Open Deformed Cage" , this , this , SLOT(open_deformed_cage()) );
        DetailedAction * saveDeformedModel = new DetailedAction( QIcon("./icons/save.png") , "Save mesh" , "Save Mesh" , this , this , SLOT(saveDeformedModel()) );
        DetailedAction * saveDeformedCage = new DetailedAction( QIcon("./icons/save-cage.png") , "Save Cage" , "Save Cage" , this , this , SLOT(saveDeformedCage()) );
        DetailedAction * help = new DetailedAction( QIcon("./icons/help.png") , "HELP" , "HELP" , this , this , SLOT(help()) );

        DetailedAction * saveCamera = new DetailedAction( QIcon("./icons/save-camera.png") , "Save camera" , "Save camera" , this , this , SLOT(saveCamera()) );
        DetailedAction * openCamera = new DetailedAction( QIcon("./icons/open-camera.png") , "Open camera" , "Open camera" , this , this , SLOT(openCamera()) );
        DetailedAction * saveSnapShotPlusPlus = new DetailedAction( QIcon("./icons/save-snapshot.png") , "Save snapshot" , "Save snapshot" , this , this , SLOT(saveSnapShotPlusPlus()) );

        DetailedAction * quit = new DetailedAction( QIcon("./icons/quit.png") , "Quit" , "Quit" , this , qApp, SLOT (closeAllWindows()));


        selectionToolCombo = new QComboBox;
        selectionToolCombo->addItem(QIcon("./icons/select-rect.svg"),"",QVariant(0));
        selectionToolCombo->addItem(QIcon("./icons/select-face.svg"),"",QVariant(1));
        selectionToolCombo->setCurrentIndex(0);
        connect( selectionToolCombo , SIGNAL(currentIndexChanged(int)) , this , SLOT(selectionToolChanged(int)) );
        nbSelectionTools = 2;

        QComboBox *updateMethod = new QComboBox;
#ifdef ALLOW_TRI_MVC
        updateMethod->addItem(tr("MVC"));
#endif
#ifdef ALLOW_QMVC
        updateMethod->addItem(tr("QMVC"));
#endif
#ifdef ALLOW_QMVC_MEC
        updateMethod->addItem(tr("QMVC_MEC"));
#endif
#ifdef ALLOW_SMVC
        updateMethod->addItem(tr("SMVC"));
#endif
#ifdef ALLOW_GC
        updateMethod->addItem(tr("GC"));
#endif
        updateMethod->setCurrentIndex(0);
        setMethod( updateMethod->currentIndex() );
        connect( updateMethod , SIGNAL(currentIndexChanged(int)) , this , SLOT(methodChanged(int)) );

        QComboBox *updateMode = new QComboBox;
        updateMode->addItem(tr("Interactive"));
        updateMode->addItem(tr("RealTime"));
        updateMode->setCurrentIndex(1);
        setMode( updateMode->currentIndex() );
        updateMode->setToolTip(QString("vertex update mode"));
        connect( updateMode , SIGNAL(currentIndexChanged(int)) , this , SLOT(modeChanged(int)) );
        QComboBox *normalupdateMode = new QComboBox;
        normalupdateMode->addItem(tr("None"));
        normalupdateMode->addItem(tr("Interactive"));
        normalupdateMode->addItem(tr("RealTime"));
        normalupdateMode->setCurrentIndex(2);
        normalupdateModeChanged( normalupdateMode->currentIndex() );
        normalupdateMode->setToolTip(QString("vertex normals update mode"));
        connect( normalupdateMode , SIGNAL(currentIndexChanged(int)) , this , SLOT(normalupdateModeChanged(int)) );

        // Add them :
        QToolBar *toolBar = new QToolBar;
        toolBar->addAction (open_mesh);
        toolBar->addAction (open_cage);
//        toolBar->addAction( computeMEC );
        toolBar->addAction (open_deformed_cage);
        toolBar->addAction (saveDeformedModel);
        toolBar->addAction (saveDeformedCage);
        toolBar->addSeparator ();
        toolBar->addAction (openCamera);
        toolBar->addAction (saveCamera);
        toolBar->addAction (saveSnapShotPlusPlus);
        toolBar->addSeparator ();
        toolBar->addWidget (new QLabel (tr (" Selection: "), this));
        toolBar->addWidget (selectionToolCombo);
        toolBar->addSeparator ();
        toolBar->addWidget (new QLabel (tr (" Coordinates: "), this));
        toolBar->addWidget (updateMethod);
        toolBar->addWidget (new QLabel (tr (" Vertex update: "), this));
        toolBar->addWidget (updateMode);
        toolBar->addWidget (new QLabel (tr (" Normal update: "), this));
        toolBar->addWidget (normalupdateMode);
        toolBar->addSeparator ();
        toolBar->addAction (help);
        toolBar->addSeparator ();
        toolBar->addAction (quit);
        toolBar->setIconSize(QSize(75, 75));
        layout->addRow( toolBar );
        /*QHBoxLayout * updatesBoxes = new QHBoxLayout;
        updatesBoxes->addWidget(updateMethod);
        updatesBoxes->addWidget(updateMode);
        updatesBoxes->addWidget(normalupdateMode);
        layout->addRow( updatesBoxes );
        */
        //controls->show ();
    }


    void playMovie( )
    {
        // record camera parameters:

        // we revolve arount the point recorded by qglviewer (the user can change it easily):
        revolvePoint = camera()->revolveAroundPoint();
        position = camera()->position();
        distanceToRevolvePoint = (revolvePoint - position).norm();

        // record the camera orientation, so that we can restore it afterwards:
        upVector = camera()->upVector();
        rightVector = camera()->rightVector();
        viewDirection = camera()->viewDirection();

        animationTimer.playAnimation = true;
        animationTimer.restart();

        startAnimation();
        setAnimationPeriod(10);
    }


    void draw()
    {
        if (envMappingMode == true) {
            glEnable (GL_TEXTURE_GEN_S);
            glEnable (GL_TEXTURE_GEN_T);
            glEnable (GL_TEXTURE_GEN_R);
            glEnable (GL_TEXTURE_CUBE_MAP_EXT);
        } else {
            glDisable (GL_TEXTURE_GEN_S);
            glDisable (GL_TEXTURE_GEN_T);
            glDisable (GL_TEXTURE_GEN_R);
            glDisable (GL_TEXTURE_CUBE_MAP_EXT);
        }

        // set camera parameters based on the timer :
        animationTimer.recomputeU();
        if( animationTimer.isPlaying() ) {
            double uAnim = animationTimer.getU();
            if(uAnim > 1.0) {
                animationTimer.playAnimation = false;
                camera()->setPosition( position );
                camera()->setViewDirection( viewDirection );
                camera()->setUpVector( upVector );
                camera()->computeModelViewMatrix();
                stopAnimation();
            }
            else{
                // find the corresponding rotation angle:
                double video_degre = (double)(m_video_degre_Start) + uAnim * (double)(m_video_degre_End - m_video_degre_Start);

                // camera transform:
                double localAngle = video_degre * 2.0 *  M_PI / 360.0;
                double cosLocalAngle = std::cos(localAngle);
                double sinLocalAngle = std::sin(localAngle);

                double lowAngleRevolution = m_video_lowAngleRevolutionDegrees * M_PI / 180.0;
                double cosLowAngleRevolution = cos(lowAngleRevolution);
                double sinLowAngleRevolution = sin(lowAngleRevolution);

                qglviewer::Vec localViewDirectionBeforeLowAngle = ( cosLocalAngle * viewDirection - sinLocalAngle * rightVector );
                qglviewer::Vec localViewDirection = cosLowAngleRevolution * localViewDirectionBeforeLowAngle - sinLowAngleRevolution * upVector;
                qglviewer::Vec localPosition = revolvePoint - distanceToRevolvePoint * localViewDirection;
                qglviewer::Vec localUpVector = sinLowAngleRevolution * localViewDirectionBeforeLowAngle + cosLowAngleRevolution * upVector;

                camera()->setPosition( localPosition );
                camera()->setViewDirection( localViewDirection );
                camera()->setUpVector( localUpVector );
                camera()->computeModelViewMatrix();
            }
        }

        glDisable(GL_LIGHTING);
        glColor3f(0,0,1);
        cm_interface.drawCageSelectedVertices();


//        glColor3f(1,0,0);
//        cm_interface.drawNANVertices();

        glEnable( GL_LIGHTING );
        if(coordinate_id_to_show < 0) {
            glColor3f( 0.8 , 0.8 , 0.8 );
            if (envMappingMode == true) {
                glEnable (GL_TEXTURE_GEN_S);
                glEnable (GL_TEXTURE_GEN_T);
                glEnable (GL_TEXTURE_GEN_R);
                glEnable (GL_TEXTURE_CUBE_MAP_EXT);
            }
            cm_interface.drawModifiedModel( );
            {
                glDisable (GL_TEXTURE_GEN_S);
                glDisable (GL_TEXTURE_GEN_T);
                glDisable (GL_TEXTURE_GEN_R);
                glDisable (GL_TEXTURE_CUBE_MAP_EXT);
            }
        }
        else {
            scalarTexture->enableTexturesAndBind();
            cm_interface.drawModifiedModelWithTexture( scalarTexture , coordinate_id_to_show );
            scalarTexture->unBindAndDisableTextures();
        }
        glDisable( GL_LIGHTING );

        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glColor3f(0.f,0.f,1.f);
        cm_interface.drawCage( 0 );

        glColor3f(0.f,1.f,0.f);
        cm_interface.drawcage_quads_cuts();
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(1,1);
        glEnable( GL_BLEND );
        glColor4f( 0.2 , 0.2 , 0.8 , 0.25 );
        cm_interface.drawCage( 10 );
        glDisable(GL_POLYGON_OFFSET_FILL);
        glDisable( GL_BLEND );

        if(coordinate_id_to_show >= 0 && coordinate_id_to_show < (int)(cm_interface.get_cage_vertices().size())) {
            glPointSize(6);
            glColor3f(1,0,0);
            glDisable(GL_LIGHTING);
        glBegin(GL_POINTS);
        point3d pCageVertex = cm_interface.get_cage_vertices()[coordinate_id_to_show];
        glVertex3f( pCageVertex[0] , pCageVertex[1] , pCageVertex[2] );
        glEnd();
        }

        _selection_rectangle->draw();

        glDisable( GL_DEPTH_TEST );
        manipulator->draw();
        glEnable( GL_DEPTH_TEST );
    }

    void pickBackgroundColor()
    {
        QColor _bc = QColorDialog::getColor( this->backgroundColor(), this);
        if( _bc.isValid() )
        {
            this->setBackgroundColor( _bc );
            this->updateGL();
        }
    }

    void adjustCamera( point3d const & bb , point3d const & BB )
    {
        point3d const & center = ( bb + BB )/2.f;
        setSceneCenter (qglviewer::Vec( center[0] , center[1] , center[2]));
        setSceneRadius ((BB - bb ).norm()/2.f);
        showEntireScene();
    }


    void showControls()
    {
        // Show controls :
        controls->close();
        controls->show();
    }



    void init()
    {
        makeCurrent();
        initializeOpenGLFunctions();

        setMouseTracking(true);// Needed for MouseGrabber.

        QColor initBackgroundColor (255,255,255);
        initBackgroundColor.setAlpha (255);
        setBackgroundColor(initBackgroundColor);

        // textures:
        scalarTexture = new ScalarTextureHandler;
        scalarTexture->initTexture(this);

        // Lights:
        GLTools::initLights();
        GLTools::setSunsetLight();
        GLTools::setDefaultMaterial();

        //
        glShadeModel(GL_SMOOTH);
        glFrontFace(GL_CCW); // CCW ou CW

        glEnable (GL_POINT_SMOOTH);
        glPointSize (8.f);
        glHint (GL_POINT_SMOOTH_HINT, GL_NICEST);
        glEnable (GL_LINE_SMOOTH);
        glLineWidth (2.5f);
        glHint (GL_LINE_SMOOTH_HINT, GL_NICEST);
        glHint (GL_POLYGON_SMOOTH_HINT, GL_NICEST);

        glEnable(GL_DEPTH);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LEQUAL);

        glEnable(GL_CLIP_PLANE0);

        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        glEnable(GL_COLOR_MATERIAL);

        context()->format().setSampleBuffers(true);

        initEnvMap ();

        //
        setSceneCenter( qglviewer::Vec( 0 , 0 , 0 ) );
        setSceneRadius( 10.f );
        showEntireScene();
    }


    QString helpString() const
    {
        QString text("<h2>QMVC Viewer</h2>");
        text += "<p>";
        text += "This application is Reference implementation of the research paper:<br>";
        text += "  <b>Mean value coordinates for quad cages in 3D</b><br>";
        text += "  <i>Jean-Marc Thiery, Pooran Memari and Tamy Boubekeur</i><br>";
        text += "  ACM Transactions on Graphics - Proc. SIGGRAPH Asia 2018<br>";
        text += "  <a href=\"https://www.telecom-paristech.fr/~boubek/papers/QMVC<br>\" Go to the project page</a><br>";
        text += "<br>";
        text += "<h3>Basics</h3>";
        text += "<p>";
        text += "<ul>";
        text += "<li>H :   make this help appear</li>";
        text += "<li>Ctrl + mouse right button double click :   choose background color</li>";
        text += "</ul>";
        text += "<ul>";
        text += "<li>Shift + mouse left: select cage vertices</li>";
        text += "<li>Shift + Ctrl + mouse left: unselect cage vertices</li>";
        text += "<li>Shift + mouse right click: manipulate the selected vertices</li>";
        text += "</ul>";
        text += "<h3>User guide</h3>";
        text += "<p>";
        text += "<b>FIRST</b> open a mesh, <b>THEN</b> open a cage. OBJ and OFF files are supported.";
        text += "</p>";
        text += "<p>";
        text += "When you load a cage, there will be a latence time due to the computation of the cage coordinates.";
        text += "</p>";
        text += "<p>";
        text += "To deform the mesh, use the 'Rectangle' selection tool to select cage vertices, by using the mouse while keeping 'Shift' pressed (discussed before).";
        text += "To disable the manipulation tool, right click on it.";
        text += "</p>";
        return text;
    }

    void keyPressEvent( QKeyEvent * event )
    {
        if( event->key() == Qt::Key_H )
        {
            help();
        }
        else if( event->key() == Qt::Key_Space )
        {
            playMovie();
        }
        else if( event->key() == Qt::Key_E )
        {
            envMappingMode = !envMappingMode;
            update();
            return;
        }
        else if( event->key() == Qt::Key_T  && ( event->modifiers() & Qt::CTRL ) )
        {
            bool ok;
            QString text = QInputDialog::getText(this, tr(""), tr("title:"), QLineEdit::Normal,this->windowTitle(), &ok);
            if (ok && !text.isEmpty())
            {
                this->setWindowTitle( text );
            }
        }
        else if( event->key() == Qt::Key_Tab )
        {
            selectionToolCombo->setCurrentIndex( (selectionToolCombo->currentIndex() + 1) % nbSelectionTools );
        }
        else if( event->key() == Qt::Key_Up )
        {
            if( (event->modifiers() & Qt::ControlModifier) )
                ++coordinate_id_to_show;
            update();
        }
        else if( event->key() == Qt::Key_Down )
        {
            if( (event->modifiers() & Qt::ControlModifier) )
                coordinate_id_to_show = std::max( coordinate_id_to_show-1 , -1 );
            update();
        }

    }


    void mouseDoubleClickEvent( QMouseEvent * e )
    {
        if( (e->modifiers() & Qt::ControlModifier)  &&  (e->button() == Qt::RightButton) )
        {
            pickBackgroundColor();
            return;
        }

        if( (e->modifiers() & Qt::ControlModifier)  &&  (e->button() == Qt::LeftButton) )
        {
            showControls();
            return;
        }

        QGLViewer::mouseDoubleClickEvent( e );
    }

    void mousePressEvent(QMouseEvent* e )
    {
        if( ( e->modifiers() & Qt::ShiftModifier ) )
        {
            if( selectionTool == 0 )
            {
                if( _selection_rectangle->isInactive() )
                {
                    _selection_rectangle->activate();
                }
                _selection_rectangle->mousePressEvent( e , camera() );
                updateGL();
                return;
            }
        }

        QGLViewer::mousePressEvent(e);
    }

    void mouseMoveEvent(QMouseEvent* e  )
    {
        if( ( e->modifiers() & Qt::ShiftModifier ) )
        {
            if( selectionTool == 0 )
            {
                if( _selection_rectangle->isInactive() )
                {
                    _selection_rectangle->activate();
                }
                _selection_rectangle->mouseMoveEvent( e , camera() );
                updateGL();
                return;
            }
        }

        QGLViewer::mouseMoveEvent(e);
    }

    void mouseReleaseEvent(QMouseEvent* e  )
    {
        if( ( e->modifiers() & Qt::ShiftModifier ) )
        {
            if( selectionTool == 0 )
            {
                if( _selection_rectangle->isInactive() )
                {
                    _selection_rectangle->activate();
                }
                _selection_rectangle->mouseReleaseEvent( e , camera() );
                updateGL();
                return;
            }
        }

        QGLViewer::mouseReleaseEvent(e);
    }



    void saveCameraInFile(const QString &filename){
        std::ofstream out (filename.toUtf8());
        if (!out)
            exit (EXIT_FAILURE);
        // << operator for point3 causes linking problem on windows
        out << camera()->position()[0] << " \t" << camera()->position()[1] << " \t" << camera()->position()[2] << " \t" " " <<
                                          camera()->viewDirection()[0] << " \t" << camera()->viewDirection()[1] << " \t" << camera()->viewDirection()[2] << " \t" << " " <<
                                          camera()->upVector()[0] << " \t" << camera()->upVector()[1] << " \t" <<camera()->upVector()[2] << " \t" <<" " <<
                                          camera()->fieldOfView();
        out << std::endl;

        out.close ();
    }

    void openCameraFromFile(const QString &filename){

        std::ifstream file;
        file.open(filename.toStdString().c_str());

        qglviewer::Vec pos;
        qglviewer::Vec view;
        qglviewer::Vec up;
        float fov;

        file >> (pos[0]) >> (pos[1]) >> (pos[2]) >>
                                                    (view[0]) >> (view[1]) >> (view[2]) >>
                                                                                           (up[0]) >> (up[1]) >> (up[2]) >>
                                                                                                                            fov;

        camera()->setPosition(pos);
        camera()->setViewDirection(view);
        camera()->setUpVector(up);
        camera()->setFieldOfView(fov);

        camera()->computeModelViewMatrix();
        camera()->computeProjectionMatrix();

        updateGL();
    }


public slots:

    void open_mesh()
    {
        point3d bb , BB;
        QString fileName = QFileDialog::getOpenFileName(this,
                                                        tr ("Choose a mesh to load"),
                                                        "",
                                                        "3D Surface Mesh Formats (*.obj *.off)");
        if ( !fileName.isNull() ) {                 // got a file name
            if(fileName.endsWith(QString(".off"))) {
                //BasicConvertor::OFF2OBJ(fileName.toStdString() , (fileName + QString(".obj")).toStdString());
                cm_interface.open_OFF_mesh(fileName.toStdString() , bb , BB);
            }
            else if(fileName.endsWith(QString(".obj"))) {
                //BasicConvertor::OBJ2OFF(fileName.toStdString() , (fileName + QString(".off")).toStdString());
                cm_interface.open_OBJ_mesh(fileName.toStdString() , bb , BB);
            }
            adjustCamera( bb , BB );
            manipulator->setDisplayScale( (BB - bb).norm() / 10.f );
            std::cout << "Bounding box: " << bb << "   -->  " << BB << std::endl;
            std::cout << cm_interface.get_mesh_vertices().size() << " mesh vertices" << std::endl;
            updateGL();
        }
    }

    void saveDeformedModel()
    {
        QString fileName = QFileDialog::getSaveFileName(this,
                                                        tr ("Choose a filename to save the deformed mesh"),
                                                        "",
                                                        "3D Surface Mesh Formats (*.obj *.off)");
        if ( !fileName.isNull() ) {                 // got a file name
            cm_interface.saveDeformedModel(fileName.toStdString());
        }
    }

    void saveDeformedCage()
    {
        QString fileName = QFileDialog::getSaveFileName(this,
                                                        tr ("Choose a filename to save the deformed cage"),
                                                        "",
                                                        "3D Surface Mesh Formats (*.obj *.off)");
        if ( !fileName.isNull() ) {                 // got a file name
            cm_interface.saveDeformedCage(fileName.toStdString());
        }
    }


    void openCamera(){
        QString fileName = QFileDialog::getOpenFileName(this,"","*.cam");
        if ( !fileName.isNull() ) {                 // got a file name
            openCameraFromFile(fileName);
        }
    }
    void saveCamera(){
        QString fileName = QFileDialog::getSaveFileName(this,"","*.cam");
        if ( !fileName.isNull() ) {                 // got a file name
            saveCameraInFile(fileName);
        }
    }

    void saveSnapShotPlusPlus(){
        QString fileName = QFileDialog::getSaveFileName(this,"*.png","");
        if ( !fileName.isNull() ) {                 // got a file name
            setSnapshotFormat("PNG");
            setSnapshotQuality(100);
            saveSnapshot( fileName );
            saveCameraInFile( fileName+QString(".cam") );
        }
    }


    void open_cage()
    {
        QString fileName = QFileDialog::getOpenFileName(this,
                                                        tr ("Choose a binding cage open"),
                                                        "",
                                                        "3D Surface Mesh Formats (*.obj *.off)");
        if ( !fileName.isNull() ) {                 // got a file name
            if(fileName.endsWith(QString(".off"))) {
                cm_interface.open_OFF_cage(fileName.toStdString());
                //BasicConvertor::OFF2OBJ(fileName.toStdString() , fileName.toStdString() + ".obj");
            }
            else if(fileName.endsWith(QString(".obj"))) {
                cm_interface.open_OBJ_cage(fileName.toStdString());
            }
            std::cout << cm_interface.get_cage_vertices().size() << " cage vertices" << std::endl;
            updateGL();
        }
    }
    void open_deformed_cage()
    {
        QString fileName = QFileDialog::getOpenFileName(this,
                                                        tr ("Choose a deformed cage open"),
                                                        "",
                                                        "3D Surface Mesh Formats (*.obj *.off)");
        if ( !fileName.isNull() ) {                 // got a file name
            if(fileName.endsWith(QString(".off"))) {
                std::vector< point3d > newVerts;
                OFFIO::open(fileName.toStdString() , newVerts );
                cm_interface.setCageDeformedVertices(newVerts);
            }
            else if(fileName.endsWith(QString(".obj"))) {
                std::vector< point3d > newVerts;
                OBJIO::open(fileName.toStdString() , newVerts );
                cm_interface.setCageDeformedVertices(newVerts);
            }
            updateGL();
        }
    }



    void selectSlot( QRectF const & zone )
    {
        cm_interface.select( zone );
    }


    void unselectSlot( QRectF const & zone )
    {
        cm_interface.unselect( zone );
    }



    void computeManipulatorForSelectionSlot()
    {
        cm_interface.computeManipulatorForSelection( manipulator );
    }


    void cageChanged()
    {
        cm_interface.cageChanged( manipulator );
    }



    void manipulatorReleased()
    {
        cm_interface.manipulatorReleased( );
    }



    void modeChanged( int m )
    {
        cm_interface.setMode( m );
    }
    void normalupdateModeChanged( int m )
    {
        cm_interface.setNormalUpdateMode( m );
    }
    void methodChanged( int m )
    {
        cm_interface.setMethod( m );
        cm_interface.updateModel();
        cm_interface.update_mesh_vertex_normals();
        updateGL();
    }


    void setMode( int m )
    {
        cm_interface.setMode( m );
    }
    void setMethod( int m )
    {
        cm_interface.setMethod( m );
    }


    void selectionToolChanged( int t )
    {
        selectionTool = t;
    }


//    void computeMEC() {
//        cm_interface.clampNegativeCoordinatesAndUseMEC();
//        cm_interface.updateModel();
//        cm_interface.update_mesh_vertex_normals();
//    }
};







#endif // CAGEMANIP_H
