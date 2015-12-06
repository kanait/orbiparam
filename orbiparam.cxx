////////////////////////////////////////////////////////////////////
//
// $Id: $
//
// Copyright (c) 2004-2015 by Takashi Kanai. All rights reserved.
//
////////////////////////////////////////////////////////////////////

#include "envDep.h"

#include "GL/glew.h"
#if defined(WIN32)
#include "GL/wglew.h"
#endif

#include <GLFW/glfw3.h>

#include "mydef.h"

#include <Point3.h>
#include <Vector3.h>
#ifdef VM_INCLUDE_NAMESPACE
using namespace kh_vecmath;
#endif // VM_INCLUDE_NAMESPACE

#include "GLPanel.hxx"
#include "GLMaterial.hxx"
#include "VWIO.hxx"
#include "PNGImage.hxx"
#include "PNGImage.hxx"

int width = 800;
int height = 800;
// int width = 1024;
// int height = 1024;

static float myLight[] = {
  0.0f, 0.0f, 100.0f, 1.0f,
  0.8f, 0.8f, 0.8f, 1.0f,
  1.0f, 1.0f, 1.0f, 1.0f,
  0.8f, 0.8f, 0.8f, 1.0f,
  1.0f, // 0.1f,
  0.0f,  // 0.05f
  0.0f  // 0.05f
};

static float myMatl[] = {
  0.2f, 0.2f, 0.2f, 1.0f, 
  //  0.6f, 0.8f, 0.6f, 1.0f, 
  0.8f, 0.8f, 0.6f, 1.0f, 
  // 0.6f, 0.6f, 0.8f, 1.0f, 
  //  0.6f, 0.8f, 0.8f, 1.0f, 
  0.0f, 0.0f, 0.0f, 1.0f, 
  0.8f, 0.8f, 0.8f, 1.0f, 
  80.0f
};

GLPanel pane;
GLPanel pane2d;

// for picking
GLint viewport[4];
GLdouble MV[16];
GLdouble P[16];

bool pngflag = false;
// keyboard
bool shift_key_pressed = false;
bool control_key_pressed = false;
// mouse
bool left_button_pressed = false;
bool right_button_pressed = false;

#define RENDER_VBO 0
#define RENDER_VAR 1
unsigned short render_mode = RENDER_VBO;

// current shader for shading
unsigned short shading_program = PHONG_SHADING;

// #include <sstream>

// #include "strutil.h"
//#include "gzfileopen.h"

//GzFileOpen gzfile;

#if 0
/*
** シェーダ
*/
static GLuint vertShader;
static GLuint fragShader;
static GLuint gl2Program;
#endif

////////////////////////////////////////////////////////////////////////////////////

#include "nvtimer.h"

timer fps(10);
char buf[BUFSIZ];
char txt[BUFSIZ];
float max_fps = 0.0f;
bool first = true;

////////////////////////////////////////////////////////////////////////////////////

#include "SMFRIO.hxx"
#include "OFFRIO.hxx"
#include "MeshR.hxx"
#include "GLMeshVBO.hxx"

MeshR   meshR;
int tex_id;
MeshR   pmeshR;

GLMeshVBO glmeshvbo;
GLMeshVBO glmeshvbop;

bool isDrawParamMesh = false;

////////////////////////////////////////////////////////////////////////////////////

#include "MyMesh.hxx"
#include "Orbifold.hxx"
#include "Param.hxx"
#include "ToMeshR.hxx"

MyMesh  mesh;
Orbifold orbi;
std::vector<MyMesh::VertexHandle> cs_vertices;
bool isCalculated = false;

////////////////////////////////////////////////////////////////////////////////////

void checkGLErrors(char *s)
{
  GLenum error;
  while ((error = glGetError()) != GL_NO_ERROR) {
    fprintf(stderr, "%s: error - %s\n", s, (char *) gluErrorString(error));
  }
}

////////////////////////////////////////////////////////////////////////////////////

void display()
{
  if ( isDrawParamMesh == false )
    {
      pane.clear( width, height );
      pane.setView();

#if 1
      glMatrixMode(GL_PROJECTION);
      glGetIntegerv(GL_VIEWPORT, viewport); 
      glGetDoublev(GL_PROJECTION_MATRIX, P);
      glMatrixMode(GL_MODELVIEW);
      glGetDoublev(GL_MODELVIEW_MATRIX, MV);
#endif

      pane.setLight();

      pane.changeProgram( shading_program );

      if ( shading_program == PHONG_TEXTURE )
        {
          // ::glEnable( GL_TEXTURE_2D );
          ::glBindTexture( GL_TEXTURE_2D, meshR.texID() );
        }

      glmeshvbo.drawShading();

      if ( shading_program == PHONG_TEXTURE )
        {
          ::glBindTexture( GL_TEXTURE_2D, 0 );
          // ::glDisable( GL_TEXTURE_2D );
        }

      if ( glmeshvbo.isDrawWireframe() )
        {
          pane.changeProgram( WIREFRAME );
          // glmeshvbo.drawWireframe();

          // cone singularity points
          glDisable( GL_LIGHTING );
          glColor3f( 1.0f, .0f, .0f );
          glPointSize( 5.0f );
          glBegin(GL_POINTS);
          for ( int i = 0; i < cs_vertices.size(); ++i )
            {
              MyMesh::Point p = mesh.point( cs_vertices[i] );
              glVertex3d( p[0], p[1], p[2] );
            }
          glEnd();

          // boundary path
          glColor3f( .0f, .0f, .0f );
          glLineWidth( 3.0f );
          std::vector<std::vector<MyMesh::VertexHandle> >& path = orbi.path();
          for ( int i = 0; i < path.size() / 2; ++i ) // search for only two paths
            {
              // cone singularities are not concerned
              glBegin(GL_LINE_STRIP);
              for ( int j = 0; j < path[i].size(); ++j )
                {
                  MyMesh::Point p = mesh.point( path[i][j] );
                  glVertex3d( p[0], p[1], p[2] );
                }
              glEnd();
            }
        }
      pane.finish();
      //   checkGLErrors("display");
    }
  else // draw 2d
    {
      pane2d.clear( width, height );
      pane2d.setView();

      // pane2d.setLight();
      // pane2d.changeProgram( PHONG_SHADING );

      pane.changeProgram( WIREFRAME );
      glDisable( GL_LIGHTING );
 
      glmeshvbop.drawShading();

      if ( glmeshvbop.isDrawWireframe() )
        {
          pane.changeProgram( WIREFRAME );
          glDisable( GL_LIGHTING );

          std::vector<std::vector<MyMesh::VertexHandle> >& path = orbi.path();
          // boundary path
          glColor3f( .0f, .0f, 1.0f );
          glLineWidth( 5.0f );
          for ( int i = 0; i < path.size(); ++i ) // search for only two paths
            {
              // cone singularities are not concerned
              glBegin(GL_LINE_STRIP);
              for ( int j = 0; j < path[i].size(); ++j )
                {
                  MyMesh::TexCoord2D t = mesh.texcoord2D( path[i][j] );
                  glVertex3d( t[0], t[1], 0.0 );
                }
              glEnd();
            }
          // cone singularity points
          glColor3f( 1.0f, .0f, .0f );
          glPointSize( 5.0f );
          glBegin(GL_POINTS);
          for ( int i = 0; i < path.size(); ++i ) // search for only two paths
            {
              for ( int j = 0; j < path[i].size(); ++j )
                {
                  if ( (j == 0) || (j == path[i].size()-1) )
                    {
                      MyMesh::TexCoord2D t = mesh.texcoord2D( path[i][j] );
                      glVertex3d( t[0], t[1], 0.0 );
                    }
                }
            }
          glEnd();
        }

      pane2d.finish();
    }

}

////////////////////////////////////////////////////////////////////////////////////

static void error_callback(int error, const char* description)
{
  fputs(description, stderr);
}

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
  // ESC
  if ( (key == GLFW_KEY_ESCAPE) && (action == GLFW_PRESS) )
    glfwSetWindowShouldClose(window, GL_TRUE);

  // q
  else if ( (key == GLFW_KEY_Q) && (action == GLFW_PRESS) )
    glfwSetWindowShouldClose(window, GL_TRUE);

  // 3 PHONG_SHADING
  else if ( (key == GLFW_KEY_3) && (action == GLFW_PRESS) )
    {
      shading_program = PHONG_SHADING;
      pane.changeProgram( shading_program );
      isDrawParamMesh = false;
    }

  // 4 GOURAND_SHADING
  else if ( (key == GLFW_KEY_4) && (action == GLFW_PRESS) )
    {
      shading_program = GOURAND_SHADING;
      pane.changeProgram( shading_program );
      isDrawParamMesh = false;
    }

  // 5 PHONG_TEXTURE
  else if ( (key == GLFW_KEY_5) && (action == GLFW_PRESS) )
    {
      shading_program = PHONG_TEXTURE;
      pane.changeProgram( shading_program );
      isDrawParamMesh = false;
    }

  // 6 Parameter Mesh
  else if ( (key == GLFW_KEY_6) && (action == GLFW_PRESS) )
    {
      // pane.initViewParameters( width, height );
      // // pane.setViewPoint( .5f, .5f, -2.5f );
      // // pane.setLookPoint( .5f, .5f, .0f );
      // pane.initView();
      // shading_program = PHONG_SHADING;
      // pane.changeProgram( shading_program );
      isDrawParamMesh = true;
    }

  // w
  else if ( (key == GLFW_KEY_W) && (action == GLFW_PRESS) )
    {
      if ( glmeshvbo.isDrawWireframe() == false )
        {
          glmeshvbo.setIsDrawWireframe( true );
          glmeshvbop.setIsDrawWireframe( true );
        }
      else
        {
          glmeshvbo.setIsDrawWireframe( false );
          glmeshvbop.setIsDrawWireframe( false );
        }
    }
  
  // x
  else if ( (key == GLFW_KEY_X) && (action == GLFW_PRESS) )
    {
#if 1
      if ( cs_vertices.size() != 3 )
        {
          cout << "cannot compute orbifold parameterization. " << endl;
          return;
        }
#endif

      // // for sphere50
      // cs_vertices.push_back( mesh.vertex_handle(2) );
      // cs_vertices.push_back( mesh.vertex_handle(23) );
      // cs_vertices.push_back( mesh.vertex_handle(17) );
      
      // parameterization based on Euclidean orbifold
      std::vector<double> paramx;
      std::vector<double> paramy;
      orbi.setCSVertices( cs_vertices );
      Param param;
      // param.applyParam_Orbifold( mesh, orbi, paramx, paramy, BICGSTAB, MVW );
      // param.applyParam_Orbifold( mesh, orbi, paramx, paramy, SPARSELU, MVW );
      param.applyParam_Orbifold( mesh, orbi, paramx, paramy, SPARSELU, COTW );

#if 0
      // store param to vertex
      MyMesh::VertexIter v_it, v_end(mesh.vertices_end());
      for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
        {
          MyMesh::VertexHandle vh = *v_it;
          int id = vh.idx();
          MyMesh::Point p( paramx[id], paramy[id], 0.0 );
          mesh.set_point( vh, p );
        }
#endif

      // store param to texture 2d
      mesh.request_vertex_texcoords2D();
      MyMesh::VertexIter v_it, v_end(mesh.vertices_end());
      for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
        {
          MyMesh::VertexHandle vh = *v_it;
          int id = vh.idx();
          MyMesh::TexCoord2D t( paramx[id], paramy[id] );
          mesh.set_texcoord2D( vh, t );
        }

      // write mesh to output.obj
      std::cout << "save mesh ... ";
      OpenMesh::IO::Options opt;
      opt += OpenMesh::IO::Options::VertexTexCoord;
      // opt += OpenMesh::IO::Options::VertexColor;
      if ( !OpenMesh::IO::write_mesh(mesh, "output.obj", opt) )
        {
          std::cerr << "Cannot write mesh to file. " << std::endl;
          return;
        }
      std::cout << "done." << std::endl;

      // cs_vertices.clear();
      isCalculated = true;

      // meshR recreation
      meshR.clear();
      ToMeshR tomeshr;
      tomeshr.apply( mesh, meshR );
      meshR.createVertexNormals();
      meshR.setTexID( tex_id );

      // pmeshR creation
      tomeshr.apply_param( mesh, pmeshR );

      // GLMeshVBO recreation
      glmeshvbo.clear();
      glmeshvbo.setMesh( meshR );

      // GLMeshVBOp creation
      glmeshvbop.setMesh( pmeshR );
    }

  // i
  else if ( (key == GLFW_KEY_I) && (action == GLFW_PRESS) )
    {
      cout << "output to .vw file ... " << endl;
      VWIO vw_out;
      vw_out.outputToFile( "tmp.vw", pane.manip() );
      cout << "done." << endl;
    }

  // p
  else if ( (key == GLFW_KEY_P) && (action == GLFW_PRESS) )
    {
      pngflag = true;
    }

  // shift
  else if ( (key == GLFW_KEY_LEFT_SHIFT ) && (action == GLFW_PRESS) )
    {
      shift_key_pressed = true;
    }
  else if ( (key == GLFW_KEY_LEFT_SHIFT ) && (action == GLFW_RELEASE) )
    {
      shift_key_pressed = false;
    }

  // control
  else if ( (key == GLFW_KEY_LEFT_CONTROL ) && (action == GLFW_PRESS) )
    {
      control_key_pressed = true;
    }
  else if ( (key == GLFW_KEY_LEFT_CONTROL ) && (action == GLFW_RELEASE) )
    {
      control_key_pressed = false;
    }
}

static void mousebutton_callback(GLFWwindow* window, int button, int action, int mods)
{
  double xd, yd;
  glfwGetCursorPos( window, &xd, &yd );
  pane.setScreenXY( (int) xd, (int) yd );

  int x = (int) xd;
  int y = (int) yd;

  if ( (button == GLFW_MOUSE_BUTTON_1) && (action == GLFW_PRESS) )
    {
      // cout << "left button pressed." << endl;
      left_button_pressed = true;
      pane.startRotate();
      pane.startZoom();
      pane.startMove();

      if ( isCalculated == false )
        {
          int oldX = x; 
          int oldY = y; 
          int window_y = (height - y);
          float norm_y = float(window_y)/float(height/2.0);
          int window_x = x ;
          float norm_x = float(window_x)/float(width/2.0);

          float winZ=0;
          glReadPixels( x, height-y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );
          double objX=0, objY=0, objZ=0;
          gluUnProject(window_x,window_y, winZ,  MV,  P, viewport, &objX, &objY, &objZ);
          MyMesh::Point q( objX, objY, objZ );

          MyMesh::VertexIter v_it, v_end(mesh.vertices_end());
          for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
            {
              MyMesh::Point p = mesh.point( *v_it );

              if( (p-q).norm()<0.01)
                {
                  int sel_i = v_it->idx();
                  cout << "vt " << sel_i << " selected. " << endl;
                  cs_vertices.push_back( *v_it );
                  // printf("Intersected at %d\n",i);
                  // printf("Pt [ %3.3f,%3.3f,%3.3f ]\n",X[i].x(), X[i].y(), X[i].z());
                  break;
                }
            }
        }
    }
  else if ( (button == GLFW_MOUSE_BUTTON_1) && (action == GLFW_RELEASE) )
    {
      left_button_pressed = false;
      pane.finishRMZ();
    }
  else if ( (button == GLFW_MOUSE_BUTTON_2) && (action == GLFW_PRESS) )
    {
      // cout << "right button pressed." << endl;
      right_button_pressed = true;
    }
  else if ( (button == GLFW_MOUSE_BUTTON_2) && (action == GLFW_RELEASE) )
    {
      right_button_pressed = false;
    }
}

static void cursorpos_callback(GLFWwindow* window, double xd, double yd )
{
  int x = (int) xd;
  int y = (int) yd;

  if ( left_button_pressed && !shift_key_pressed && !control_key_pressed )
    {
      pane.updateRotate( x, y );
    }
  else if ( left_button_pressed && shift_key_pressed && !control_key_pressed )
    {
      pane.updateZoom( x, y );
    }
  else if ( left_button_pressed && !shift_key_pressed && control_key_pressed )
    {
      pane.updateMove( x, y );
    }
}

// mouse wheel
static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
  //cout << yoffset << endl;
  pane.updateWheelZoom( yoffset );
}

// window resize
static void windowsize_callback(GLFWwindow* window, int w, int h )
{
  width = w;
  height = h;
  pane.changeSize( w, h );
}

////////////////////////////////////////////////////////////////////////////////////

int main( int argc, char **argv )
{
  //char filename[BUFSIZ];
  if ( argc != 3 )
    {
      std::cerr << "Usage: " << argv[0] << " in.off|obj|ply in.png ." << std::endl;
      exit(1);
    }

  mesh.request_vertex_colors();
  OpenMesh::IO::Options options;
  options += OpenMesh::IO::Options::VertexColor;
  if ( !OpenMesh::IO::read_mesh(mesh, argv[1], options) )
    {
      std::cerr << "Error: Cannot read mesh from " << argv[1] << std::endl;
      return 1;
    }

  cout << "options.vertex_has_color() = " << options.vertex_has_color() << endl;
  cout << "mesh.has_vertex_colors() = " << mesh.has_vertex_colors() << endl;

  // MyMesh::VertexIter v_it, v_end(mesh.vertices_end());
  // for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
  //   {
  //     MyMesh::VertexHandle vh = *v_it;
  //     MyMesh::Color c = mesh.color( vh );
  //     printf("%u %u %u\n", c[0], c[1], c[2] );
  //   }

  orbi.init( mesh );

  ToMeshR tomeshr;
  tomeshr.apply( mesh, meshR );
  meshR.createVertexNormals();

  meshR.printInfo();

  // GLGW initialization
  glfwSetErrorCallback(error_callback);
  if ( !glfwInit() ) return EXIT_FAILURE;

  // glfwWindowHint(GLFW_SAMPLES, 4);
  // glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  // glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
  // glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  // glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  GLFWwindow* window = glfwCreateWindow( width, height, "Orbifold Parameterization", NULL, NULL );
  if ( !window )
    {
      glfwTerminate();
      return EXIT_FAILURE;
    }
  glfwMakeContextCurrent( window );
  glfwSetKeyCallback( window, key_callback );
  glfwSetMouseButtonCallback( window, mousebutton_callback );
  glfwSetCursorPosCallback( window, cursorpos_callback );
  glfwSetScrollCallback( window, scroll_callback );
  glfwSetWindowSizeCallback( window, windowsize_callback );

  pane.init( width, height );

  // Initialize some state for window's rendering context.
  pane.initGL();
  pane.initGLEW();

  std::vector<unsigned char> img;
  int format, w, h;
  tex_id = pane.loadTexture( argv[2], img, &format, &w, &h );
  pane.assignTexture( tex_id, img, format, w, h );

  if ( !pane.initShader() ) return -1;
  
  pane.setViewPoint( .0f, .0f, 2.5f );
  pane.setIsGradientBackground( false );
  //pane.setLightParameters( 0, myLight );

  // glfwSwapInterval(1);
  glfwSwapInterval(0);
#if defined(WIN32)
  //wglSwapIntervalEXT(0);
#endif

  glmeshvbo.setMesh( meshR );
  glmeshvbo.setIsSmoothShading( true );
  glmeshvbo.setIsDrawWireframe( true );

  // pane 2d
  pane2d.init( width, height );
  pane2d.initGL();
  pane2d.initGLEW();
  pane2d.setViewPoint( .5f, .5f, -2.5f );
  pane2d.setLookPoint( .5f, .5f, .0f );
  pane2d.setIsGradientBackground( false );

  glmeshvbop.setIsSmoothShading( true );
  glmeshvbop.setIsDrawWireframe( true );

  // GLFW rendering process
  while ( !glfwWindowShouldClose(window) )
    {
      display();

      // for measuring fps
      fps.frame();
      if ( fps.timing_updated() )
        {
          float f = fps.get_fps();
          if ( max_fps < f ) max_fps = f;
          sprintf( buf,"%.3f fps - max %.3f fps", f, max_fps );
        }
      sprintf( txt, "Orbifold Parameterization - %s", buf );
      glfwSetWindowTitle( window, txt );

      if ( pngflag )
        {
          PNGImage pi( width, height, false );
          pi.capture_and_write("screen.png");
          pngflag = false;
        }

      glfwSwapBuffers(window);
      glfwPollEvents();
    }

  glfwDestroyWindow(window);
  glfwTerminate();

  return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////////

