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

int width = 512;
int height = 512;
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

#include <sstream>

#include "strutil.h"
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
#include "GLMeshR.hxx"

MeshR   meshR;

GLMeshVBO glmeshvbo;
GLMeshR glmeshr;

////////////////////////////////////////////////////////////////////////////////////

#include "MyMesh.hxx"
#include "ShortestPathDijkstra.hxx"
#include "ToMeshR.hxx"

MyMesh  mesh;
std::vector<MyMesh::VertexHandle> path;
ShortestPathDijkstra sp;

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
  pane.clear( width, height );
  pane.setView();

  float current_light[3];
  pane.getRealLightPosition( current_light );
  pane.setLightPos( 0, current_light[0], current_light[1], current_light[2], 1.0f );

  pane.setLight();

  pane.changeProgram( shading_program );
  if ( render_mode == RENDER_VBO )
    {
      glmeshvbo.drawShading();
    }
  else if ( render_mode == RENDER_VAR )
    {
      //       glmeshr.drawColor();
      glmeshr.drawShading();
    }

  if ( glmeshvbo.isDrawWireframe() )
    {
      pane.changeProgram( WIREFRAME );
      if ( render_mode == RENDER_VBO )
        {
          glmeshvbo.drawWireframe();
        }
      else if ( render_mode == RENDER_VAR )
        {
          //       glmeshr.drawColor();
          glmeshr.drawWireframe();
        }

      // path
      glDisable( GL_LIGHTING );
      glColor3f( .0f, .0f, .0f );
      glLineWidth( 3.0 );
      glBegin(GL_LINE_STRIP);
      for ( int i = 0; i < path.size(); ++i )
        {
          MyMesh::Point p = mesh.point( path[i] );
          glVertex3d( p[0], p[1], p[2] );
        }
      glEnd();
    }

  pane.finish();
  //   checkGLErrors("display");
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

  // 1
  else if ( (key == GLFW_KEY_1) && (action == GLFW_PRESS) )
    render_mode = RENDER_VBO;

  // 2
  else if ( (key == GLFW_KEY_2) && (action == GLFW_PRESS) )
    render_mode = RENDER_VAR;

  // 3 PHONG_SHADING
  else if ( (key == GLFW_KEY_3) && (action == GLFW_PRESS) )
    {
      shading_program = PHONG_SHADING;
      pane.changeProgram( shading_program );
    }

  // 4 GOURAND_SHADING
  else if ( (key == GLFW_KEY_4) && (action == GLFW_PRESS) )
    {
      shading_program = GOURAND_SHADING;
      pane.changeProgram( shading_program );
    }

  // new path
  else if ( (key == GLFW_KEY_N) && (action == GLFW_PRESS) )
    {
      cout << "new path: " << endl;
      int si, ei;
      cin >> si >> ei;

      path.clear();
      MyMesh::VertexHandle sv = mesh.vertex_handle(si);
      MyMesh::VertexHandle ev = mesh.vertex_handle(ei);
      sp.apply(sv, ev, path);
    }

  // w
  else if ( (key == GLFW_KEY_W) && (action == GLFW_PRESS) )
    {
      if ( glmeshvbo.isDrawWireframe() == false )
        {
          glmeshvbo.setIsDrawWireframe( true );
        }
      else
        {
          glmeshvbo.setIsDrawWireframe( false );
        }
    }
  
  // s
  else if ( (key == GLFW_KEY_S) && (action == GLFW_PRESS) )
    {
      cout << "output to file ... " << endl;
      SMFRIO rio_output( meshR );
      rio_output.outputToFile( "tmp.smf" );
      cout << "done." << endl;
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

  if ( (button == GLFW_MOUSE_BUTTON_1) && (action == GLFW_PRESS) )
    {
      // cout << "left button pressed." << endl;
      left_button_pressed = true;
      pane.startRotate();
      pane.startZoom();
      pane.startMove();
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
  // char filename[BUFSIZ];
  if ( (argc < 2) || (argc > 2) )
    {
      std::cerr << "Usage: " << argv[0] << " in.off|ply|obj ." << std::endl;
      exit(1);
    }

  if ( ! OpenMesh::IO::read_mesh(mesh, argv[1]) )
    {
      std::cerr << "Error: Cannot read mesh from " << argv[1] << std::endl;
      return 1;
    }

  sp.init( mesh );
  MyMesh::VertexHandle sv = mesh.vertex_handle(0);
  MyMesh::VertexHandle ev = mesh.vertex_handle(1000);
  sp.apply(sv, ev, path);

  for ( int i = 0; i < path.size(); ++i)
    {
      cout << path[i].idx() << endl;
    }

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

  GLFWwindow* window = glfwCreateWindow( width, height, "VBO GLFW", NULL, NULL );
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

  if ( !pane.initShader() ) return -1;
  
  // Point3f p( .0f, .0f, 25.0f );
  // Vector3f v( .0f, .0f, -25.0f );
  Point3f p( .0f, .0f, 2.5f );
  Vector3f v( .0f, .0f, -2.5f );
  pane.setViewPoint( p );
  pane.setViewVector( v ); 
  pane.setIsGradientBackground( false );
  //pane.setLightParameters( 0, myLight );

  // glfwSwapInterval(1);
  glfwSwapInterval(0);
#if defined(WIN32)
  //wglSwapIntervalEXT(0);
#endif

  VWIO vw_in;
  if ( argc == 3 )
    {
      vw_in.inputFromFile( argv[2], pane.manip() );
    }

  glmeshr.setMesh( meshR );
  glmeshr.setIsSmoothShading( true );
  // glmeshr.setIsSmoothShading( false );
  // glmeshr.setIsDrawWireframe( true );
  //glmeshr.setMaterial( myMatl );
  //glmeshr.setMaterial( 0 );

  glmeshvbo.setMesh( meshR );
  glmeshvbo.setIsSmoothShading( true );
  //glmeshvbo.setMaterial( myMatl );
  //glmeshvbo.setMaterial( 0 );

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
      sprintf( txt, "VBO GLFW - %s", buf );
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

#if 0  
  mesh.add_property(elen);
  MyMesh::EdgeIter e_it, e_end(mesh.edges_end());
  for (e_it=mesh.edges_begin(); e_it!=e_end; ++e_it)
    {
      MyMesh::HalfedgeHandle h = mesh.halfedge_handle( *e_it, 0 );
      MyMesh::VertexHandle v0 = mesh.to_vertex_handle( h );
      MyMesh::VertexHandle v1 = mesh.from_vertex_handle( h );
      MyMesh::Point p0 = mesh.point( v0 );
      MyMesh::Point p1 = mesh.point( v1 );
      Elen& el = mesh.property( elen, *e_it );
      el.setLen( (p1 - p0).length() );
    }
 
  //
  MyMesh::VertexIter v_it, v_end(mesh.vertices_end());
  for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
    {
      MyMesh::VertexHandle vh = *v_it;
      MyMesh::VertexEdgeIter ve_it;
      double err = .0;
      for ( ve_it = mesh.ve_iter( vh ); ve_it.is_valid(); ++ve_it )
        {
          MyMesh::EdgeHandle eh = *ve_it;
          Elen& el = mesh.property( elen, eh );
          err += el.len();
        }
      pq.push( err, vh );
    }

  while ( !(pq.empty()) )
    {
      cout << "err " << pq.top_err() << " v " << pq.top_vh().idx() << endl;
      pq.pop();
    }
#endif
}

