# orbiparam

orbiparam is a sample implementation of the following paper:

Noam Aigerman, Yaron Lipman: ``Orbifold Tutte embeddings'', ACM Transactions on Graphics, 34(6), Article No.190, pp.1-12, 2015.

Only type I is implemented in this implementation.

## Getting Started

At first, you may try using the binary release (x64),
which is available [here](https://github.com/kanait/orbiparam/releases/tag/v1.0).
Uncompress the zip file and then execute orbiparam_bimba100K_x64_ao.bat.

On the display window, click the left mouse button three times to select three points on the surface of a mesh, and press the 'x' key to compute parameterization. If the process successfully finishes, a shortest path (black line) passing through the three points will be displayed.

Press the '5' key to show the texture mapping of the mesh using the computed 2D texture parameters.

Press the '6' key to show the 2D parameterization.

Press the 'w' key to hide/appear the shortest path.

Press the 's' key to save an OBJ file with the computed 2D texture parameters.

## Compilation

The following libraries are required for successfully compiling this software.

### [OpenMesh](https://www.openmesh.org) (Version 8.1 or higher)
### [GLFW](https://www.glfw.org/)
### [GLEW](http://glew.sourceforge.net/)
### [Eigen](https://gitlab.com/libeigen/eigen)
### [stb](https://github.com/nothings/stb)
### [render_Eigen](https://github.com/kanait/render_Eigen)

First, execute the "git clone" command with the "--recursive" option. This will clone the repository and automatically install Eigen, stb, and render_Eigen as submodules.

```
% git clone https://github.com/kanait/orbiparam.git --recursive
```
### ubuntu

OpenMesh is installed from source.
```
% tar zxf OpenMesh-x.x.tar.gz
% cd OpenMesh-x.x
% mkdir build
% cd build
% cmake .. && sudo make all install
```

GLFW and GLEW are installed by apt
```
% sudo apt -y install libglfw3-dev libglew-dev
```

Then execute in orbiparam directory as follows:

```
% cd orbiparam
% mkdir build
% cd build
% cmake ..
% make
```
Once the compilation process is completed successfully, an executable named "orbiparam" will be created. 

### Windows

To compile the software on Visual Studio 2019, a solution file named "orbiparam-vs2019.sln" is provided. For Visual Studio 2022, there is also a solution file named "orbiparam-vs2022.sln".

To successfully compile the software, you will need the following libraries for Windows: OpenMesh, GLFW, and GLEW.
Make sure you have these libraries installed on your system. Then, set the include and library directories appropriately in the project files of Visual Studio 2019 and Visual Studio 2022. This can be done in the project's properties, where you can specify the paths for including headers and linking libraries.

## Authors

* **[Takashi Kanai](https://graphics.c.u-tokyo.ac.jp/hp/en/)** - The University of Tokyo

## License

This software is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
