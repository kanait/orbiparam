# orbiparam

orbiparam is a sample implementation of the following paper:

Noam Aigerman, Yaron Lipman: ``Orbifold Tutte embeddings'', ACM Transactions on Graphics, 34(6), Article No.190, pp.1-12, 2015.

Only type I is implemented in this implementation.

## Getting Started

At first, you may try to use binary release (x64), 
which is available from [here](https://github.com/kanait/orbiparam/releases/tag/v1.0).
Uncompress zip file and then execute orbiparam_bimba100K_x64_ao.bat.

On display window, click left mouse button three times to select three points on the surface of a mesh, and press 'x' key to compute parameterization. If the process is successfully finished, a shortest path (black line) passing three points is displayed.

press '5' key to show texture mapping of a mesh by using computed 2D texture parameters.

press '6' key to show 2D parameterization.

press 'w' key to hide/appear a shortest path.

press 's' key to save a obj file with computed 2D texture parameters.

## Compilation

The following libraries are required for successfully compiling this software.

### [OpenMesh](https://www.openmesh.org) (Version 8.1 or higher)
### [GLFW](https://www.glfw.org/)
### [GLEW](http://glew.sourceforge.net/)
### [Eigen](https://gitlab.com/libeigen/eigen)
### [stb](https://github.com/nothings/stb)
### [render_Eigen](https://github.com/kanait/render_Eigen)

First you execute "git clone" with with "--recursive" option. Then Eigen, stb, render_Eigen are installed automatically as submodules.

```
% git clone https://github.com/kanait/orbiparam.git --recursive
```
### ubuntu

OpenMesh is installed from source.
```
% tar zxf OpenMesh-8.1.tar.gz
% cd OpenMesh-8.1
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
Then an executable "orbiparam" is created if compilation completed successfully.

### Windows

A solution file (orbiparam-vs2019.sln) is provided for compiling on VS2019.
For compilation, windows library for OpenMesh, GLFW and GLEW are required and then set their include and library directories appropriately on the property of a VS2019 project file.

## Authors

* **[Takashi Kanai](https://graphics.c.u-tokyo.ac.jp/hp/en/)** - The University of Tokyo

## License

This software is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
