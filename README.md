# DSFPH

- [Introduction](#Introduction)
- [Compile](#Compilation)
- [Explanation](#Explain)


<a name="Introduction"></a>
## Introduction


This repository contains the code and the dependencies associated with a computer graphics project. It implements Divergence Free Smooth Particle Hydrodynamics (DSFPH), an advanced fuild simulation technique.
The code creates a 3D scene with 500 000 water particles in it and simulates it in a rocky environnment. There are also somme camera and birds animations in the scene.


<br>

<a name="Compilation"></a>
## Compile the project


This project is designed for Linux distros. In order to compile it, you need to install the development versions of glfw3 and openGL, and make.
To compile it, just place yourself at the root of the repo and :

```bash
make
```

You can then execute it :
```bash
./pgm
```


<br>

<a name="Explain"></a>
## Explanations


DSFPH comes from: https://animation.rwth-aachen.de/media/papers/2015-SCA-DFSPH.pdf
The scene code can be found under `scenes/DFSPH` folder.
All libraries except openGL and GLFW3 are contained within the `third_party` folder.