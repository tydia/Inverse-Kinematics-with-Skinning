# Inverse Kinematics with Skinning
### Implemented on macOS 10.15.4, tested on both macOS and Ubuntu 20.04
### [link to project description page](http://barbic.usc.edu/cs520-s20/assign3/) 

## Prerequisite: 
1. [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) library. This is mainly a header only library, so you don't need to install it.
2. [ADOL-C](https://github.com/coin-or/ADOL-C) library. Please visit the link for installation instructions.

## What is this project? 
The code as a whole is a complete character animation pipeline. Our job in this project is to implement the missing main functionalities to complete the pipeline. Specifically, for character animation, I implemented forward kinematics (FK), inverse kinematics (IK), and skinning. For IK, two methods are implemented for comparison: damped least squares method and the pseudo-inverse method. Damped least squares produces more stable simulation results, whereas pseudo-inverse method is less stable. However, pseudo-inverse method is more mathematically pleasing compared to damped least squares method. For skinning, Linear Blend Skinning and Dual Quaternion Skinning methods are implemented.


## Some Demos
1. Simulation with Damped Least Squares method:  
  [![](http://img.youtube.com/vi/Wtzs6BYJqIo/0.jpg)](http://www.youtube.com/watch?v=Wtzs6BYJqIo "IK with Damped Least Squares")
2. Simulation with Pseudo Inverse method:  
   [![](http://img.youtube.com/vi/bKWI_KLOr10/0.jpg)](http://www.youtube.com/watch?v=bKWI_KLOr10 "IK with Pseudo Inverse")
