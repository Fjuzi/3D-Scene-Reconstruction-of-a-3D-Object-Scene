# 3D-Scene-Reconstruction-of-a-3D-Object-Scene
3D Scene reconstruction from multiple images and camera calibration
This repo contains a 3D reconstruction of a 3D object/scen. It includes the calibration of a camera, feature extraction and matching of feature points between the views, to computation the fundamental matrix between views to ob a 3D point cloud reconstruction. And finally it includes the representatin of the object geometric model over this point cloud.

![alt text](https://github.com/Fjuzi/fjuzi.github.io/blob/master/dist/assets/img/vmmc.gif)

The project is divided into multiple sections:

### Section 1: Obtention of the intrinsic parameters of a camera
In this section we calibrate the camera, which will be used later on. We get the intrinsic parameters via a pyhsical checkerboard.

### Section 2: Finding local matches between several views of an object.
In this section several views of a scene will be captured. Afterwards for pairs of views, detection, description and matching of features is presented. Several methods are investigated for detector and descriptor combinations:

- Difference of Hessian (DoH) + SIFT
- SURF + SURF
- Kaze + Kaze
- SIFT + DSP-SIFT

Finally the perfomrance in evaluated via fundamental matrix estimation

### Section 3:  3D reconstruction and calibration

From the intrinsic parameters of the camera from Section 1, and the feature point matches between images, the final 3D reconstruction is obtained.



The whole report can be read [at this pdf](https://github.com/Fjuzi/3D-Scene-Reconstruction-of-a-3D-Object-Scene/blob/master/3D%20Reconstruction%20of%20a%20Scene.pdf).
