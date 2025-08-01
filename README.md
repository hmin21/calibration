# IMPORTANT 
Please cite [1]	M. Han, F. Lei, W. Shi, S. Lu, and X. Li, “Uniaxial MEMS-based 3D reconstruction using pixel refinement,” Opt. Express, vol. 31, no. 1, pp. 536-554, 2023.

The folder of "Code" provides the calibration with pixel refiment and several reconstruction examples.

The folder of "Data" provides the parameters of camera calibration and the absolute phase of several reconstruction examples.

Please place the two folders (i.e., "Code" and "Data") in same root to order not to change the file path in the codes.

# Help documentation of file "Calibration.m"
The file of "Calibration.m"  can be used for the system calibration with the pixel refinement.

##Data preparation

1.cameraParams.mat:  This file is generated by MATLAB camera calibration tool kit, and is generated by software analysis through passing in the calibration plate image collected by the camera.

2.\Images\Phase:  This folder should store their phase maps in the order of the calibration board images, using the numbers 1,2,3..

3.You can change the number of placement positions of the calibration board N, the number of corners in the calibration board M, and the resolution(img_width x img_ height) of the picture according to the actual selection 

##Coordinate data pre acquisition

The first section（Line1-58） of this code file solves the coordinate value of the intersection of the ray of each pixel point with the plane where the calibration plate is located, according to the idea of pixel-by-pixel calibration, which is saved in the ‘camera_coordinate’ variable;
Line6-17：Load camera internal parameters and RT matrix of each calibration board;
Line19: Load the world coordinates of corners;
Line21-25: Load the phase map corresponding to each calibration board；
Line26-36: Estimating the equation of each calibration plate plane in camera coordinate system by least square method;
Line37-49: Solve the equation of ray corresponding to each pixel;
Line50-58: Coordinate of intersection point obtained by solving equation of joint ray and plane;

##Pixel-wise calibration process

Line59-80: This part uses the least square method to calibrate the parameters to be obtained, and finally stores them in 'Calibration_Table.mat'.

##Pixel-wise calibration error analysis

Line82-100: Evaluate the calibration errors of different pixels and store them in the relevant variables；

##Error visualization

Line102-115: Draw the frequency distribution histogram of errors, and grade the pixels with different errors. Finally, save the grades in 'Label_Table.mat'.


# Help documentation of file "Pixel_Refinement_Reconstruction.m"
The file can be used to reconstruct standard sphere with pixel with different grades, which are generated from the calibration with pixel refinement.
Users can perform it directly.

# Help documentation of file "Reconstruction.m"
The file can be used to reconstruct five complex surfaces with absolute phase.
Users can execute it directly.

# Help documentation of file "multi_frequency.m"
The file of "multi_frequency.m"  can be used to obtain the absolute phase with phase shifting plus multi-frequency algorithm. The absolute phase has been provided in the folder of "Data".

##Data preparation

1. The naming rule for images at the kth position, the ith frequency and the jth phase shift is k_ i_ j. bmp, which are stored in the /Image/ folder

##Details

You can select different frequency number m and phase shift number n, and you can also select the position where you want to calculate the phase as required. The final generated phase map will be saved in the form of a. mat file, and the distribution map on a line will be drawn at the end.
