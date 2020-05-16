2020 Spring CSCI 520, Assignment 3

Name: Richard Chien

<Coding Environment>

Visual Studio 2017  | x64 | Debug

<Description of what you have accomplished>

The default skinning algorithm is Linear Blending Skinning (LBS).

================================== FINISHED PARTS ======================================

* Skinning (skinning.cpp)

1) Finished the applySkinning() function by using equation:
newVertexPosition = TotalSum(jointWeight * jointTransformMatrix * restVertexPosition)

* Forward Kinematics (FK.cpp)

1) Finished the computeLocalAndGlobalTransforms() function
To implement it:
Firstly, get the local transform matrices. Use this equation:
M(Local) = M(JointOritation) * M(EulertoRotation) *  M(Translations)
Secondly, go recursively to get the globalTransform. Use this equation:
globalTransform = parentGlobalTransform * localTransform

2) Finished the computeSkinningTransforms() function
To implement it, use this equation:
skinTransform = globalTransform * invRestTransform

* Inverse Kinematics (IK.cpp)

1) Finished the forwardKinematicsFunction() function
Use the same logic as computeLocalAndGlobalTransforms()

2) Finished the train_adolc() function
To implement it:
Firstly, define adol_c inputs.
Secondly, use forwardKinematicsFunction() as the computed function f to get the outputs.

3) Finished the doIK() function
To implement it:
Solve A x = b using Eigen
here A = (Jt*J + alpha*I)
b = Jt * (change of handle global positions: delta pos)
theta (x) is n*1 vector representing the change of Euler angles we want to find. 

4) Animation
In 300 pictures, I drag IK handle and record the result in Animation folder.


================================== EXTRA CREDITS 1========================================

Dual Quaternion Skinning (DQS)

To test this part, uncomment the DQS function in applySkinning() function in skinning.cpp, and comment LBS function.

To implement this:
1) Convert the jointSkinTransforms affine transformation matrixes into rotation matrices and translation vectors. Since I found that if I used Quaternion.h of Eigen, it could only convert 3*3 rotation matrix to quaternion, instead of 4*4 transformation matrix.
2) Use Eigen::Quaternion() convert rotation matrices to quaternions 
3) Compute normalized quaternions q using the related joints' weights. Compute translation vectors using the related joints' weights. 
4) Transform quaternion q back into rotation matrices. Combine rotation matrices with new translation vectors, using affine transform. 
5) Calculate the deformed skin vertex position, using the affine transformation matrices to multiply the rest position vectors.
 
Comparison between linear blend skinning and dual quaternion skinning:
1) DQS can avoid the loss of volume problem
2) Suppose we know two vertices positions: P1 and P2. If we use linear interpolation to get the new position, the new position will be lying on the segment between P1 and P2. If we use DQS, the new position will be lying on the arc circle contains P1 and P2, which will avoid mesh shrinkage.


================================== EXTRA CREDITS 2========================================

PseudoInverse IK

To test this part, uncomment the PseudoInverse() function in doIK() function in IK.cpp, and comment TikhonovIK() function.

To implement this:
1) Primitive steps are same as TikhonovIK(), but instead of computing Ax = B, I compute Jdagger by transpose(J) * inverse(J * transpose(J)) * deltaPos.


================================== EXTRA CREDITS 3========================================

Sub-Step IK process
When the user moves the IK handle for a long distance, divide the IK process into several sub-steps to improve the solution, where each sub-step solves the IK problem on a portion of the original distance.

This part is implemented in PseudoInverse() function with timeStep = 0.2 and maxDistance = 0.04;

To implement this:
1) If the distance of IK handle to its original pos is greater than maxDistance, then I sub-step the IK by times new euler result by timeStep.