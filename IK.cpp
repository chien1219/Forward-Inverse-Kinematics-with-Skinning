#include "IK.h"
#include "FK.h"
#include "minivectorTemplate.h"
#include <Eigen/Dense>
#include <adolc/adolc.h>
#include <cassert>
#include <vector>
#if defined(_WIN32) || defined(WIN32)
  #ifndef _USE_MATH_DEFINES
    #define _USE_MATH_DEFINES
  #endif
#endif
#include <math.h>
using namespace std;

// CSCI 520 Computer Animation and Simulation
// Jernej Barbic and Yijing Li

namespace
{

// Converts degrees to radians.
template<typename real>
inline real deg2rad(real deg) { return deg * M_PI / 180.0; }

template<typename real>
Mat3<real> Euler2Rotation(const real angle[3], RotateOrder order)
{
  Mat3<real> RX = Mat3<real>::getElementRotationMatrix(0, deg2rad(angle[0]));
  Mat3<real> RY = Mat3<real>::getElementRotationMatrix(1, deg2rad(angle[1]));
  Mat3<real> RZ = Mat3<real>::getElementRotationMatrix(2, deg2rad(angle[2]));

  switch(order)
  {
    case RotateOrder::XYZ:
      return RZ * RY * RX;
    case RotateOrder::YZX:
      return RX * RZ * RY;
    case RotateOrder::ZXY:
      return RY * RX * RZ;
    case RotateOrder::XZY:
      return RY * RZ * RX;
    case RotateOrder::YXZ:
      return RZ * RX * RY;
    case RotateOrder::ZYX:
      return RX * RY * RZ;
  }
  assert(0);
}

// Performs forward kinematics, using the provided "fk" class.
// This is the function whose Jacobian matrix will be computed using adolc.
// numIKJoints and IKJointIDs specify which joints serve as handles for IK:
//   IKJointIDs is an array of integers of length "numIKJoints"
// Input: numIKJoints, IKJointIDs, fk, eulerAngles (of all joints)
// Output: handlePositions (world-coordinate positions of all the IK joints; length is 3 * numIKJoints)
template<typename real>
void forwardKinematicsFunction(
    int numIKJoints, const int * IKJointIDs, const FK & fk,
    const std::vector<real> & eulerAngles, std::vector<real> & handlePositions)
{
  // Students should implement this.
  // The implementation of this function is very similar to function computeLocalAndGlobalTransforms in the FK class.
  // The recommended approach is to first implement FK::computeLocalAndGlobalTransforms.
  // Then, implement the same algorithm into this function. To do so,
  // you can use fk.getJointUpdateOrder(), fk.getJointRestTranslation(), and fk.getJointRotateOrder() functions.
  // Also useful is the multiplyAffineTransform4ds function in minivectorTemplate.h .
  // It would be in principle possible to unify this "forwardKinematicsFunction" and FK::computeLocalAndGlobalTransforms(),
  // so that code is only written once. We considered this; but it is actually not easily doable.
  // If you find a good approach, feel free to document it in the README file, for extra credit.
	int N = fk.getNumJoints();
	vector<Mat3<real>> localTransform(N);
	vector<Mat3<real>> globalTransform(N);
	vector<Vec3<real>> localTranslation(N);
	vector<Vec3<real>> globalTranslation(N);

	for (int i = 0; i < N; i++)
	{
		real eulerAngleArr[3], jointOrientationArr[3];
		for (int j = 0; j < 3; j++)
		{
			eulerAngleArr[j] = eulerAngles[3 * i + j];
			jointOrientationArr[j] = fk.getJointOrient(i)[j];
			localTranslation[i][j] = fk.getJointRestTranslation(i)[j];
		}
		Mat3<real> eulerAngleMat = Euler2Rotation(eulerAngleArr, fk.getJointRotateOrder(i));
		Mat3<real> jointOrientationMat = Euler2Rotation(jointOrientationArr, XYZ);
		
		localTransform[i] = jointOrientationMat * eulerAngleMat;
	}
	for (int i = 0; i < N; i++)
	{
		int j = fk.getJointUpdateOrder(i);
		int parent = fk.getJointParent(j);
		if (parent == -1)
		{
			globalTransform[j] = localTransform[i];
			globalTranslation[j] = localTranslation[j];
		}
		else
		{
			multiplyAffineTransform4ds(globalTransform[parent], globalTranslation[parent], localTransform[j], localTranslation[j], globalTransform[j], globalTranslation[j]);
		}

		for (int jID = 0; jID < numIKJoints; jID++)
		{
			if (j == IKJointIDs[jID])
			{
				handlePositions[jID * 3] = globalTranslation[j][0];
				handlePositions[jID * 3 + 1] = globalTranslation[j][1];
				handlePositions[jID * 3 + 2] = globalTranslation[j][2];
			}
		}
	}
}

} // end anonymous namespaces

IK::IK(int numIKJoints, const int * IKJointIDs, FK * inputFK, int adolc_tagID)
{
  this->numIKJoints = numIKJoints;
  this->IKJointIDs = IKJointIDs;
  this->fk = inputFK;
  this->adolc_tagID = adolc_tagID;

  FKInputDim = fk->getNumJoints() * 3;
  FKOutputDim = numIKJoints * 3;

  train_adolc();
}

void IK::train_adolc()
{
  // Students should implement this.
  // Here, you should setup adol_c:
  //   Define adol_c inputs and outputs. 
  //   Use the "forwardKinematicsFunction" as the function that will be computed by adol_c.
  //   This will later make it possible for you to compute the gradient of this function in IK::doIK
  //   (in other words, compute the "Jacobian matrix" J).
  // See ADOLCExample.cpp .

	int n = FKInputDim;
	int m = FKOutputDim;

	trace_on(adolc_tagID);

	vector<adouble> x(n);
	vector<adouble> y(m);

	for (int i = 0; i < n; i++)
		x[i] <<= 0;

	forwardKinematicsFunction(numIKJoints, IKJointIDs, (*fk), x, y);

	// output is handle global positions
	vector<double> output(m);
	for (int i = 0; i < m; i++)
		y[i] >>= output[i];

	// ADOL-C tracking finished
	trace_off();
}

void IK::doIK(const Vec3d * targetHandlePositions, Vec3d * jointEulerAngles)
{
	// Choose IK Method
	TikhonovIK(targetHandlePositions, jointEulerAngles);
	//PseudoInverse(targetHandlePositions, jointEulerAngles);
}


void IK::TikhonovIK(const Vec3d * targetHandlePositions, Vec3d * jointEulerAngles)
{
	// You may find the following helpful:
	int numJoints = fk->getNumJoints(); // Note that is NOT the same as numIKJoints!

	// Students should implement this.
	// Use adolc to evalute the forwardKinematicsFunction and its gradient (Jacobian). It was trained in train_adolc().
	// Specifically, use ::function, and ::jacobian .
	// See ADOLCExample.cpp .
	//
	// Use it implement the Tikhonov IK method (or the pseudoinverse method for extra credit).
	// Note that at entry, "jointEulerAngles" contains the input Euler angles. 
	// Upon exit, jointEulerAngles should contain the new Euler angles.

	// input and output dimension
	int n = FKInputDim;
	int m = FKOutputDim;

	// define the alpha value, the larger the alpha, the more stable the method
	double alpha = 0.01;

	// input values are current Euler angles
	// output values are current handle positions
	double* input_x_values = new double[n];
	double* output_y_values = new double[m];

	for (int i = 0; i < numJoints; i++)
		for (int j = 0; j < 3; j++)
			input_x_values[i * 3 + j] = jointEulerAngles[i][j];


	::function(adolc_tagID, m, n, input_x_values, output_y_values);

	// You can call ::jacobian(adolc_tagID, ...) as many times as you like to ask ADOL-C to evalute the jacobian matrix of f on different x:
	double* jacobianMatrix = new double[m * n]; // We store the matrix in row-major order.
	double** jacobianMatrixEachRow = new double*[m];
	for (int i = 0; i < m; i++)
		jacobianMatrixEachRow[i] = &jacobianMatrix[i * n]; // pointer array where each pointer points to one row of the jacobian matrix

	::jacobian(adolc_tagID, m, n, input_x_values, jacobianMatrixEachRow); // each row is the gradient of one output component of the function

	// solve (Jt*J + alpha*I)*theta = Jt * b using Eigen, note as Ax = B
	// b is delta pos, theta (x) is n*1 vector representing the change of Euler angles we look for.
	Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);
	Eigen::MatrixXd J(m, n), Jt(n, m);
	Eigen::MatrixXd A(n, n);
	Eigen::VectorXd B(m);

	for (int r = 0; r < m; r++)
		for (int c = 0; c < n; c++)
			J(r, c) = jacobianMatrix[r * n + c];

	Jt = J.transpose();
	A = Jt * J + alpha * I;

	Eigen::VectorXd deltaPos(m);
	for (int i = 0; i < numIKJoints; i++)
		for (int j = 0; j < 3; j++)
			deltaPos[i * 3 + j] = targetHandlePositions[i][j] - output_y_values[i * 3 + j];

	B = Jt * deltaPos;

	Eigen::VectorXd x = A.ldlt().solve(B);

	vector<Vec3d> deltaEuler(numJoints);
	for (int i = 0; i < numJoints; i++)
		deltaEuler[i] = { x[i * 3], x[i * 3 + 1], x[i * 3 + 2] };

	// Upodate final euler angles
	for (int i = 0; i < numJoints; i++)
		jointEulerAngles[i] += deltaEuler[i];
}


void IK::PseudoInverse(const Vec3d * targetHandlePositions, Vec3d * jointEulerAngles)
{
	int numJoints = fk->getNumJoints();

	// input and output dimension
	int n = FKInputDim;
	int m = FKOutputDim;

	// define the alpha value, the larger the alpha, the more stable the method
	double alpha = 0.01;

	// sub-step to solve the IK problem on a portion of the original distance if greater than maxDistance.
	double timeStep = 1;
	double maxDistance = 0.04;

	// input values are current Euler angles
	// output values are current handle positions
	double* input_x_values = new double[n];
	double* output_y_values = new double[m];

	for (int i = 0; i < numJoints; i++)
		for (int j = 0; j < 3; j++)
			input_x_values[i * 3 + j] = jointEulerAngles[i][j];


	::function(adolc_tagID, m, n, input_x_values, output_y_values);

	// You can call ::jacobian(adolc_tagID, ...) as many times as you like to ask ADOL-C to evalute the jacobian matrix of f on different x:
	double* jacobianMatrix = new double[m * n]; // We store the matrix in row-major order.
	double** jacobianMatrixEachRow = new double*[m];
	for (int i = 0; i < m; i++)
		jacobianMatrixEachRow[i] = &jacobianMatrix[i * n]; // pointer array where each pointer points to one row of the jacobian matrix

	::jacobian(adolc_tagID, m, n, input_x_values, jacobianMatrixEachRow); // each row is the gradient of one output component of the function

	Eigen::MatrixXd J(m, n), Jt(n, m);

	for (int r = 0; r < m; r++)
		for (int c = 0; c < n; c++)
			J(r, c) = jacobianMatrix[r * n + c];

	Jt = J.transpose();

	Eigen::VectorXd deltaPos(m);
	for (int i = 0; i < numIKJoints; i++)
	{
		for (int j = 0; j < 3; j++)
			deltaPos[i * 3 + j] = targetHandlePositions[i][j] - output_y_values[i * 3 + j];

		if (sqrt(deltaPos[i * 3] * deltaPos[i * 3] + deltaPos[i * 3 + 1] * deltaPos[i * 3 + 1] + deltaPos[i * 3 + 2] * deltaPos[i * 3 + 2]) > maxDistance)
			timeStep = 0.2;
	}
	
	Eigen::VectorXd deltaEuler = Jt * (J * Jt).inverse() * deltaPos;

	// Upodate final euler angles
	for (int i = 0; i < numJoints; i++)
		for (int j = 0; j < 3; j++)
			jointEulerAngles[i][j] += deltaEuler[i * 3 + j] * timeStep;
}