#include "IK.h"
#include "FK.h"
#include "minivectorTemplate.h"
#include <Eigen/Dense>
#include <adolc/adolc.h>
#include <cassert>
#if defined(_WIN32) || defined(WIN32)
  #ifndef _USE_MATH_DEFINES
    #define _USE_MATH_DEFINES
  #endif
#endif
#include <math.h>
using namespace std;
using namespace Eigen;

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

  int numJoints = fk.getNumJoints();
 
  Mat3<adouble> currEulerMat3, currJtMat3,                                 // store temp conversion result
                localRotations[numJoints], globalRotations[numJoints];     // store rotations
  Vec3<adouble> localTranslation[numJoints], globalTranslation[numJoints]; // store translations

  // calculate local rotation and translation
  for(int i=0; i<numJoints; i++)
  {
    // get joint orientation and euler angle as adouble for adol-c to track
    // as input
    adouble adolcJointAngle[3] = {fk.getJointOrient(i)[0], 
                                  fk.getJointOrient(i)[1], 
                                  fk.getJointOrient(i)[2]};
    adouble adolcEulerAngle[3] = {eulerAngles[3 * i + 0],
                                  eulerAngles[3 * i + 1],
                                  eulerAngles[3 * i + 2]};

    // calculate local 3x3 rotation matrix and store for later use
    currJtMat3 = Euler2Rotation(adolcJointAngle, XYZ);
    currEulerMat3 = Euler2Rotation(adolcEulerAngle, fk.getJointRotateOrder(i));
    localRotations[i] = currJtMat3 * currEulerMat3;

    // get local translation 3-vector for later use
    localTranslation[i] = {fk.getJointRestTranslation(i)[0],
                           fk.getJointRestTranslation(i)[1],
                           fk.getJointRestTranslation(i)[2]};
  }
  
  // calculate handle positions
  for(int i=0; i<numJoints; i++)
  {
    // get current joint
    int currJoint = fk.getJointUpdateOrder(i);

    // recursively compute global transformation for calculating handle positions
    // here we want rotations and translations to be seprated for updating IK
    // joint positions, so we will use multiplyAffineTransform4ds to do so, which
    // is equivalent to multiplying two 4x4 transformation matrices
    if (fk.getJointParent(currJoint) == -1) {
      // M(g)root = M(l)root b/c root has no parent
      globalRotations[currJoint] = localRotations[currJoint];
      globalTranslation[currJoint] = localTranslation[currJoint];
    } else {
      // M(g)currJoint = M(g)currJoint'sParent * M(l)currJoint
      multiplyAffineTransform4ds(globalRotations[fk.getJointParent(currJoint)],   // R1
                                 globalTranslation[fk.getJointParent(currJoint)], // t1
                                 localRotations[currJoint],                       // R2
                                 localTranslation[currJoint],                     // t2
                                 globalRotations[currJoint],                      // Rout
                                 globalTranslation[currJoint]);                   // tout
    }
  }

  // set IK joint positions in world coord
  for(int i=0; i<numIKJoints; i++)
  {
      handlePositions[3 * i + 0] = globalTranslation[IKJointIDs[i]][0];
      handlePositions[3 * i + 1] = globalTranslation[IKJointIDs[i]][1];
      handlePositions[3 * i + 2] = globalTranslation[IKJointIDs[i]][2];
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

  // set input and output dimension, which were calculated in IK::IK()
  int n = FKInputDim;  // input dimension = # euler angles = numJoints *3(angles)
  int m = FKOutputDim; // output dimension = # IK joints * 3(positions)

  // tell adolc to trace below computations
  trace_on(adolc_tagID);
    // setup input and output
    vector<adouble> x(n);
    for(int i = 0; i < n; i++)
      x[i] <<= 0.0;
    vector<adouble> y(m);

    // tell adolc we want to use the FK function, y, handlePositions, is the output
    forwardKinematicsFunction(numIKJoints, IKJointIDs, *fk, x, y);

    // get output IK handle positions calculated by adolc
    vector<double> output(m);
    for(int i = 0; i < m; i++)
      y[i] >>= output[i];
  trace_off();
}

/**********************************************************************************/
/*                      Damped Least Squares Implementation                       */
/**********************************************************************************/
void IK::doIK(const Vec3d * targetHandlePositions, Vec3d * jointEulerAngles)
{
  // Students should implement this.
  // Use adolc to evalute the forwardKinematicsFunction and its gradient (Jacobian). It was trained in train_adolc().
  // Specifically, use ::function, and ::jacobian .
  // See ADOLCExample.cpp .
  //
  // Use it implement the Tikhonov IK method (or the pseudoinverse method for extra credit).
  // Note that at entry, "jointEulerAngles" contains the input Euler angles. 
  // Upon exit, jointEulerAngles should contain the new Euler angles.

  // prepare input and output arrays for adol-c
  int numJoints = fk->getNumJoints();
  double input_x_values[FKInputDim], output_y_values[FKOutputDim];
  // for the first pass, use the input jointEulerAngles
  for(int i = 0; i < numJoints; i++)
  {
    input_x_values[3 * i + 0] = jointEulerAngles[i][0];
    input_x_values[3 * i + 1] = jointEulerAngles[i][1];
    input_x_values[3 * i + 2] = jointEulerAngles[i][2];
  }

  // let adol-c evaluate forwardKinematicsFunction for different joint euler angles
  ::function(adolc_tagID, FKOutputDim, FKInputDim, input_x_values, output_y_values);

  // let adol-c evaluate Jacobian
  double jacobianMatrix[FKInputDim * FKOutputDim];
  // create pointer array where each pointer points to one row of the jacobian matrix
  double * jacobianMatrixEachRow[FKOutputDim];
  for(int i = 0; i < FKOutputDim; i++) {
    jacobianMatrixEachRow[i] = & jacobianMatrix[i * FKInputDim];
  }
  // each row is the gradient of one output component of the function
  ::jacobian(adolc_tagID, FKOutputDim, FKInputDim, input_x_values, jacobianMatrixEachRow);

  // get Jacobian matrix (adol-c stores matrix in row-major order)
  MatrixXd J(FKOutputDim, FKInputDim);
  for(int m = 0; m < FKOutputDim; m++)
    for(int n = 0; n < FKInputDim; n++)
    {
      J(m, n) = jacobianMatrix[m * FKInputDim + n];
    }
  // get J^T
  MatrixXd J_transpose(FKInputDim, FKOutputDim);
  J_transpose = J.transpose();
  // creata an nxn identity matrix
  MatrixXd I(FKInputDim, FKInputDim);
  I = MatrixXd::Identity(FKInputDim, FKInputDim);

  // set "punish" param alpha, it avoids theta change to become unreasonabally
  // large
  double alpha = 0.01;

  // now we get matrix A, goal is to solve Ax=b
  // A is square since J is mxn, J^T is nxm -> J^T * J is nxn, also, I is nxn
  MatrixXd A(FKInputDim, FKInputDim);
  A = J_transpose * J + alpha * I;

  // now define the b, rhs
  // make a posChanges vector to store IK handle positions changes
  // since FKOutputDim = 3 * numIKJoints, we can fill in posChanges as below
  VectorXd posChanges(FKOutputDim);
  for(int i = 0; i < numIKJoints; i++)
  {
    // get position change, posChange = targetHandlePos - currentPos
    posChanges[3 * i + 0] = targetHandlePositions[i][0] - output_y_values[3 * i + 0];
    posChanges[3 * i + 1] = targetHandlePositions[i][1] - output_y_values[3 * i + 1];
    posChanges[3 * i + 2] = targetHandlePositions[i][2] - output_y_values[3 * i + 2];
  }

  // make a b with length FKOutputDim and calculate rhs
  VectorXd b(FKOutputDim);
  b = J_transpose * posChanges;

  // now solve the linear system for Euler angles change
  VectorXd x = A.ldlt().solve(b);

  // finally, update jointEulerAngles with solution x (Euler angles change) for next pass
  for(int i = 0; i < numJoints; i++)
  {
    jointEulerAngles[i] += Vec3d(x[3 * i + 0], x[3 * i + 1], x[3 * i + 2]);
  }
}

// below is the extra implementation for doIK with pseudo inverse method
// uncomment it and comment out above doIk() if you'd like to test it

/**********************************************************************************/
/*                    Pseudo-inverse Method Implementation                        */
/**********************************************************************************/
// void IK::doIK(const Vec3d * targetHandlePositions, Vec3d * jointEulerAngles)
// {
//   // Students should implement this.
//   // Use adolc to evalute the forwardKinematicsFunction and its gradient (Jacobian). It was trained in train_adolc().
//   // Specifically, use ::function, and ::jacobian .
//   // See ADOLCExample.cpp .
//   //
//   // Use it implement the Tikhonov IK method (or the pseudoinverse method for extra credit).
//   // Note that at entry, "jointEulerAngles" contains the input Euler angles. 
//   // Upon exit, jointEulerAngles should contain the new Euler angles.

//   // prepare input and output arrays for adol-c
//   int numJoints = fk->getNumJoints();
//   double input_x_values[FKInputDim], output_y_values[FKOutputDim];
//   // for the first pass, use the input jointEulerAngles
//   for(int i = 0; i < numJoints; i++)
//   {
//     input_x_values[3 * i + 0] = jointEulerAngles[i][0];
//     input_x_values[3 * i + 1] = jointEulerAngles[i][1];
//     input_x_values[3 * i + 2] = jointEulerAngles[i][2];
//   }

//   // let adol-c evaluate forwardKinematicsFunction for different joint euler angles
//   ::function(adolc_tagID, FKOutputDim, FKInputDim, input_x_values, output_y_values);

//   // let adol-c evaluate Jacobian
//   double jacobianMatrix[FKInputDim * FKOutputDim];
//   // create pointer array where each pointer points to one row of the jacobian matrix
//   double * jacobianMatrixEachRow[FKOutputDim];
//   for(int i = 0; i < FKOutputDim; i++) {
//     jacobianMatrixEachRow[i] = & jacobianMatrix[i * FKInputDim];
//   }
//   // each row is the gradient of one output component of the function
//   ::jacobian(adolc_tagID, FKOutputDim, FKInputDim, input_x_values, jacobianMatrixEachRow);

//   // get Jacobian matrix (adol-c stores matrix in row-major order)
//   MatrixXd J(FKOutputDim, FKInputDim);
//   for(int m = 0; m < FKOutputDim; m++)
//     for(int n = 0; n < FKInputDim; n++)
//     {
//       J(m, n) = jacobianMatrix[m * FKInputDim + n];
//     }
//   // get J^T
//   MatrixXd J_transpose(FKInputDim, FKOutputDim);
//   J_transpose = J.transpose();

//   // calculate J_dagger with pseudo inverse in process, it has the same dimension with J
//   MatrixXd J_dagger(FKInputDim, FKOutputDim);
//   J_dagger = J_transpose * (J * J_transpose).inverse();

//   // creata an nxn identity matrix as A, since lfh only has delta theta
//   MatrixXd A(FKInputDim, FKInputDim);
//   A = MatrixXd::Identity(FKInputDim, FKInputDim);

//   // now define the b, rhs
//   // make a posChanges vector to store IK handle positions changes
//   // since FKOutputDim = 3 * numIKJoints, we can fill in posChanges as below
//   VectorXd posChanges(FKOutputDim);
//   for(int i = 0; i < numIKJoints; i++)
//   {
//     // get position change, posChange = targetHandlePos - currentPos
//     posChanges[3 * i + 0] = targetHandlePositions[i][0] - output_y_values[3 * i + 0];
//     posChanges[3 * i + 1] = targetHandlePositions[i][1] - output_y_values[3 * i + 1];
//     posChanges[3 * i + 2] = targetHandlePositions[i][2] - output_y_values[3 * i + 2];
//   }

//   // make a b with length FKOutputDim and calculate rhs
//   VectorXd b(FKOutputDim);
//   b = J_dagger * posChanges;

//   // now solve the linear system for Euler angles change
//   VectorXd x = A.ldlt().solve(b);

//   // finally, update jointEulerAngles with solution x (Euler angles change) for next pass
//   for(int i = 0; i < numJoints; i++)
//   {
//     jointEulerAngles[i] += Vec3d(x[3 * i + 0], x[3 * i + 1], x[3 * i + 2]);
//   }
// }