#include "skinning.h"
#include "vec3d.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
using namespace std;
using namespace Eigen;

// CSCI 520 Computer Animation and Simulation
// Jernej Barbic and Yijing Li

Skinning::Skinning(int numMeshVertices, const double * restMeshVertexPositions,
    const std::string & meshSkinningWeightsFilename)
{
  this->numMeshVertices = numMeshVertices;
  this->restMeshVertexPositions = restMeshVertexPositions;

  cout << "Loading skinning weights..." << endl;
  ifstream fin(meshSkinningWeightsFilename.c_str());
  assert(fin);
  int numWeightMatrixRows = 0, numWeightMatrixCols = 0;
  fin >> numWeightMatrixRows >> numWeightMatrixCols;
  assert(fin.fail() == false);
  assert(numWeightMatrixRows == numMeshVertices);
  int numJoints = numWeightMatrixCols;

  vector<vector<int>> weightMatrixColumnIndices(numWeightMatrixRows);
  vector<vector<double>> weightMatrixEntries(numWeightMatrixRows);
  fin >> ws;
  while(fin.eof() == false)
  {
    int rowID = 0, colID = 0;
    double w = 0.0;
    fin >> rowID >> colID >> w;
    weightMatrixColumnIndices[rowID].push_back(colID);
    weightMatrixEntries[rowID].push_back(w);
    assert(fin.fail() == false);
    fin >> ws;
  }
  fin.close();

  // Build skinning joints and weights.
  numJointsInfluencingEachVertex = 0;
  for (int i = 0; i < numMeshVertices; i++)
    numJointsInfluencingEachVertex = std::max(numJointsInfluencingEachVertex, (int)weightMatrixEntries[i].size());
  assert(numJointsInfluencingEachVertex >= 2);

  // Copy skinning weights from SparseMatrix into meshSkinningJoints and meshSkinningWeights.
  meshSkinningJoints.assign(numJointsInfluencingEachVertex * numMeshVertices, 0);
  meshSkinningWeights.assign(numJointsInfluencingEachVertex * numMeshVertices, 0.0);
  for (int vtxID = 0; vtxID < numMeshVertices; vtxID++)
  {
    vector<pair<double, int>> sortBuffer(numJointsInfluencingEachVertex);
    for (size_t j = 0; j < weightMatrixEntries[vtxID].size(); j++)
    {
      int frameID = weightMatrixColumnIndices[vtxID][j];
      double weight = weightMatrixEntries[vtxID][j];
      sortBuffer[j] = make_pair(weight, frameID);
    }
    sortBuffer.resize(weightMatrixEntries[vtxID].size());
    assert(sortBuffer.size() > 0);
    sort(sortBuffer.rbegin(), sortBuffer.rend()); // sort in descending order using reverse_iterators
    for(size_t i = 0; i < sortBuffer.size(); i++)
    {
      meshSkinningJoints[vtxID * numJointsInfluencingEachVertex + i] = sortBuffer[i].second;
      meshSkinningWeights[vtxID * numJointsInfluencingEachVertex + i] = sortBuffer[i].first;
    }

    // Note: When the number of joints used on this vertex is smaller than numJointsInfluencingEachVertex,
    // the remaining empty entries are initialized to zero due to vector::assign(XX, 0.0) .
  }
}

/**********************************************************************************/
/*                    Linear Blend Skinning Implementation                        */
/**********************************************************************************/
// void Skinning::applySkinning(const RigidTransform4d * jointSkinTransforms, double * newMeshVertexPositions) const
// {
//   // Students should implement this
//   // Formula: currNewVertPosVec = sum_overRelaventJoints(jointWight_j * jointTrans * restMeshVertPos)
//   Vec4d currRestMeshVertPosVec, currNewVertPosVec(0,0,0,0);
//   int currInd;
//   for(int i=0; i<numMeshVertices; i++)
//   {
//     // get restMeshVertPos. Elevated from actual caculation for loop beloww
//     // b/c they remain the same
//     currRestMeshVertPosVec = Vec4d(restMeshVertexPositions[3 * i + 0],
//                                    restMeshVertexPositions[3 * i + 1],
//                                    restMeshVertexPositions[3 * i + 2],
//                                    1);

//     // calculate current joint position based on the formula
//     for(int j=0; j<numJointsInfluencingEachVertex; j++)
//     {
//       currInd = numJointsInfluencingEachVertex * i + j;
//       currNewVertPosVec +=   meshSkinningWeights[currInd]
//                            * jointSkinTransforms[meshSkinningJoints[currInd]]
//                            * currRestMeshVertPosVec;
//     }

//     newMeshVertexPositions[3 * i + 0] = currNewVertPosVec[0];
//     newMeshVertexPositions[3 * i + 1] = currNewVertPosVec[1];
//     newMeshVertexPositions[3 * i + 2] = currNewVertPosVec[2];

//     // clear currNewVertPosVec before next joint weight calculationg
//     currNewVertPosVec = Vec4d(0,0,0,0);
//   }
// }

/**********************************************************************************/
/*                    Dual Quaternion Skinning Implementation                     */
/**********************************************************************************/
void Skinning::applySkinning(const RigidTransform4d * jointSkinTransforms, double * newMeshVertexPositions) const
{
  // Students should implement this
  // Formula: currNewVertPosVec = sum_overRelaventJoints(jointWight_j * dual_quaternion(q0,q1))

  int currInd;
  Mat3d currRotation;
  Vec3d currTranslation;
  for(int i=0; i<numMeshVertices; i++)
  {
    // get restMeshVertPos (x)
    Vector3d currRestMeshVertPosVec = {restMeshVertexPositions[3 * i + 0],
                                       restMeshVertexPositions[3 * i + 1],
                                       restMeshVertexPositions[3 * i + 2]};

    Quaterniond currNewQ0(0,0,0,0), currNewQ1(0,0,0,0);
    // calculate current joint position based on the formula
    for(int j=0; j<numJointsInfluencingEachVertex; j++)
    {
      currInd = numJointsInfluencingEachVertex * i + j;

      // extract rotation and translation
      currRotation = jointSkinTransforms[meshSkinningJoints[currInd]].getLinearTrans();
      currTranslation = jointSkinTransforms[meshSkinningJoints[currInd]].getTranslation();

      // form q0
      // convert rotation to eigen matrix for quaternion calculation
      Matrix3d currRotationEigenMat(3,3);
      for(int rowID = 0; rowID < 3; rowID++)
        for(int colID = 0; colID < 3; colID++)
          currRotationEigenMat(rowID,colID) = currRotation[rowID][colID];
      Quaterniond q0 = Quaterniond(currRotationEigenMat);

      // check if angle < 0
      if(q0.dot(currNewQ0) < 0) {
        q0 = q0.coeffs() * -1;
      }

      // form q1
      Quaterniond t = Quaterniond(0, currTranslation[0], currTranslation[1], currTranslation[2]);
      t = t.coeffs() * 0.5;
      Quaterniond q1 = t * q0;

      // summing up dual quaternion
      currNewQ0.coeffs() += meshSkinningWeights[currInd] * q0.coeffs();
      currNewQ1.coeffs() += meshSkinningWeights[currInd] * q1.coeffs();
    }

    Quaterniond c0, c_eps;
    c0 = currNewQ0.coeffs() / currNewQ0.norm();
    c_eps = currNewQ1.coeffs() / currNewQ0.norm();

    // calculate final rotation and translation
    Matrix3d currFinalRotation = c0.toRotationMatrix();
    Vector3d currFinalTranslation = {2 * (-c_eps.coeffs()[0]*c_eps.coeffs()[1] + c_eps.coeffs()[1]*c_eps.coeffs()[0] - c_eps.coeffs()[2]*c_eps.coeffs()[3] + c_eps.coeffs()[3]*c_eps.coeffs()[2]),
                                     2 * (-c_eps.coeffs()[0]*c_eps.coeffs()[2] + c_eps.coeffs()[1]*c_eps.coeffs()[3] + c_eps.coeffs()[2]*c_eps.coeffs()[0] - c_eps.coeffs()[3]*c_eps.coeffs()[1]),
                                     2 * (-c_eps.coeffs()[0]*c_eps.coeffs()[3] - c_eps.coeffs()[1]*c_eps.coeffs()[2] + c_eps.coeffs()[2]*c_eps.coeffs()[1] + c_eps.coeffs()[3]*c_eps.coeffs()[0])};

    // calculate new vertex position
    Vector3d currNewVertPosVec = currFinalTranslation + currFinalRotation * currRestMeshVertPosVec;

    newMeshVertexPositions[3 * i + 0] = currNewVertPosVec[0];
    newMeshVertexPositions[3 * i + 1] = currNewVertPosVec[1];
    newMeshVertexPositions[3 * i + 2] = currNewVertPosVec[2];
  }
}