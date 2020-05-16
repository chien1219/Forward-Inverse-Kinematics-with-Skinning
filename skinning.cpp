#include "skinning.h"
#include "vec3d.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
#include <Eigen/Geometry> 

using namespace std;

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

void Skinning::applySkinning(const RigidTransform4d * jointSkinTransforms, double * newMeshVertexPositions, FK * fk) const
{
	// Choose Skinning Method
	LinearBlendingSkinning(jointSkinTransforms, newMeshVertexPositions);
	//DualQuaternionsSkinning(jointSkinTransforms, newMeshVertexPositions, fk);
}

void Skinning::LinearBlendingSkinning(const RigidTransform4d * jointSkinTransforms, double * newMeshVertexPositions) const
{
	// pi = ∑ wj* Mj * p'i
	
	for (int vtxID = 0; vtxID < numMeshVertices; vtxID++)
	{
		Mat4d sumMatrix;
		Vec4d restPos(restMeshVertexPositions[3 * vtxID + 0], restMeshVertexPositions[3 * vtxID + 1], restMeshVertexPositions[3 * vtxID + 2], 1);
		Vec4d newPos(0, 0, 0, 0);

		for (int j = 0; j < numJointsInfluencingEachVertex; j++)
		{
			int index = vtxID * numJointsInfluencingEachVertex + j;
			int infJointID = meshSkinningJoints[index];
			newPos += meshSkinningWeights[index] * jointSkinTransforms[infJointID] * restPos;
		}
		newMeshVertexPositions[3 * vtxID + 0] = newPos[0];
		newMeshVertexPositions[3 * vtxID + 1] = newPos[1];
		newMeshVertexPositions[3 * vtxID + 2] = newPos[2];
	}
}

void Skinning::DualQuaternionsSkinning(const RigidTransform4d * jointSkinTransforms, double * newMeshVertexPositions, FK * fk) const
{
	int N = fk->getNumJoints();

	vector<Eigen::Quaterniond> getJointQuaternion(N);
	vector<Vec3d> getJointTranslationVec(N);
	vector<Vec3d> getNewTranslationVec(numMeshVertices);
	vector<Eigen::Quaterniond> getNewQuaternion(numMeshVertices);

	// Convert rigid joint transform matrices to quaternion  Mat3d -> Eigen::Matrix3d
	for (int i = 0; i < N; i++)
	{
		// Get translation vectors
		getJointTranslationVec[i] = jointSkinTransforms->getTranslation();

		Eigen::Matrix3d getTransformMat3d;
		Mat3d linearMat3d = jointSkinTransforms[i].getLinearTrans();

		for (int rowID = 0; rowID < 3; rowID++)
			for (int colID = 0; colID < 3; colID++)
				getTransformMat3d(rowID, colID) = linearMat3d[rowID][colID];

		getJointQuaternion[i] = getTransformMat3d;
		getJointQuaternion[i].normalized();
	}

	for (int i = 0; i < numMeshVertices; i++)
	{
		// Quaternions
		getNewQuaternion[i] = Eigen::Quaterniond(0.0, 0.0, 0.0, 0.0);
		// Translation vectors
		getNewTranslationVec[i] = Vec3d(0.0, 0.0, 0.0);

		for (int j = 0; j < numJointsInfluencingEachVertex; j++)
		{
			// Compute new quaternions
			getNewQuaternion[i].coeffs() += getJointQuaternion[meshSkinningJoints[i * numJointsInfluencingEachVertex + j]].coeffs() * meshSkinningWeights[i * numJointsInfluencingEachVertex + j];

			// Compute new translation vector
			getNewTranslationVec[i] += meshSkinningWeights[i * numJointsInfluencingEachVertex + j] * getJointTranslationVec[meshSkinningJoints[i * numJointsInfluencingEachVertex + j]];
		}
		getNewQuaternion[i].normalized();
	}

	// Convert quaternions to rotation matrices  Eigen::Matrix3d -> Mat3d
	for (int i = 0; i < numMeshVertices; i++)
	{
		Eigen::Matrix3d getNewRotationMat3d = getNewQuaternion[i].toRotationMatrix();
		Mat3d newRotationMat;
		RigidTransform4d newTransformMat;

		for (int rowID = 0; rowID < 3; rowID++)
			for (int colID = 0; colID < 3; colID++)
				newRotationMat[rowID][colID] = getNewRotationMat3d(rowID, colID);

		newTransformMat = AffineTransform4d(newRotationMat, getNewTranslationVec[i]);

		Vec4d restVertexMeshVec4d = Vec4d(restMeshVertexPositions[3 * i + 0], restMeshVertexPositions[3 * i + 1], restMeshVertexPositions[3 * i + 2], 1.0);
		Vec4d newVertexMeshVec4d = newTransformMat * restVertexMeshVec4d;

		for (int j = 0; j < 3; j++)
			newMeshVertexPositions[3 * i + j] = newVertexMeshVec4d[j];
	}
}