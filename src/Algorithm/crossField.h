#pragma once
#include "..\MeshViewer\MeshDefinition.h"
#include <Eigen\Sparse>

#include "..\Toolbox\Math\mathFunctions.h"
#include "..\Toolbox\stream\dprint.h"
#include "..\Toolbox\stream\filesOperator.h"
#include "..\Toolbox\Math\BoolVector.h"
#include "StatisticsMostValues.h"
/*
generate crossfield for quad layout
should satisfy the following properties
1. align with principal direction at least in most regions
2. singularities should not appear in where smaller principal curvature is almost zero
3. smooth in most regions but at singularities
4. may decided later in umblic regions

*/
namespace HDP
{
	class crossField
	{
	public:
		crossField(){}
		crossField(Mesh* mesh_, std::string file_name_ = "yanyisheshou") : mesh(mesh_), file_name(file_name_) {};
		crossField(const crossField& cf)
		{
			mesh = cf.mesh;
			file_name = cf.file_name;
			position = cf.position; normal = cf.normal; faceBase = cf.faceBase;
			constraintId = cf.constraintId; constraintVector = cf.constraintVector;
			crossfield = cf.crossfield; matching = cf.matching; singularity = cf.singularity;
			updateField = cf.updateField; cur = cf.cur;
		}
		~crossField() {};

	public:
		Mesh* mesh = nullptr;
		std::string file_name;
		typedef std::complex<double>  COMPLEX;


		Eigen::Matrix3Xd position; const Eigen::Matrix3Xd& getPosition() { return position; } void setPosition();
		Eigen::Matrix3Xd normal;   const Eigen::Matrix3Xd& getNormal() { return normal; }     void setNormal();
		Eigen::Matrix3Xd faceBase; const Eigen::Matrix3Xd& getFaceBase() { return faceBase; } void setFaceBase();
		void initMeshInfo();

		std::vector<int> constraintId;
		Eigen::Matrix3Xd constraintVector;
		void setCurvatureConstraint();
		void setBoundaryConstraint();
		void setOuterConstraint(BoolVector& cons_flag, Eigen::Matrix3Xd& cons_direction, bool curvature_constraint = true);

		Eigen::Matrix3Xd crossfield;  Eigen::Matrix3Xd& getCrossField() { return crossfield; }   void setField();

		std::vector<int> matching;    std::vector<int>& getMatching() { return matching; }       void setMatching();
		std::vector<std::vector<int>> singularity; std::vector<std::vector<int>>& getSingularity() { return singularity; } void setSingularity();
		void initFieldInfo(); bool updateField = false;

		void read_field();
		void write_field();

		double axisDif(HalfedgeHandle hh);
		std::vector<std::pair<double, double>> cur;
	};
}


