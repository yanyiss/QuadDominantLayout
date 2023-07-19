#pragma once
#include <QObject>
#include "crossField.h"
#include "../src/Dependency/HLBFGS/HLBFGS.h"
namespace HDP
{
	void func_callback(const size_t N, const std::vector<double>& x, double& f, std::vector<double>& g, void* user_supply);

	using namespace Eigen;
	//template <int dimension>
	//int dimension = 3;
	class deformer : public QObject
	{
		Q_OBJECT
	public:
		deformer(Mesh& m) 
		{
			mesh = &m;
			rotation_on_triangle.resize(mesh->n_faces());
			triangle_area.resize(mesh->n_faces());
			cf = new crossField(mesh);
		}
		~deformer(){ if (cf) { delete cf; cf = nullptr; } }

	public:
		Mesh* mesh;
		crossField* cf = nullptr;
		double w = 100.0;
		int dimension = 3;
		std::vector<Matrix3d> rotation_on_triangle;
		std::vector<double> triangle_area;

		

		inline Vector3d Vec3d_to_Vector3d(OpenMesh::Vec3d& pos) { return Vector3d(pos[0], pos[1], pos[2]); }
		void update_triangle_area() { for (auto tf : mesh->faces()) { triangle_area[tf.idx()] = mesh->calc_face_area(tf); } }

		inline int another_layer(int current_layer, int plus)
		{
			return ((current_layer % 4) + plus) > 3 ? current_layer + plus - 4 : current_layer + plus;
		}

		void fill_rotation_matrix_dimension(const Matrix3Xd& crossfield, int x_index, int y_index, MatrixXd& rotation_matrix);

		void fill_LC_connection_rotation_matrix(HalfedgeHandle h, MatrixXd& rotation_matrix);

		void compute_rotation_matrix(HalfedgeHandle h, MatrixXd& rotation_matrix);

		void field_energy_evaluation(const std::vector<double>& x, double& f, std::vector<double>& g);
		;
		void update_field();

		void compute_rotation();

		void compute_position();

		void run();
	};


}

