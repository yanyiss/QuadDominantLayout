#include "deformer.h"
namespace HDP
{
	void func_callback(const size_t N, const std::vector<double>& x, double& f, std::vector<double>& g, void* user_supply)
	{
		auto ptr_this = static_cast<deformer*>(user_supply);
		ptr_this->field_energy_evaluation(x, f, g);
	}

	void deformer::fill_rotation_matrix_dimension(const Matrix3Xd& crossfield, int x_index, int y_index, MatrixXd& rotation_matrix)
	{
		rotation_matrix.resize(dimension, dimension);
		rotation_matrix.setIdentity();
		Vector3d z = crossfield.col(x_index).cross(crossfield.col(y_index));
		for (int i = 0; i < 3; ++i)
		{
			rotation_matrix(i, 0) = crossfield(i, x_index);
			rotation_matrix(i, 1) = crossfield(i, y_index);
			rotation_matrix(i, 2) = z(i);
		}
	}

	void deformer::fill_LC_connection_rotation_matrix(HalfedgeHandle h, MatrixXd& rotation_matrix)
	{
		//rotation_matrix rotate h.opp().face() to the plane of h.face()
		rotation_matrix.resize(dimension, dimension);
		rotation_matrix.setIdentity();
		Vec3d fn = mesh->calc_face_normal(mesh->face_handle(h));
		Vec3d gn = mesh->calc_face_normal(mesh->opposite_face_handle(h));
		Vec3d en = mesh->calc_edge_vector(h).normalized();
		double cos_theta = fn.dot(gn);
		double sin_theta = sin(conservativeArcCos(fn.dot(gn)) * (fn.cross(gn).dot(en) > 0 ? -1.0 : 1.0));
		Matrix3d nnT; nnT << en[0] * en[0], en[0] * en[1], en[0] * en[2],
			en[1] * en[0], en[1] * en[1], en[1] * en[2],
			en[2] * en[0], en[2] * en[1], en[2] * en[2];
		Matrix3d anti_symmetric; anti_symmetric << 0.0, -en[2], en[1], en[2], 0.0, -en[0], -en[1], en[0], 0.0;
		//Rodrigues' rotation formula
		rotation_matrix.block(0, 0, 3, 3) = cos_theta * Matrix3d::Identity() + (1 - cos_theta) * nnT + sin_theta * anti_symmetric;
	}

	void deformer::compute_rotation_matrix(HalfedgeHandle h, MatrixXd& rotation_matrix)
	{
		MatrixXd mat_i, mat_j;
		int fid4 = mesh->face_handle(h).idx() * 4;
		int gid4 = mesh->opposite_face_handle(h).idx() * 4;
		auto& matching = cf->getMatching();
		fill_rotation_matrix_dimension(cf->getCrossField(), fid4, another_layer(gid4, matching[h.idx()]), mat_i);
		fill_rotation_matrix_dimension(cf->getCrossField(), gid4, another_layer(fid4, 4 - matching[h.idx()]), mat_j);
		fill_LC_connection_rotation_matrix(h, rotation_matrix);
		rotation_matrix *= mat_j * mat_i.inverse();
	}

	void deformer::field_energy_evaluation(const std::vector<double>& x, double& f, std::vector<double>& g)
	{
		f = 0;
		g.clear();
		g.resize(dimension * dimension * mesh->n_faces(), 0);

		for (auto tf : mesh->faces())
		{
			int fid = tf.idx();
			int fid_plus = fid * dimension * dimension;
			for (auto tfh : mesh->fh_range(tf))
			{
				if (tfh.edge().is_boundary())
					continue;
				int gid = tfh.opp().face().idx();
				int gid_plus = gid * dimension * dimension;
				MatrixXd F;
				compute_rotation_matrix(tfh, F);
				MatrixXd F_inv = F.inverse();
				for (int i = 0; i < dimension; ++i)
				{
					for (int j = 0; j < dimension; ++j)
					{
						int fpidj = fid_plus + i * dimension + j;
						double sum = x[fpidj];
						for (int k = 0; k < dimension; ++k)
						{
							sum -= x[gid_plus + i * dimension + k] * F(k, j);
						}
						f += sum * sum;

						g[fpidj] += sum;
						for (int p = 0; p < dimension; ++p)
						{
							double sum = 0;
							for (int k = 0; k < dimension; ++k)
							{
								sum += x[fid_plus + i * dimension + k] * F_inv(k, p);
							}
							g[fpidj] += F_inv(j, p) * (sum - x[gid_plus + i * dimension + p]);
						}
					}
				}
			}
			for (int i = 0; i < dimension; ++i)
			{
				for (int j = 0; j < dimension; ++j)
				{
					int fpidj = fid_plus + i * dimension + j;
					double sum = 0;
					for (int k = 0; k < dimension; ++k)
					{
						sum += x[fid_plus + i * dimension + k] * x[fid_plus + j * dimension + k];
					}
					f += w * sum * sum - 2.0 * w * x[fpidj] * x[fpidj];

					for (int p = 0; p < dimension; ++p)
					{
						if (p != i)
						{
							double sum = 0;
							for (int k = 0; k < dimension; ++k)
							{
								sum += x[fid_plus + i * dimension + k] * x[fid_plus + p * dimension + k];
							}
							g[fpidj] += w * x[fid_plus + p * dimension + j] * sum;
						}
						else
						{
							//g[fpidj] += 2.0 * w * x[fpidj] * x[fid_plus + i * dimension + p] * x[fid_plus + i * dimension + p];
							double sum = 0;
							for (int k = 0; k < dimension; ++k)
							{
								sum += x[fid_plus + i * dimension + k] * x[fid_plus + i * dimension + k];
							}
							g[fpidj] += 2.0 * w * x[fpidj] * sum;
						}
					}
					g[fpidj] -= 2.0 * w * x[fpidj];
					g[fpidj] *= 2.0;

				}
				f += w;
				/*for (int k = 0; k < dimension; ++k)
				{
					f -= 2 * w * x[fid_plus + i * dimension + k] * x[fid_plus + i * dimension + k];
				}*/
			}
		}
	}
	
	void deformer::update_field()
	{
		cf->initMeshInfo();
		cf->setBoundaryConstraint();
		cf->setField();
		cf->initFieldInfo();
		//cf->write_field();
	}

	void deformer::compute_rotation()
	{
		auto& crossfield = cf->getCrossField();
		HLBFGS solver;
		int vn = dimension * dimension * mesh->n_faces();
		solver.set_number_of_variables(vn);
		solver.set_verbose(true);
		solver.set_func_callback(func_callback, 0, 0, 0, 0);
		std::vector<double> v;
		v.resize(vn, 0);
		for (int i = 0; i < mesh->n_faces(); ++i)
		{
			for (int j = 0; j < dimension; ++j)
			{
				for (int k = 0; k < dimension; ++k)
				{
					v[i * dimension * dimension + j * dimension + k] = (j == k) ? 1.2 : 0.2;
				}
			}
		}
		solver.optimize_without_constraints(&v[0], 100, this);
		for (int i = 0; i < mesh->n_faces(); ++i)
		{
			for (int j = 0; j < dimension; ++j)
			{
				for (int k = 0; k < dimension; ++k)
				{
					rotation_on_triangle[i](j, k) = v[i * dimension * dimension + j * dimension + k];
				}
			}
			JacobiSVD<MatrixXd> svd(rotation_on_triangle[i], ComputeFullU | ComputeFullV);
			rotation_on_triangle[i] = svd.matrixU() * svd.matrixV().transpose();
			for (int j = 0; j < 4; ++j)
			{
				crossfield.col(i * 4 + j) = rotation_on_triangle[i] * crossfield.col(i * 4 + j);
			}
		}
		for (auto& si : cf->getSingularity())
		{
			for (auto& s : si)
			{
				for (auto tf : mesh->vf_range(mesh->vertex_handle(s)))
				{
					dprint(rotation_on_triangle[tf.idx()]);
					dprint();
				}
			}
		}
	}

	void deformer::compute_position()
	{
		double fixed_weight = std::accumulate(triangle_area.begin(), triangle_area.end(), 0.0);
		fixed_weight *= 12.0 / mesh->n_faces();

		std::vector<Triplet<double>> triplets; triplets.reserve(mesh->n_edges() * 2 + mesh->n_vertices());
		MatrixX3d right; right.resize(mesh->n_vertices(), 3);
		for (auto tv : mesh->vertices())
		{
			dprint(tv.idx());
			int vid = tv.idx();
			Vector3d tvpos = Vec3d_to_Vector3d(mesh->point(tv));
			double sum = 0;
			Vector3d weight_vec; weight_vec.setZero();
			for (auto tvoh : mesh->voh_range(tv))
			{
				if (tvoh.edge().is_boundary())
					continue;
				int fid = tvoh.face().idx();
				int gid = tvoh.opp().face().idx();
				double weight = triangle_area[fid] + triangle_area[gid];
				triplets.emplace_back(vid, tvoh.to().idx(), -weight);
				sum += weight;
				weight_vec += (triangle_area[fid] * rotation_on_triangle[fid] + triangle_area[gid] * rotation_on_triangle[gid])
					* (tvpos - Vec3d_to_Vector3d(mesh->point(tvoh.to())));
			}
			if (vid == 0)
			{
				sum += fixed_weight;
				weight_vec += fixed_weight * tvpos;
			}
			triplets.emplace_back(vid, vid, sum);
			right.row(vid) << weight_vec(0), weight_vec(1), weight_vec(2);
		}

		SparseMatrix<double> laplace; laplace.resize(mesh->n_vertices(), mesh->n_vertices());
		laplace.setFromTriplets(triplets.begin(), triplets.end());
		SimplicialCholesky<SparseMatrix<double>> solver;
		solver.compute(laplace);
		right = solver.solve(right);
		for (auto tv : mesh->vertices())
		{
			int vid = tv.idx();
			mesh->set_point(tv, OpenMesh::Vec3d(right(vid, 0), right(vid, 1), right(vid, 2)));
			dprint(mesh->point(tv));
		}
	}

	void deformer::run()
	{
		update_field();
		update_triangle_area();
		compute_rotation();
		//compute_position();
	}
}