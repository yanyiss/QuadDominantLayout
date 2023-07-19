#include <iostream>
#include <algorithm>
#include <cmath>

#include <qapplication.h>

#include <OpenMesh/Core/Utils/vector_cast.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include "MeshViewerWidget.h"
#include "..\src\Dependency\Common\CommonDefinitions.h"

using namespace Qt;

MeshViewerWidget::MeshViewerWidget(QWidget* parent)
 : QGLViewerWidget(parent)
{
	mesh.request_vertex_status();
	mesh.request_edge_status();
	mesh.request_face_status();

	mesh.request_face_normals();
	mesh.request_vertex_normals();

	mesh_vector.clear(); mesh_vector_index = -1;

	draw_BBox_OK = false;
	draw_mesh_boundary_ok = false;
	first_init = true;
}

MeshViewerWidget::MeshViewerWidget(QGLFormat& _fmt, QWidget* _parent)
: QGLViewerWidget(_fmt, _parent)
{
	mesh.request_vertex_status();
	mesh.request_edge_status();
	mesh.request_face_status();

	mesh.request_face_normals();
	mesh.request_vertex_normals();

	mesh_vector.clear(); mesh_vector_index = -1;

	draw_BBox_OK = false;
	draw_mesh_boundary_ok = false;
	first_init = true;
}

MeshViewerWidget::~MeshViewerWidget()
{
}

void MeshViewerWidget::updateMeshCenter()
{
	typedef Mesh::Point Point;
	Mesh::VertexIter vIt = mesh.vertices_begin();
	Mesh::VertexIter vEnd = mesh.vertices_end();
	bbMin = bbMax = OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh.point(vIt));

	size_t count = 0;
	for(; vIt != vEnd; ++vIt, ++count)
	{
		bbMin.minimize( OpenMesh::vector_cast<OpenMesh::Vec3d>( mesh.point(vIt) ) );
		bbMax.maximize( OpenMesh::vector_cast<OpenMesh::Vec3d>( mesh.point(vIt) ) );
	}

	Mesh::EdgeIter e_it = mesh.edges_begin();
	Mesh::EdgeIter e_end = mesh.edges_end();
	double aveLen = 0.0; double maxLen = 0.0; double minLen = mesh.calc_edge_length(e_it);
	double e_len = 0.0;
	for(; e_it != e_end; ++e_it)
	{
		double e_len = mesh.calc_edge_length(e_it);
		if( e_len > maxLen )
		{
			maxLen = e_len;
		}
		else if(e_len < minLen )
		{
			minLen = e_len;
		}
		aveLen += e_len;
	}

	// set center and radius and box's radius.
	/*OpenMesh::Vec3d c = (bbMin + bbMax)*0.5;
	for (vIt = mesh.vertices_begin(); vIt != vEnd; ++vIt, ++count)
	{
		OpenMesh::Vec3d p = mesh.point(vIt) - c;
		mesh.set_point(vIt, p);
	}*/

	if( first_init )
	{
		set_scene_pos( (bbMin+bbMax)*0.5, (bbMin-bbMax).norm()*0.5 );
		first_init = false;
	}
	else
	{
		set_scene_pos((bbMin + bbMax)*0.5, (bbMin - bbMax).norm()*0.5);
	}
	
	
	printf("BoundingBox:\nX : [ %f , %f ]\n",bbMin[0],bbMax[0]);
	printf("Y : [ %f , %f ]\n",bbMin[1],bbMax[1]);
	printf("Z : [ %f , %f ]\n",bbMin[2],bbMax[2]);
	printf("Diag length of BBox : %f\n", (bbMax-bbMin).norm() );
	printf("Edge Length : Max : %f; Min : %f; AVG : %f\n",maxLen,minLen, aveLen/mesh.n_edges() );
}

void MeshViewerWidget::updateMeshNormals()
{
	mesh.update_face_normals();
	mesh.update_vertex_normals();
}

bool MeshViewerWidget::openMesh(const char* filename)
{
	file = filename;
	clearAllMesh();
	bool read_OK = OpenMesh::IO::read_mesh( mesh, filename );

	printf("%s\n", filename);
	if ( read_OK )
	{
		initMesh();
		/*Mesh::EdgeHandle eh0 = mesh.edge_handle(960);
		flip_openmesh(eh0, mesh);
		Mesh::EdgeHandle eh1 = mesh.edge_handle(85);
		flip_openmesh(eh1, mesh);*/
		/*FILE* f_de = fopen("b.de", "w");
		for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
		{
			if (mesh.is_boundary(v_it))
			{
				OpenMesh::Vec3d& p = mesh.point(v_it);
				fprintf(f_de, "%d %15.14f %15.14f\n", v_it.handle().idx(), p[0], p[1]);
			}
		}
		fclose(f_de);*/
		// loading done
		mesh_vector.push_back(mesh); mesh_vector_index = 0;
		return true;
	}
	return false;
}

void MeshViewerWidget::initMesh()
{
	mesh.request_vertex_status();
	mesh.request_edge_status();
	mesh.request_face_status();

	mesh.request_face_normals();
	mesh.request_vertex_normals();
	printBasicMeshInfo();
	updateMesh();
	ifUpdateMesh = true;
	flag = true;
	regularflag = true;
	if_has_field = false;
	if_has_ppl = false;
}

void MeshViewerWidget::printBasicMeshInfo()
{
	if (mesh.n_vertices() == 0)
		printf("No Mesh\n");

	checkMeshMode();

	QString meshType;
	switch(meshMode())
	{
	case TRIANGLE:
		printf("Triangle Mesh.\n");
		break;
	case QUAD:
		printf("Quadrilateral Mesh.\n");
		break;
	default:
		printf("General Mesh.\n");
		break;
	}

	printf("Information of the input mesh:\nVertex : %d;\nFace : %d;\nEdge : %d, HalfEdge : %d\n",
		mesh.n_vertices(),mesh.n_faces(),mesh.n_edges(),mesh.n_halfedges());;

	//
	/*for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		OpenMesh::Vec3d& p = mesh.point(v_it);
		p[2] = 0.0;
	}*/

	//save tet mesh
	/*int nf = mesh.n_faces(); int nv = mesh.n_vertices(); int ne = mesh.n_edges();
	std::vector<OpenMesh::Vec3d> face_n(nf);
	std::vector<double> edge_len(ne);
	for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		std::vector<OpenMesh::Vec3d> one_face;
		for (Mesh::FaceVertexIter fv_it = mesh.fv_iter(f_it); fv_it; ++fv_it)
		{
			one_face.push_back(mesh.point(fv_it));
		}
		OpenMesh::Vec3d n = OpenMesh::cross(one_face[1] - one_face[0], one_face[2] - one_face[0]);
		face_n[f_it->idx()] = n.normalize();
	}
	for (Mesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
	{
		edge_len[e_it->idx()] = mesh.calc_edge_length(e_it);
	}
	std::vector<OpenMesh::Vec3d> v_v(nv);
	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		double ave_edge_len = 0.0; OpenMesh::Vec3d vn(0, 0, 0); int edge_count = 0;
		for (Mesh::VertexOHalfedgeIter voh_it = mesh.voh_iter(v_it); voh_it; ++voh_it)
		{
			ave_edge_len += edge_len[mesh.edge_handle(voh_it).idx()]; edge_count += 1;
			Mesh::FaceHandle fh = mesh.face_handle(voh_it);
			if (fh != Mesh::InvalidFaceHandle)
			{
				vn += face_n[fh.idx()];
			}
		}
		vn.normalize();
		v_v[v_it->idx()] = mesh.point(v_it) - vn*ave_edge_len / edge_count;
	}
	std::vector<std::vector<int>> all_tets(nf*6);
	for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		std::vector<int> one_face; 
		for (Mesh::FaceVertexIter fv_it = mesh.fv_iter(f_it); fv_it; ++fv_it)
		{
			one_face.push_back( fv_it->idx() );
		}
		int f_id = f_it->idx();
		std::vector<int> one_tet(4);
		one_tet[0] = one_face[0]; one_tet[1] = one_face[2]; one_tet[2] = one_face[1]; one_tet[3] = one_face[0] + nv;
		all_tets[6 * f_id + 0] = one_tet;
		one_tet[0] = one_face[0]; one_tet[1] = one_face[2]; one_tet[2] = one_face[1]; one_tet[3] = one_face[1] + nv;
		all_tets[6 * f_id + 1] = one_tet;
		one_tet[0] = one_face[0]; one_tet[1] = one_face[2]; one_tet[2] = one_face[1]; one_tet[3] = one_face[2] + nv;
		all_tets[6 * f_id + 2] = one_tet;
		one_tet[0] = one_face[0] + nv; one_tet[1] = one_face[1] + nv; one_tet[2] = one_face[2] + nv; one_tet[3] = one_face[0];
		all_tets[6 * f_id + 3] = one_tet;
		one_tet[0] = one_face[0] + nv; one_tet[1] = one_face[1] + nv; one_tet[2] = one_face[2] + nv; one_tet[3] = one_face[1];
		all_tets[6 * f_id + 4] = one_tet;
		one_tet[0] = one_face[0] + nv; one_tet[1] = one_face[1] + nv; one_tet[2] = one_face[2] + nv; one_tet[3] = one_face[2];
		all_tets[6 * f_id + 5] = one_tet;
	}

	FILE* f_tet = fopen("tet.tet", "w");
	fprintf(f_tet, "Vertices %d", nv+nv);
	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		OpenMesh::Vec3d p = mesh.point(v_it);
		fprintf(f_tet, "\n%.19f %.19f %.19f", p[0], p[1], p[2]);
	}
	for (int i = 0; i < nv;++i)
	{
		fprintf(f_tet, "\n%.19f %.19f %.19f", v_v[i][0], v_v[i][1], v_v[i][2]);
	}
	fprintf(f_tet, "\nTetrahedra %d", nf*6);
	for (int i = 0; i < all_tets.size(); ++i)
	{
		fprintf(f_tet, "\n%d %d %d %d", all_tets[i][0] + 1, all_tets[i][1] + 1, all_tets[i][2] + 1, all_tets[i][3] + 1);
	}
	fclose(f_tet);*/
}


bool MeshViewerWidget::saveMesh(const char* filename)
{
	return OpenMesh::IO::write_mesh(mesh, filename);
}
bool MeshViewerWidget::saveScreen(const char* filePath)
{
	QImage image = grabFrameBuffer();
	image.save(filePath);
	return true;
}

void MeshViewerWidget::checkMeshMode()
{
	Mesh::FaceIter fIt = mesh.faces_begin();
	Mesh::FaceIter fEnd = mesh.faces_end();
	Mesh::FaceEdgeIter fe_it;
	int count = 1;
	int meshType[3] = {0};
	for(fIt; fIt != fEnd; ++fIt)
	{
		fe_it = mesh.fe_iter(fIt);
		while(--fe_it)
		{
			++count;
		}
		if(count == 4)
		{
			meshType[1]++;
		}
		else if(count == 3)
		{
			meshType[0]++;
		}
		else
		{
			meshType[2]++;
		}
		count = 1;
	}
	int faceNum = mesh.n_faces();
	if(meshType[0] == faceNum)//triangle
	{
		setMeshMode(TRIANGLE);
	}
	else if(meshType[1] == faceNum)//no 
	{
		setMeshMode(QUAD);
	}
	else
	{
		setMeshMode(N_MESH_MODES);
	}
}

void MeshViewerWidget::updateIndices()
{
	//Mesh::ConstFaceIter f_it(mesh.faces_sbegin()), 
	//	f_end(mesh.faces_end());
	//Mesh::ConstFaceVertexIter  fv_it;

	////Indices.clear(); VIndices.clear();
	//int PolygonClass = 0 ;
	//switch(meshMode())
	//{
	//case TRIANGLE:
	//	PolygonClass = 3;
	//	break;
	//case QUAD:
	//	PolygonClass = 4;
	//	break;
	//default:
	//	return;
	//}
	//
	//Indices.reserve(mesh.n_faces()*PolygonClass);
	//for (; f_it!=f_end; ++f_it)
	//{
	//	fv_it = mesh.cfv_iter(f_it); 
	//	for( fv_it; fv_it; ++fv_it)
	//	{
	//		Indices.push_back((fv_it).handle().idx());
	//	}
	//}

	//VIndices.resize(mesh.n_vertices());
	//for(int i=0;i<mesh.n_vertices();++i)
	//{
	//	VIndices[i] = i;
	//}
}

void MeshViewerWidget::draw_scene(int drawmode)
{
	QFont Text_Font("Courier", 12);
	glViewport ( 0,0, width(),height());
	glMatrixMode( GL_PROJECTION );
	glLoadMatrixd( &ProjectionMatrix[0] );
	glMatrixMode( GL_MODELVIEW );
	glLoadMatrixd( &ModelViewMatrix[0] );

	/*double r = (bbMin-bbMax).norm()*0.5;
	OpenMesh::Vec3d c = (bbMin+ bbMax)*0.5;
	OpenMesh::Vec3d x(1.0*r,0,0); OpenMesh::Vec3d y(0,1.0*r,0); OpenMesh::Vec3d z(0,0,1.0*r);
	OpenMesh::Vec3d temp;
	glDisable(GL_LIGHTING);
	glLineWidth(2.0);
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3dv( c.data() );
	temp = (c + x);
	glVertex3dv( temp.data() );
	glEnd();
	renderText(temp[0] + 0.01*r, temp[1], temp[2] , "X", Text_Font );

	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINES);
	glVertex3dv( c.data() );
	temp = (c + y);
	glVertex3dv( temp.data() );
	glEnd();
	renderText(temp[0], temp[1] + 0.01*r, temp[2] , "Y", Text_Font );

	glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_LINES);
	glVertex3dv( c.data() );
	temp = (c + z);
	glVertex3dv( temp.data() );
	glEnd();
	renderText(temp[0], temp[1], temp[2] + 0.01*r , "Z", Text_Font );*/

	if(draw_BBox_OK)
	{
		OpenMesh::Vec3d temp0 = bbMin;
		OpenMesh::Vec3d temp1;
		glLineWidth(2.0);
		glColor3f(1.0, 1.0, 0.0);
		glBegin(GL_LINES);
		temp1 = bbMin; temp1[0] = bbMax[0];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );
		temp1 = bbMin; temp1[1] = bbMax[1];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );
		temp1 = bbMin; temp1[2] = bbMax[2];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );

		temp0 = bbMin; temp0[0] = bbMax[0];
		temp1 = bbMax; temp1[1] = bbMin[1];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );

		temp0 = bbMin; temp0[0] = bbMax[0];
		temp1 = bbMax; temp1[2] = bbMin[2];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );

		temp0 = bbMin; temp0[1] = bbMax[1];
		temp1 = bbMax; temp1[2] = bbMin[2];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );

		temp0 = bbMin; temp0[1] = bbMax[1];
		temp1 = bbMax; temp1[0] = bbMin[0];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );

		temp0 = bbMin; temp0[2] = bbMax[2];
		temp1 = bbMax; temp1[1] = bbMin[1];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );

		temp0 = bbMin; temp0[2] = bbMax[2];
		temp1 = bbMax; temp1[0] = bbMin[0];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );

		temp0 = bbMax;
		temp1 = bbMax; temp1[0] = bbMin[0];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );
		temp1 = bbMax; temp1[1] = bbMin[1];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );
		temp1 = bbMax; temp1[2] = bbMin[2];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );
		glEnd();
	}

	if(draw_mesh_boundary_ok)
	{
		glLineWidth(2.0);
		glColor3f(1.0, 0.5, 0.0);
		glBegin(GL_LINES);
		for(Mesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
		{
			if( mesh.is_boundary(e_it ) )
			{
				Mesh::HalfedgeHandle heh0 = mesh.halfedge_handle(e_it, 0);
				glVertex3dv( mesh.point(mesh.to_vertex_handle(heh0)).data() );
				glVertex3dv( mesh.point(mesh.from_vertex_handle(heh0)).data() );
			}
		}
		glEnd();

		/*FILE* f_bde = fopen("bde.de", "w");
		for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
		{
			if (mesh.is_boundary(v_it))
			{
				OpenMesh::Vec3d& np = mesh.point(v_it);
				fprintf(f_bde, "%d %20.19f %20.19f 0.0\n", v_it.handle().idx(), np[0], np[1]);
			}
		}
		fclose(f_bde);*/
	}

	draw_scene_mesh(drawmode);
}

void MeshViewerWidget::draw_scene_mesh(int drawmode)
{
	if(mesh.n_vertices() == 0) { return; }

	switch (drawmode)
	{
	case WIRE_FRAME:
		glDisable(GL_LIGHTING);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		draw_mesh_wireframe();
		//draw_meshpointset();
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		break;
	case HIDDEN_LINES:
		glDisable(GL_LIGHTING);
		draw_mesh_hidden_lines();
		break;
	case SOLID_FLAT:
		glEnable(GL_LIGHTING);
		glShadeModel(GL_FLAT);
		draw_mesh_solidflat();
		//draw_meshpointset();
		glDisable(GL_LIGHTING);
		break;
	case FLAT_POINTS:
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.5f,2.0f);
		glEnable(GL_LIGHTING);
		glShadeModel(GL_FLAT);
		draw_mesh_solidflat();
		glDisable(GL_POLYGON_OFFSET_FILL);
		//draw_meshpointset();
		glDisable(GL_LIGHTING);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		draw_mesh_wireframe();
		//draw_feature();
		//draw_feature1();
		//draw_meshpointset();
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		break;
	case SOLID_SMOOTH:
		glEnable(GL_LIGHTING);
		glShadeModel(GL_SMOOTH);
		draw_mesh_solidsmooth();
		//draw_meshpointset();
		glDisable(GL_LIGHTING);
		break;
	case POINT_SET:
		glDisable(GL_LIGHTING);
		draw_mesh_pointset();
		break;
	case CHECKBOARD:
		draw_CrossField();

		glEnable(GL_LIGHTING);
		glShadeModel(GL_FLAT);
		draw_mesh_solidflat();
		//draw_meshpointset();
		glDisable(GL_LIGHTING);
		break;
	case DIAGONAL_MESH:
		//glEnable(GL_POLYGON_OFFSET_FILL);
		//glPolygonOffset(1.5f, 2.0f);
		//glEnable(GL_LIGHTING);
		//glShadeModel(GL_FLAT);
		//draw_AnisotropicMesh();
		//glDisable(GL_POLYGON_OFFSET_FILL);
		////draw_meshpointset();
		//glDisable(GL_LIGHTING);
		//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		//draw_mesh_wireframe();
		////draw_meshpointset();
		//glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		draw_RegularCrossField();

		glEnable(GL_LIGHTING);
		glShadeModel(GL_FLAT);
		draw_mesh_solidflat();
		//draw_meshpointset();
		glDisable(GL_LIGHTING);
		break;
	default:
		break;
	}
}

void MeshViewerWidget::draw_mesh_wireframe()
{
	glLineWidth(1);
	//glColor3f(0.753, 0.753, 0.753);
	//glColor3f(0.0, 0.0, 0.0);
	glColor3f(0.0, 0.0, 0.25);
	//if(meshMode() != TRIANGLE && meshMode() != QUAD)
	{
		Mesh::ConstFaceIter fIt(mesh.faces_begin()),
			fEnd(mesh.faces_end());
		Mesh::ConstFaceVertexIter fvIt;
		for (; fIt != fEnd; ++fIt)
		{
			//GL::glNormal(dualMesh.normal(f_it));
			fvIt = mesh.cfv_iter(fIt); 
			glBegin(GL_POLYGON);
			for( fvIt; fvIt; ++fvIt )
			{
				glVertex3dv( mesh.point(fvIt).data() );
			}
			glEnd();
		}
		//std::cout<<"OK";
		//glColor3d(1., 0., 0.);
		//glBegin(GL_POLYGON);
		//glVertex3dv(mesh.point(mesh.vertex_handle(19953)).data());
		//glVertex3dv(mesh.point(mesh.vertex_handle(19952)).data());
		//glVertex3dv(mesh.point(mesh.vertex_handle(37708)).data());
		//glVertex3dv(mesh.point(mesh.vertex_handle(37707)).data());
		//glEnd();
	}

	/*OpenMesh::Vec3d pos = mesh.point(mesh.vertex_handle(0));
	GLdouble  winX, winY, winZ;
	GLint     viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	gluProject(pos[0],pos[1],pos[2],&ModelViewMatrix[0][0],&ProjectionMatrix[0][0],viewport,&winX,&winY,&winZ);
	int x = (long)winX;
	int y = viewport[3]-(long)winY;
	QString str = "fxum";
	render_text(x,y,str);*/
	/*else
	{
		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(3, GL_DOUBLE, 0, mesh.points());
		switch(meshMode())
		{
		case TRIANGLE:
			glDrawElements(GL_TRIANGLES,Indices.size(),GL_UNSIGNED_INT,&Indices[0]);
			break;
		case QUAD:
			glDrawElements(GL_QUADS,Indices.size(),GL_UNSIGNED_INT,&Indices[0]);
			break;
		default:
			break;
		}

		glDisableClientState(GL_VERTEX_ARRAY);
	}*/

}

void MeshViewerWidget::draw_mesh_hidden_lines() const
{
	Mesh::ConstFaceIter f_it(mesh.faces_begin());
	Mesh::ConstFaceIter f_end(mesh.faces_end());
	Mesh::ConstFaceVertexIter fv_it;
	glLineWidth(2.0);
	glColor3f(0.0, 1.0, 1.0);

	glDrawBuffer(GL_NONE);
	glDepthRange(0.01, 1.0);
	switch(meshMode())
	{
	case TRIANGLE:
		glBegin(GL_TRIANGLES);
		for (; f_it!=f_end; ++f_it)
		{
			fv_it = mesh.cfv_iter(f_it); 
			for(fv_it;fv_it;++fv_it)
			{
				glVertex3dv(&mesh.point(fv_it)[0]);
			}

		}
		glEnd();
		break;
	case QUAD:
		glBegin(GL_QUADS);
		for (; f_it!=f_end; ++f_it)
		{
			fv_it = mesh.cfv_iter(f_it); 
			for(fv_it;fv_it;++fv_it)
			{
				glVertex3dv(&mesh.point(fv_it)[0]);
			}
		}
		glEnd();
		break;
	default:
		for (; f_it!=f_end; ++f_it)
		{
			fv_it = mesh.cfv_iter(f_it); 
			glBegin(GL_POLYGON);
			for(fv_it;fv_it;++fv_it)
			{
				glVertex3dv(&mesh.point(fv_it)[0]);
			}
			glEnd();
		}
		break;
	}

	glDrawBuffer(GL_BACK);
	glDepthRange(0.0, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glDepthFunc(GL_LEQUAL);

	switch(meshMode())
	{
	case TRIANGLE:
		glBegin(GL_TRIANGLES);
		for (f_it = mesh.faces_begin(); f_it!=f_end; ++f_it)
		{
			fv_it = mesh.cfv_iter(f_it); 
			for(fv_it;fv_it;++fv_it)
			{
				glVertex3dv(&mesh.point(fv_it)[0]);
			}
		}
		glEnd();
		break;
	case QUAD:
		glBegin(GL_QUADS);
		for (f_it = mesh.faces_begin(); f_it!=f_end; ++f_it)
		{
			fv_it = mesh.cfv_iter(f_it); 
			for(fv_it;fv_it;++fv_it)
			{
				glVertex3dv(&mesh.point(fv_it)[0]);
			}
		}
		glEnd();
		break;
	default:
		for (f_it = mesh.faces_begin(); f_it!=f_end; ++f_it)
		{
			fv_it = mesh.cfv_iter(f_it); 
			glBegin(GL_POLYGON);
			for(fv_it;fv_it;++fv_it)
			{
				glVertex3dv(&mesh.point(fv_it)[0]);
			}
			glEnd();
		}
		break;
	}

	glDepthFunc(GL_LESS);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void MeshViewerWidget::draw_mesh_solidflat() const
{
	Mesh::ConstFaceIter fIt(mesh.faces_begin()),
		fEnd(mesh.faces_end());
	Mesh::ConstFaceVertexIter fvIt;

	GLfloat mat_a[] = { 0.7f, 0.7f, 0.7f, 1.0f };
	GLfloat mat_d[] = { 0.88f, 0.84f, 0.76f, 1.0f };
	GLfloat mat_s[] = { 0.4f, 0.4f, 0.4f, 1.0f };
	GLfloat shine[] = { 120.0f };

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_a);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_d);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_s);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shine);

	switch(meshMode())
	{
	case TRIANGLE:
		glBegin(GL_TRIANGLES);
		for (fIt; fIt != fEnd; ++fIt) 
		{
			glNormal3dv(mesh.normal(fIt).data());
			fvIt = mesh.cfv_iter(fIt.handle());
			for(fvIt; fvIt; ++fvIt)
			{
				glVertex3dv(mesh.point(fvIt).data());
			}
		}
		glEnd();
		break;
	case QUAD:
		glBegin(GL_QUADS);
		for (;fIt != fEnd; ++fIt) 
		{
			glNormal3dv(&mesh.normal(fIt)[0]);
			fvIt = mesh.cfv_iter(fIt.handle());
			for(fvIt; fvIt; ++fvIt)
			{
				glVertex3dv(&mesh.point(fvIt)[0]);
			}
		}
		glEnd();
		break;
	default:
		for (; fIt != fEnd; ++fIt) 
		{
			glBegin(GL_POLYGON);
			glNormal3dv(&mesh.normal(fIt)[0]);
			fvIt = mesh.cfv_iter(fIt.handle());
			for(fvIt; fvIt; ++fvIt)
			{
				glVertex3dv(&mesh.point(fvIt)[0]);
			}
			glEnd();
		}
		break;
	}
}

void MeshViewerWidget::draw_mesh_solidsmooth() const 
{
	bool drawOK = false;
	glLoadName(mesh.n_vertices());

	//glEnableClientState(GL_VERTEX_ARRAY);
	//glVertexPointer(3, GL_DOUBLE, 0, mesh.points());
	//glEnableClientState(GL_NORMAL_ARRAY);
	//glNormalPointer(GL_DOUBLE, 0, mesh.vertex_normals());
	//switch(meshMode())
	//{
	//case TRIANGLE:
	//	glDrawElements(GL_TRIANGLES,Indices.size(),GL_UNSIGNED_INT,&Indices[0]);
	//	drawOK = true;
	//	break;
	//case QUAD:
	//	glDrawElements(GL_QUADS,Indices.size(),GL_UNSIGNED_INT,&Indices[0]);
	//	drawOK = true;
	//	break;
	//default:
	//	break;
	//}
	//glDisableClientState(GL_VERTEX_ARRAY);
	//glDisableClientState(GL_NORMAL_ARRAY);

	//if(drawOK) return;

	Mesh::ConstFaceIter fIt(mesh.faces_begin()),
		fEnd(mesh.faces_end());
	Mesh::ConstFaceVertexIter fvIt;

	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_DOUBLE, 0, mesh.points());
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_DOUBLE, 0, mesh.vertex_normals());

	for (; fIt != fEnd; ++fIt)
	{
		glBegin(GL_POLYGON);
		fvIt = mesh.cfv_iter(fIt.handle());
		for(fvIt; fvIt ; --fvIt)
		{
			glArrayElement(fvIt.handle().idx());
		}
		glEnd();
	}


	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
}

void MeshViewerWidget::draw_mesh_pointset() const 
{
	glPolygonMode(GL_FRONT_AND_BACK, GL_POINTS);

	glColor3f(0.0, 1.0, 1.0);
	glPointSize(10);
	Mesh::VertexIter v_it = mesh.vertices_begin();
	glBegin(GL_POINTS);
	for(v_it;v_it != mesh.vertices_end();++v_it)
	{
		glVertex3dv(mesh.point(v_it).data());
	}
	glEnd();
	/*glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_DOUBLE, 0, mesh.points());
	glDrawElements(GL_POINTS,mesh.n_vertices(),GL_UNSIGNED_INT,&VIndices[0]);
	glDisableClientState(GL_VERTEX_ARRAY);*/

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

}

//#include "..\Algorithm\dualLoop.h"
#include "..\Toolbox\stream\dprint.h"
void MeshViewerWidget::draw_CrossField()
{
	if (if_has_ppl)
	{
		if (cf && if_has_field) { delete cf; cf = nullptr; }
		cf = &(ppl->gd.cf);
		cf->initFieldInfo();
		crossfield = cf->getCrossField();
		double avgLen = 0.2 * calc_mesh_ave_edge_length(&mesh);
		for (auto tf : mesh.faces())
		{
			OpenMesh::Vec3d c = mesh.calc_centroid(tf);
			int i = tf.idx() * 4;
			Eigen::Vector3d vc(c[0], c[1], c[2]);
			Eigen::Vector3d temp = crossfield.col(i + 1);
			crossfield.col(i) = vc + crossfield.col(i) * avgLen;
			crossfield.col(i + 1) = vc + crossfield.col(i + 2) * avgLen;
			crossfield.col(i + 2) = vc + temp * avgLen;
			crossfield.col(i + 3) = vc + crossfield.col(i + 3) * avgLen;
		}
	}
	else if (!if_has_field)
	{
		if_has_field = true;
		if (cf) { delete cf; cf = nullptr; }
		cf = new HDP::crossField(&mesh);

		cf->initMeshInfo();
		cf->setBoundaryConstraint();
		cf->setField();
		cf->initFieldInfo();
		crossfield = cf->getCrossField();
		double avgLen = 0.2 * calc_mesh_ave_edge_length(&mesh);
		for (auto tf : mesh.faces())
		{
			OpenMesh::Vec3d c = mesh.calc_centroid(tf);
			int i = tf.idx() * 4;
			Eigen::Vector3d vc(c[0], c[1], c[2]);
			Eigen::Vector3d temp = crossfield.col(i + 1);
			crossfield.col(i) = vc + crossfield.col(i) * avgLen;
			crossfield.col(i + 1) = vc + crossfield.col(i + 2) * avgLen;
			crossfield.col(i + 2) = vc + temp * avgLen;
			crossfield.col(i + 3) = vc + crossfield.col(i + 3) * avgLen;
		}
	}
	glLineWidth(1);
	glColor3d(0.9, 0.1, 0.1);
	glBegin(GL_LINES);
	for (int i = 0; i < crossfield.cols(); i += 2)
	{
		glVertex3dv(crossfield.col(i).data());
		glVertex3dv(crossfield.col(i + 1).data());
	}
	glEnd();
	draw_scene_mesh(FLAT_POINTS);

	auto& sing = cf->getSingularity();
	glPointSize(8);
	glBegin(GL_POINTS);
	glColor3d(1.0, 0, 0);
	for (auto s : sing[0])
		glVertex3dv(mesh.point(mesh.vertex_handle(s)).data());
	glColor3d(0, 1.0, 0);
	for (auto s : sing[1])
		glVertex3dv(mesh.point(mesh.vertex_handle(s)).data());
	glColor3d(0, 0, 1.0);
	for (auto s : sing[2])
		glVertex3dv(mesh.point(mesh.vertex_handle(s)).data());
	glEnd();
}

void MeshViewerWidget::draw_RegularCrossField()
{
	double ff = std::sqrt(3);
	double data[60] = {
		3,0,0,
		2,1,0,
		1,2,0,
		0,3,0,
		2,0,0,
		1,1,0,
		0,2,0,
		1,0,0,
		0,1,0,
		0,0,0,

		2,0,1,
		1,1,1,
		0,2,1,
		1,0,1,
		0,1,1,
		0,0,1,

		1,0,2,
		0,1,2,
		0,0,2,

		0,0,3
	};
	std::vector<OpenMesh::Vec3d> pos; pos.reserve(20);
	for (int i = 0; i < 20; i++)
	{
		pos.emplace_back(data[3 * i], data[3 * i + 1], data[3 * i + 2]);
	}
	glLineWidth(10);
	dprint("f");
	glBegin(GL_LINES);
	//���߿�
	glColor3d(0, 0, 0);
	glVertex3dv(pos[0].data()); glVertex3dv(pos[1].data());
	glVertex3dv(pos[1].data()); glVertex3dv(pos[2].data());
	glVertex3dv(pos[2].data()); glVertex3dv(pos[3].data());
	glVertex3dv(pos[3].data()); glVertex3dv(pos[6].data());
	glVertex3dv(pos[6].data()); glVertex3dv(pos[8].data());
	glVertex3dv(pos[8].data()); glVertex3dv(pos[9].data());
	glVertex3dv(pos[9].data()); glVertex3dv(pos[7].data());
	glVertex3dv(pos[7].data()); glVertex3dv(pos[4].data());
	glVertex3dv(pos[4].data()); glVertex3dv(pos[0].data());
	glVertex3dv(pos[0].data()); glVertex3dv(pos[10].data());
	glVertex3dv(pos[10].data()); glVertex3dv(pos[16].data());
	glVertex3dv(pos[16].data()); glVertex3dv(pos[19].data());
	glVertex3dv(pos[3].data()); glVertex3dv(pos[12].data());
	glVertex3dv(pos[12].data()); glVertex3dv(pos[17].data());
	glVertex3dv(pos[17].data()); glVertex3dv(pos[19].data());
	glVertex3dv(pos[9].data()); glVertex3dv(pos[15].data());
	glVertex3dv(pos[15].data()); glVertex3dv(pos[18].data());
	glVertex3dv(pos[18].data()); glVertex3dv(pos[19].data());

	//�����ϱ�
	glColor3d(1, 0, 0);
	glVertex3dv(pos[1].data()); glVertex3dv(pos[4].data());
	glVertex3dv(pos[4].data()); glVertex3dv(pos[5].data());
	glVertex3dv(pos[5].data()); glVertex3dv(pos[1].data());
	glVertex3dv(pos[2].data()); glVertex3dv(pos[5].data());
	glVertex3dv(pos[5].data()); glVertex3dv(pos[6].data());
	glVertex3dv(pos[6].data()); glVertex3dv(pos[2].data());
	glVertex3dv(pos[5].data()); glVertex3dv(pos[7].data());
	glVertex3dv(pos[7].data()); glVertex3dv(pos[8].data());
	glVertex3dv(pos[8].data()); glVertex3dv(pos[5].data());
	glColor3d(0, 1, 0);
	glVertex3dv(pos[1].data()); glVertex3dv(pos[10].data());
	glVertex3dv(pos[10].data()); glVertex3dv(pos[11].data());
	glVertex3dv(pos[11].data()); glVertex3dv(pos[1].data());
	glVertex3dv(pos[2].data()); glVertex3dv(pos[11].data());
	glVertex3dv(pos[11].data()); glVertex3dv(pos[12].data());
	glVertex3dv(pos[12].data()); glVertex3dv(pos[2].data());
	glVertex3dv(pos[11].data()); glVertex3dv(pos[16].data());
	glVertex3dv(pos[16].data()); glVertex3dv(pos[17].data());
	glVertex3dv(pos[17].data()); glVertex3dv(pos[11].data());
	glColor3d(0, 0, 1);
	glVertex3dv(pos[4].data()); glVertex3dv(pos[10].data());
	glVertex3dv(pos[10].data()); glVertex3dv(pos[13].data());
	glVertex3dv(pos[13].data()); glVertex3dv(pos[4].data());
	glVertex3dv(pos[18].data()); glVertex3dv(pos[13].data());
	glVertex3dv(pos[13].data()); glVertex3dv(pos[16].data());
	glVertex3dv(pos[16].data()); glVertex3dv(pos[18].data());
	glVertex3dv(pos[13].data()); glVertex3dv(pos[15].data());
	glVertex3dv(pos[15].data()); glVertex3dv(pos[7].data());
	glVertex3dv(pos[7].data()); glVertex3dv(pos[13].data());
	glColor3d(1, 1, 0);
	glVertex3dv(pos[6].data()); glVertex3dv(pos[12].data());
	glVertex3dv(pos[12].data()); glVertex3dv(pos[14].data());
	glVertex3dv(pos[14].data()); glVertex3dv(pos[6].data());
	glVertex3dv(pos[18].data()); glVertex3dv(pos[14].data());
	glVertex3dv(pos[14].data()); glVertex3dv(pos[17].data());
	glVertex3dv(pos[17].data()); glVertex3dv(pos[18].data());
	glVertex3dv(pos[14].data()); glVertex3dv(pos[15].data());
	glVertex3dv(pos[15].data()); glVertex3dv(pos[8].data());
	glVertex3dv(pos[8].data()); glVertex3dv(pos[14].data());

	//�ڲ���
	glColor3d(0, 1, 1);
	glVertex3dv(pos[5].data()); glVertex3dv(pos[11].data());
	glVertex3dv(pos[5].data()); glVertex3dv(pos[14].data());
	glVertex3dv(pos[5].data()); glVertex3dv(pos[13].data());
	glVertex3dv(pos[11].data()); glVertex3dv(pos[13].data());
	glVertex3dv(pos[11].data()); glVertex3dv(pos[14].data());
	glVertex3dv(pos[13].data()); glVertex3dv(pos[14].data());

	glEnd();

	//if (!loop_gen_init)
	//{
	//	lg->InitializeGraphWeight();
	//}

	//static bool flag = true;
	//static std::vector<int> loop;
	//static Eigen::Matrix3Xd crossfield;
	//static double avglen = calc_mesh_ave_edge_length(&mesh) * 0.2;
	//static Eigen::Vector3d ev;
	//if (flag)
	//{
	//	LoopGen::LoopGen lg(mesh);
	//	lg.InitializeField();
	//	crossfield = lg.cf->getCrossField();
	//	auto position = lg.cf->getPosition();
	//	ev = position.col(20918) - position.col(29054);
	//	//dprint(crossfield.col(3077 * 4).dot(ev.normalized()), crossfield.col(48306 * 4).dot(ev.normalized()), crossfield.col(48306 * 4 + 1).dot(ev.normalized()));
	//	//dprint("fheu",crossfield.col(3077 * 4), crossfield.col(48306 * 4), crossfield.col(48306 * 4 + 1));
	//	for (auto tf : mesh.faces())
	//	{
	//		OpenMesh::Vec3d c = mesh.calc_centroid(tf);
	//		int i = tf.idx() * 4;
	//		Eigen::Vector3d vc(c[0], c[1], c[2]);
	//		Eigen::Vector3d temp = crossfield.col(i + 1);
	//		crossfield.col(i) = vc + crossfield.col(i) * avglen;
	//		crossfield.col(i + 1) = vc;// +crossfield.col(i + 2) * avglen;
	//		crossfield.col(i + 2) = vc + temp * avglen;
	//		crossfield.col(i + 3) = vc + crossfield.col(i + 3) * avglen;
	//	}
	//	//lg.cf->write_field();
	//	lg.InitializeGraphWeight();
	//	lg.FieldAligned_PlanarLoop(mesh.vertex_handle(9017), loop, 0);
	//	flag = false;
	//}
	////dprint((crossfield.col(3077 * 4)- crossfield.col(3077 * 4+1)).normalized().dot(ev.normalized()), (crossfield.col(48306 * 4)- crossfield.col(48306 * 4+1)).normalized().dot(ev.normalized()), crossfield.col(48306 * 4 + 1).dot(ev.normalized()));

	////dprint("iuek",(crossfield.col(3077 * 4) - crossfield.col(3077 * 4 + 1)).normalized(), (crossfield.col(48306 * 4) - crossfield.col(48306 * 4 + 1)).normalized()); glColor3d(0.1, 0.1, 0.9);
	//glLineWidth(10);
	//glBegin(GL_LINES);
	//for (int i = 0; i < loop.size() - 1; ++i)
	//{
	//	glVertex3dv(mesh.point(mesh.vertex_handle(loop[i])).data());
	//	glVertex3dv(mesh.point(mesh.vertex_handle(loop[i + 1])).data());
	//}
	//glEnd();

	//

	//glColor3d(0.1, 0.2, 0.8);
	//glPointSize(10);
	//glBegin(GL_POINTS);
	//glVertex3dv(mesh.point(mesh.vertex_handle(9017)).data());
	//glEnd();

	//glColor3d(0.1, 0.9, 0.1);
	//glBegin(GL_POLYGON);
	//for (auto tfv : mesh.fv_range(mesh.face_handle(mesh.voh_begin(mesh.vertex_handle(9017)).handle())))
	//{
	//	glVertex3dv(mesh.point(tfv).data());
	//}
	//glEnd();
}
