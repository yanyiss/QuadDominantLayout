#pragma once
#include "deformer.h"
namespace HDP
{
	class gui_data
	{
	public:
		gui_data() {}
		~gui_data() {}
	public:
		crossField cf;
		void set_crossfield(crossField &field) { cf = field; }
	};

	class pipeline
	{
	public:
		pipeline(Mesh& m) { mesh = &m; initMeshStatusAndNormal(*mesh); }
		~pipeline() {}

	public:
		Mesh* mesh;
		gui_data gd;

		void run();
	};

}

