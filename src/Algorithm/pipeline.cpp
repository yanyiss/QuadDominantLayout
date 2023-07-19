#include "pipeline.h"
namespace HDP
{
	void pipeline::run()
	{
		deformer dfm(*mesh);
		dfm.update_field();
		dfm.run();
		gd.set_crossfield(*dfm.cf);
	}
	
}