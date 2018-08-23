#include "MLScene.h"

#include <fstream>	
#include <iomanip>

#include "MLBody.h"
#include "MLCommon.h"
#include "MLWorld.h"
#include "MLAssembler.h"

MLScene::MLScene():
t(0.0),
h(1e-2),
grav(0.0, 0.0, 0.0)
{

}

MLScene::~MLScene()
{
	
}
