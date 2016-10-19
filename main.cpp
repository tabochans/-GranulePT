#include"gHeader.h"
#include"Facet.h"
#include"OBVH.h"
#include"Scene.h"
#include"Pathtracing.h"
#include"GPathTracing.h"
#include"IBL.h"
#include<ctime>
#include<vector>
#include<fstream>

using namespace PTUtility;

int main(){

	time_t start = time(NULL);


	//OBVH構造でのシーン作成
	Scene scene(new OBVH());

	//OBJ読み込み
	scene.CreateScene("OBJ\\ES.obj");
	//scene.CreateScene("OBJ\\ES.obj");

	////いちおーBVHをOBJとして書き出しておく
	{
		//std::ofstream file;
		//file.open("BVH.obj");
		//scene.WriteBVH(file);
		//file.close();
	}


	////カメラの設定
	Vec3 cp(-2.3, 0.75, 2.1);
	cp *= 0.90;
	float fc = 0.9f;
	Camera cm(fc * Vec3(-0, 0.5f, 3.0f), Vec3(0.01f, -0.44, -1.71f), Vec3(0, 1, 0), 1, 1, 1.0f);
	//Camera cm(fc * Vec3(-0, 1.5f, 3.0f), Vec3(0.01f, -0.94, -1.71f), Vec3(0, 1, 0), 1, 1, 1.0f);
	
	
	//Camera cm(Vec3(-0.07, -0.2, 2.5), Vec3(0.1f, -0.0, -1.71f), Vec3(0, 1, 0), 1, 1, 1.0f);
	//Camera cm(cp, -cp, Vec3(0, 1, 0), 1980, 1080, 1.0f);

	//Camera cm(Vec3(0.1, -0.0, 2.5), Vec3(-0.21f, -0.0, -1.71f), Vec3(0, 1, 0), 1, 1, 1.0f);


	//Image based lighting
	//IBLTexture ibl("IBL\\IBLS.png");
	IBLTexture ibl("IBL\\IBLI3.jpg");
	//IBLDLight ibl(Vec3(1, 0, 0), 6.0f, Vec3::Ones());

	const int sc = 1;
	//Pathtracing pt(512 * sc, 512 * sc);
	GPathTrace pt(512, 512);
	//pt.RenderSceneT(scene, cm, &ibl, 1, 1, time(NULL) - start);
	pt.RenderScene(scene, cm, &ibl, 1, 2);
	
	
	pt.WriteImage("00VV.bmp");





	return 0;
}