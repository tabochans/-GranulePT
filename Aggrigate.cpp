#include"Aggrigate.h"


using namespace PTUtility;


int Aggrigate::GetSphereFromRay(
	const Vec3& VoxelCenter,
	const Ray& ray,
	std::vector<OcTree::SC>& result){

	result.clear();

	//for (int i = 0; i < m_Pack._Pos.size(); i++){
	//	float t;
	//	if(Intersect::Ray_Sphere(
	//		ray, VoxelCenter + (m_VoxelGrid->GetVoxelWidth()*1.001 * (m_Pack._Pos[i] - 0.5 * Vec3::Ones())), 
	//		m_Pack._R * m_VoxelGrid->GetVoxelWidth(), t)){

	//		SC sc;
	//		sc._Pos = VoxelCenter + (m_VoxelGrid->GetVoxelWidth()*1.001 * (m_Pack._Pos[i] - 0.5 * Vec3::Ones()));
	//		sc._t = t;
	//		sc._SID = i;

	//		result.push_back(sc);
	//	}
	//}

	m_Otree.GetSphereFromRay(Ray(ray.m_Org - VoxelCenter, ray.m_Dir), result);
	for (int i = 0; i < result.size(); i++){
		result[i]._Pos += VoxelCenter;
	}

	if (result.size() > 0){
		return result.size();
	}
	else{
		return 0;
	}

}


bool Aggrigate::LoadPack(const char* FileName){

	std::ifstream file;
	file.open(FileName);

	while (!file.eof()){
		std::string Line;
		float x, y, z;
		std::string ID("");

		std::getline(file, Line);
		std::istringstream LS(Line);
		LS >> ID;
		if (ID == "v"){
			LS >> x;
			LS >> y;
			LS >> z;
			m_Pack._Pos.push_back(Vec3(x, y, z));
		}

	}

	file.close();

	//ちょっとなんとかなんないのかね？
	m_Pack._R = 0.08617344288463569 * 0.5;
	//m_Pack._R = 0.18;
	//m_Pack._R = 0.45;

	return true;
}


Aggrigate::Aggrigate(float AggrigateWidth) : m_Traito(0.0f){

}

Aggrigate::~Aggrigate(){
	delete m_VoxelGrid;
}

bool Aggrigate::InitAggrigate(const char* FileName){

	
	LoadPack("PACK\\pp2.obj");

	m_VoxelGrid = new VoxelGrid(32);

	//ファイルから
	m_VoxelGrid->CreateVoxelGrid(FileName);

	InitCore();


	return true;
}


bool Aggrigate::InitAggrigate(const Scene& scene, int NVoxel){


	LoadPack("PACK\\pp2.obj");

	m_VoxelGrid = new VoxelGrid(5);

	//ボクセライズ
	m_VoxelGrid->CreateVoxelGridFromMesh(scene, NVoxel);
	

	std::ofstream file("VOXEL\\Vdt.obj");
	m_VoxelGrid->WriteVoxelData(file);
	file.close();

	InitCore();

	return true;
}


int Aggrigate::Intersect(const Ray& r, std::vector<const VoxelGrid::VPP>& result){

	return m_VoxelGrid->Intersect(r, result);
}


Vec3 Aggrigate::GetCenter(const VoxelGrid::VoxelData* v){
	return m_VoxelGrid->GetCenterPos(v);
}


bool Aggrigate::NextRay(
	int GrainID,
	const Ray& r,
	const Vec3& Center,
	int RandomSeed,
	Ray& result,
	float& PDF,
	RandomMT& MT){

	bool re = GetGrain(GrainID)->NextRay(r, Center, RandomSeed, result, PDF, MT);
	result.m_Dir.normalize();

	return re;
}