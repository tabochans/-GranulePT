#include"Grain.h"
#include"Scene.h"
#include "Eigen/Geometry"

using namespace PTUtility;

void iGrain::InitGrain(float VoxelWidth, const char* PhaseFileName){

	std::ifstream file;
	file.open(PhaseFileName);

	m_PPHE = new PhaseEstimater(100);

	if (!file){
		m_PPHE->EstimatePhaseFunction(this, m_Lamda_dt, m_Lamda_s[0], m_Beta);
		std::ofstream ofile;
		ofile.open(PhaseFileName);
		m_PPHE->WritePhaseFunction(ofile);
		ofile.close();
	}
	else{
		m_PPHE->LoadPhaseFunction(this, file);
		file.close();
	}

}


bool Grain_Sphere::IfIntersect(const Ray& r, const Vec3& RCenter, int RandomSeed)const{

	float t;
	return Intersect::Ray_Sphere(r, RCenter, GetGrainRadius(), t);

}


bool Grain_TranslucentSphere::NextRay(
	const Ray& r,
	const Vec3& Center,
	int RandomSeed,
	Ray& result,
	float& PDF,
	RandomMT& MT)const{

	float t;
	if (!Intersect::Ray_Sphere_In(r, Center, GetGrainRadius(), t)){
		return false;
	}

	Vec3 P = r.m_Org + (t)*r.m_Dir;
	Vec3 N = (P - Center).normalized();

	if (r.m_Dir.dot(N) > 0){
		//物体内部

		float FL = t;
		float uu = MT.genrand64_real3();
		float PSL = -(m_SSS)*log(uu + 0.00000001);

		if (FL > PSL){

			result.m_Org = r.m_Org + PSL * r.m_Dir;
			PDF = ISample::NextDirSphere(r, N, result.m_Dir, MT);
			PDF = 1.0f;

		}
		else{
			
			if (!m_Diffuse){
				//フレネル反射
				
				float nt = 1.33f;
				float nf = 1.0f;
				if (N.dot(r.m_Dir) > 0){
					nt = 1.0f;
					nf = 1.33f;
				}
				PDF = 1.0f / ISample::NextDirRefraction(r, N, nf, nt, result.m_Dir, MT);
				

				result.m_Org = P + result.m_Dir*0.0000001;
				result.m_Index = nt;

			}
			else{
				//拡散反射

				if (MT.genrand64_real1() < 0.85){
					ISample::NextDirHemSphere(r, N, result.m_Dir, MT);
				}
				else{
					ISample::NextDirHemSphere(r, -N, result.m_Dir, MT);
				}

				result.m_Org = P + result.m_Dir*0.0000001;

			}

		}




	}
	else{
		//物体外部


		if (!m_Diffuse){
			//フレネル反射
			
			float nt = 1.33f;
			float nf = 1.0f;
			if (N.dot(r.m_Dir) > 0){
				nt = 1.0f;
				nf = 1.33f;
			}
			PDF = 1.0f / ISample::NextDirRefraction(r, N, nf, nt, result.m_Dir, MT);
			

			result.m_Org = P + result.m_Dir*0.0000001;
			result.m_Index = nt;

		}
		else{
			//拡散反射

			if (MT.genrand64_real1() < 0.6){
				ISample::NextDirHemSphere(r, N, result.m_Dir, MT);
			}
			else{
				ISample::NextDirHemSphere(r, -N, result.m_Dir, MT);
			}

			result.m_Org = P + result.m_Dir*0.0000001;

		}

	}


	


	
	return true;
}

bool Grain_GlassSphere::NextRay(
	const Ray& r,
	const Vec3& Center,
	int RandomSeed,
	Ray& result,
	float& PDF,
	RandomMT& MT)const{

	float t;
	if(!Intersect::Ray_Sphere_In(r, Center, GetGrainRadius(), t)){
		return false;
	}

	Vec3 P = r.m_Org + (t)*r.m_Dir;
	Vec3 N = (P - Center).normalized();

	float nt = 1.33f;
	float nf = 1.0f;
	if (N.dot(r.m_Dir) > 0){
		nt = 1.0f;
		nf = 1.33f;
	}

	PDF = 1.0f / ISample::NextDirRefraction(r, N, nf, nt, result.m_Dir, MT);

	result.m_Org = P + result.m_Dir*0.00001;
	result.m_Index = nt;

	return true;
}


bool Grain_MetalSphere::NextRay(
	const Ray& r,
	const Vec3& Center,
	int RandomSeed,
	Ray& result,
	float& PDF,
	RandomMT& MT)const{

	float t = 0;
	if (!Intersect::Ray_Sphere(r, Center, GetGrainRadius(), t)){
		return false;
	}

	Vec3 P = r.m_Org + (t - 0.0000001)*r.m_Dir;
	Vec3 N = (P - Center).normalized();

	PDF = ISample::NextDirSpecular(r, N, result.m_Dir, MT);
	result.m_Org = P + N*0.000001;

	return true;

}


bool Grain_DiffuseSphere::NextRay(
	const Ray& r,
	const Vec3& Center,
	int RandomSeed,
	Ray& result,
	float& PDF,
	RandomMT& MT) const{

	float t = 0;
	if (!Intersect::Ray_Sphere(r, Center, GetGrainRadius(), t)){
		return false;
	}

	Vec3 P = r.m_Org + (t - 0.0000001)*r.m_Dir;
	Vec3 N = (P - Center).normalized();

	PDF = ISample::NextDirDiffuse(r, N, result.m_Dir, MT);
	result.m_Org = P + N*0.000001;

	return true;

}


bool Grain_FromMesh::NextRay(
	const Ray& r,
	const Vec3& Center,
	int RandomSeed,
	Ray& result,
	float& PDF,
	RandomMT& MT)const{

	switch (m_mat){

	case MATERIAL_TYPE::Diffuse:
		return PNextRay_DIFFUSE(r, Center, RandomSeed, result, PDF, MT);
		
	case MATERIAL_TYPE::METAL:
		return PNextRay_METAL(r, Center, RandomSeed, result, PDF, MT);
		
	case MATERIAL_TYPE::GLASS:
		return PNextRay_GLASS(r, Center, RandomSeed, result, PDF, MT);

	}

	return false;
}



bool Grain_FromMesh::PNextRay_DIFFUSE(
	const Ray& r,
	const Vec3& Center,
	int RandomSeed,
	Ray& result,
	float& PDF,
	RandomMT& MT)const{

	PDF = 1.0f;

	//らんだむろーてーと！！
	RandomMT mm(RandomSeed);
	float th = mm.genrand64_real1() * M_PI;
	float ph = mm.genrand64_real1() * 2.0f * M_PI;
	float Rot = mm.genrand64_real1() * 2.0f * M_PI;

	Eigen::Vector3f jI(sinf(th)*cosf(ph), sinf(th)*sinf(ph), cosf(th));
	Eigen::Quaternionf q(Eigen::AngleAxisf(Rot, jI));
	Eigen::Quaternionf iq = q.inverse();

	Ray HR = r;

	HR.m_Org -= Center;
	HR.m_Org /= (GetR() * GetGrainRadius());
	HR.m_Org -= GetGrainCenter();

	Eigen::Vector3f tt = q * Eigen::Vector3f(HR.m_Org[0], HR.m_Org[1], HR.m_Org[2]);
	HR.m_Org = Vec3(tt[0], tt[1], tt[2]);
	
	tt = q * Eigen::Vector3f(HR.m_Dir[0], HR.m_Dir[1], HR.m_Dir[2]);
	HR.m_Dir = Vec3(tt[0], tt[1], tt[2]);


	std::vector<AbsG::RF> ffa;
	m_scene->GetFacetFromRay(
		Ray(
		Vec3(HR.m_Org[0], HR.m_Org[1], HR.m_Org[2]), 
		Vec3(HR.m_Dir[0], HR.m_Dir[1], HR.m_Dir[2])), ffa);

	if (ffa.size() == 0){
		return false;
	}
	
	Vec3 Normal;
	int fi = -1;
	for (int i = 0; i < ffa.size(); i++){
		Normal = 
			ffa[0]._u * (ffa[0]._pF->GetNormal(1) - ffa[0]._pF->GetNormal(0)) +
			ffa[0]._v * (ffa[0]._pF->GetNormal(2) - ffa[0]._pF->GetNormal(0)) +
			ffa[0]._pF->GetNormal(0);
		Normal.normalize();
		if (Normal.dot(HR.m_Dir) < 0){
			fi = i;
			break;
		}
	}
	if (fi != -1){

		Vec3 Pos = HR.m_Org + ffa[fi]._t * HR.m_Dir;

		ISample::NextDirDiffuse(r, Normal, HR.m_Dir, MT);
		HR.m_Org = Pos;

		tt = iq * Eigen::Vector3f(HR.m_Org[0], HR.m_Org[1], HR.m_Org[2]);
		HR.m_Org = Vec3(tt[0], tt[1], tt[2]);

		tt = iq * Eigen::Vector3f(HR.m_Dir[0], HR.m_Dir[1], HR.m_Dir[2]);
		HR.m_Dir = Vec3(tt[0], tt[1], tt[2]);

		HR.m_Org += GetGrainCenter();
		HR.m_Org *= (GetR() * GetGrainRadius());
		HR.m_Org += Center;

		result.m_Dir = HR.m_Dir;
		result.m_Org = HR.m_Org + result.m_Dir * 0.000001;

	}
	else{
		return false;
	}
	



	return true;
}
bool Grain_FromMesh::PNextRay_METAL(
	const Ray& r,
	const Vec3& Center,
	int RandomSeed,
	Ray& result,
	float& PDF,
	RandomMT& MT)const{

	PDF = 1.0f;

	//らんだむろーてーと！！
	RandomMT mm(RandomSeed);
	float th = mm.genrand64_real1() * M_PI;
	float ph = mm.genrand64_real1() * 2.0f * M_PI;
	float Rot = mm.genrand64_real1() * 2.0f * M_PI;

	Eigen::Vector3f jI(sinf(th)*cosf(ph), sinf(th)*sinf(ph), cosf(th));
	Eigen::Quaternionf q(Eigen::AngleAxisf(Rot, jI));
	Eigen::Quaternionf iq = q.inverse();

	Ray HR = r;
	HR.m_Org -= Center;
	HR.m_Org /= (GetR() * GetGrainRadius());
	HR.m_Org -= GetGrainCenter();

	Eigen::Vector3f tt = q * Eigen::Vector3f(HR.m_Org[0], HR.m_Org[1], HR.m_Org[2]);
	HR.m_Org = Vec3(tt[0], tt[1], tt[2]);

	tt = q * Eigen::Vector3f(HR.m_Dir[0], HR.m_Dir[1], HR.m_Dir[2]);
	HR.m_Dir = Vec3(tt[0], tt[1], tt[2]);



	std::vector<AbsG::RF> ffa;
	m_scene->GetFacetFromRay(HR, ffa);

	if (ffa.size() == 0){
		return false;
	}

	Vec3 Normal;
	int fi = -1;
	for (int i = 0; i < ffa.size(); i++){
		Normal =
			ffa[0]._u * (ffa[0]._pF->GetNormal(1) - ffa[0]._pF->GetNormal(0)) +
			ffa[0]._v * (ffa[0]._pF->GetNormal(2) - ffa[0]._pF->GetNormal(0)) +
			ffa[0]._pF->GetNormal(0);
		Normal.normalize();
		if (Normal.dot(HR.m_Dir) < 0){
			fi = i;
			break;
		}
	}
	if (fi != -1){

		Vec3 Pos = HR.m_Org + ffa[fi]._t * HR.m_Dir;
		Vec3 re;
		ISample::NextDirSpecular(HR, Normal, re, MT);
		HR.m_Org = Pos;
		HR.m_Dir = re;

		tt = iq * Eigen::Vector3f(HR.m_Org[0], HR.m_Org[1], HR.m_Org[2]);
		HR.m_Org = Vec3(tt[0], tt[1], tt[2]);

		tt = iq * Eigen::Vector3f(HR.m_Dir[0], HR.m_Dir[1], HR.m_Dir[2]);
		HR.m_Dir = Vec3(tt[0], tt[1], tt[2]);

		HR.m_Org += GetGrainCenter();
		HR.m_Org *= (GetR() * GetGrainRadius());
		HR.m_Org += Center;

		result.m_Dir = HR.m_Dir;
		result.m_Org = HR.m_Org + result.m_Dir * 0.000001;

	}
	else{
		return false;
	}
	
	



	return true;
	
}
bool Grain_FromMesh::PNextRay_GLASS(
	const Ray& r,
	const Vec3& Center,
	int RandomSeed,
	Ray& result,
	float& PDF,
	RandomMT& MT)const{

	//らんだむろーてーと！！
	RandomMT mm(RandomSeed);
	float th = mm.genrand64_real1() * M_PI;
	float ph = mm.genrand64_real1() * 2.0f * M_PI;
	float Rot = mm.genrand64_real1() * 2.0f * M_PI;

	Eigen::Vector3f jI(sinf(th)*cosf(ph), sinf(th)*sinf(ph), cosf(th));
	Eigen::Quaternionf q(Eigen::AngleAxisf(Rot, jI));
	Eigen::Quaternionf iq = q.inverse();

	Ray HR = r;
	HR.m_Org -= Center;
	HR.m_Org /= (GetR() * GetGrainRadius());
	HR.m_Org -= GetGrainCenter();

	Eigen::Vector3f tt = q * Eigen::Vector3f(HR.m_Org[0], HR.m_Org[1], HR.m_Org[2]);
	HR.m_Org = Vec3(tt[0], tt[1], tt[2]);

	tt = q * Eigen::Vector3f(HR.m_Dir[0], HR.m_Dir[1], HR.m_Dir[2]);
	HR.m_Dir = Vec3(tt[0], tt[1], tt[2]);


	std::vector<AbsG::RF> ffa;
	m_scene->GetFacetFromRay(HR, ffa);

	if (ffa.size() == 0){
		return false;
	}

	Vec3 Normal;
	Normal =
		ffa[0]._u * (ffa[0]._pF->GetNormal(1) - ffa[0]._pF->GetNormal(0)) +
		ffa[0]._v * (ffa[0]._pF->GetNormal(2) - ffa[0]._pF->GetNormal(0)) +
		ffa[0]._pF->GetNormal(0);
	Normal.normalize();

	float N_from;
	float N_to;
	if (Normal.dot(HR.m_Dir) < 0){
		N_from = 1.0f;
		N_to = 1.2;
	}
	else{
		N_from = 1.2f;
		N_to = 1.0;
	}
	
	Vec3 Pos = HR.m_Org + ffa[0]._t * HR.m_Dir;

	Vec3 re;
	PDF = ISample::NextDirRefraction(HR, Normal, N_from, N_to, re, MT);
	HR.m_Org = Pos;
	HR.m_Dir = re;
	

	tt = iq * Eigen::Vector3f(HR.m_Org[0], HR.m_Org[1], HR.m_Org[2]);
	HR.m_Org = Vec3(tt[0], tt[1], tt[2]);

	tt = iq * Eigen::Vector3f(HR.m_Dir[0], HR.m_Dir[1], HR.m_Dir[2]);
	HR.m_Dir = Vec3(tt[0], tt[1], tt[2]);

	HR.m_Org += GetGrainCenter();
	HR.m_Org *= (GetR() * GetGrainRadius());
	HR.m_Org += Center;

	result.m_Dir = HR.m_Dir;
	result.m_Org = HR.m_Org + result.m_Dir * 0.000001;



	return true;

}



Grain_FromMesh::Grain_FromMesh(MATERIAL_TYPE type, const Vec3& Color, float Radi, const char* FileName) : m_scene(new Scene(new OBVH())), iGrain(Radi, FileName){

	m_scene->CreateScene(FileName);
	m_Alpha_s = m_Color = Color;
	m_Lamda_s = Vec3(0.01, 0.01, 0.01);
	m_Lamda_dt = 0.0f;
	m_Sigma_t = Vec3(1, 1, 1);
	m_Beta = 0.99f;
	m_PF = Grain_FromMesh::Phase;
	m_CenterPos = m_scene->GetSceneCenter();
	m_mat = type;

}
Grain_FromMesh::~Grain_FromMesh(){

	delete m_scene;
}

bool Grain_FromMesh::IfIntersect(const Ray& r, const Vec3& RCenter, int RandomSeed)const{

	float t;
	return Intersect::Ray_Sphere(r, RCenter, GetR() * 0.99, t);

}

Vec3 Grain_FromMesh::GetGrainCenter()const{
	return m_scene->GetSceneCenter();
}
float Grain_FromMesh::GetGrainRadius()const{
	float ww = m_scene->GetSceneWidth();
	return ww * 0.495;
}


void Grain_FromMesh::InitGrain(float VoxelWidth, const char* PhaseFileName){

	std::ifstream file;
	file.open(PhaseFileName);

	m_PPHE = new PhaseEstimater(50);

	if (!file){
		m_PPHE->EstimatePhaseFunction2(this, m_Lamda_dt, m_Lamda_s[0], m_Beta);
		std::ofstream ofile;
		ofile.open(PhaseFileName);
		m_PPHE->WritePhaseFunction(ofile);
		ofile.close();
	}
	else{
		m_PPHE->LoadPhaseFunction(this, file);
		file.close();
	}

}
















