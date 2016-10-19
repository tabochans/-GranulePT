#include"Otree.h"


using namespace PTUtility;

int OcTree::CreateTree(int PID){
	
	if (m_NodeArray[PID].IfLeaf){
		return 0;
	}

	float max[3][8];
	float min[3][8];

	for (int i = 0; i < 8; i++){
		for (int j = 0; j < 3; j++){
			max[j][i] = m_NodeArray[m_NodeArray[PID]._Children[i]]._AABBs.GetMaxPos()[j];
			min[j][i] = m_NodeArray[m_NodeArray[PID]._Children[i]]._AABBs.GetMinPos()[j];
		}
	}

	m_NodeArray[PID].CAABB[0][0] = _mm256_loadu_ps(min[0]);
	m_NodeArray[PID].CAABB[0][1] = _mm256_loadu_ps(min[1]);
	m_NodeArray[PID].CAABB[0][2] = _mm256_loadu_ps(min[2]);
	m_NodeArray[PID].CAABB[1][0] = _mm256_loadu_ps(max[0]);
	m_NodeArray[PID].CAABB[1][1] = _mm256_loadu_ps(max[1]);
	m_NodeArray[PID].CAABB[1][2] = _mm256_loadu_ps(max[2]);

	for (int i = 0; i < 8; i++){
		CreateTree(m_NodeArray[PID]._Children[i]);
	}

	return 0;
}

int OcTree::PGetSphereFromRay(int ID, const Ray& ray, std::vector<SC>& result){

	float u, v, t;
	u = v = t = 0.0f;
	if (m_NodeArray[ID].IfLeaf){
		//リーフノード

		int nc = m_NodeArray[ID]._Sphere.size();

		for (int i = 0; i < nc; i++){
			if (Intersect::Ray_Sphere(ray, m_SphereArray[m_NodeArray[ID]._Sphere[i]], m_SR, t)){
				SC ss;
				ss._Pos = m_SphereArray[m_NodeArray[ID]._Sphere[i]];
				ss._SID = m_NodeArray[ID]._Sphere[i];
				ss._t = t;
				result.push_back(ss);
			}
		}
		return result.size();
	}
	else{
		//下を探索

		__m256 org[3];
		__m256 idir[3];
		int sign[3];
		int mask = 0;

		for (int i = 0; i < 3; i++){

			sign[i] = ray.m_Dir[i] < 0;
			org[i] = _mm256_set1_ps(ray.m_Org[i]);

			if (abs(ray.m_Dir[i]) < 1.0f / 100000.0f){
				idir[i] = _mm256_set1_ps(100000.0f - 200000.0f * sign[i]);
			}
			else{
				idir[i] = _mm256_set1_ps(1.0f / ray.m_Dir[i]);
			}

		}


		if (Intersect::IsHit(m_NodeArray[ID].CAABB, org, idir, sign, mask)){
			for (int i = 0; i < 8; i++){
				int fg = mask & (1 << (i));
				if ((fg) > 0){
					int iid = m_NodeArray[ID]._Children[i];
					if (iid > 0){ PGetSphereFromRay(iid, ray, result); }
				}
			}
		}

		return result.size();

	}

}

int OcTree::GetSphereFromRay(const Ray& ray, std::vector<SC>& result){

	result.clear();
	int n = PGetSphereFromRay(1, ray, result);
	
	if (n > 0) { std::sort(result.begin(), result.end()); }
	
	
	return n;
}

bool OcTree::PCreateOctree(int ParentID, int depth){

	if (depth >= m_Depth){
		m_NodeArray[ParentID].IfLeaf = true;
		return false;
	}
	if (m_NodeArray[ParentID]._Sphere.size() == 0){
		m_NodeArray[ParentID].IfLeaf = true;
		return false;
	}


	//子ノード
	Node child[2][2][2];

	//境界線
	Vec3 max = m_NodeArray[ParentID]._AABBs.GetMaxPos();
	Vec3 min = m_NodeArray[ParentID]._AABBs.GetMinPos();
	Vec3 Div = (max + min) * 0.5;

	for (auto c : m_NodeArray[ParentID]._Sphere){
		Node& tn = child[m_SphereArray[c][0] < Div[0]][m_SphereArray[c][1] < Div[1]][m_SphereArray[c][2] < Div[2]];
		tn._Sphere.push_back(c);
		tn._AABBs.InSphere(m_SphereArray[c], m_SR);
	}
	m_NodeArray[ParentID]._Sphere.clear();

	int cc = 0;
	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 2; j++){
			for (int k = 0; k < 2; k++){
				child[i][j][k]._ID = m_NodeArray.size();
				m_NodeArray[ParentID]._Children[cc] = child[i][j][k]._ID;
				m_NodeArray.push_back(child[i][j][k]);
				cc++;
			}
		}
	}

	for (int i = 0; i < 8; i++){
		PCreateOctree(m_NodeArray[ParentID]._Children[i], depth + 1);
	}

	return true;
}


bool OcTree::CreateOctree(const std::vector<Vec3>& data, const float SphereRadius, int depth){

	m_NodeArray.clear();
	m_SphereArray = data;
	m_SR = SphereRadius;
	m_Depth = depth;

	//ダミーノード
	Node dn;
	dn._ID = 0;
	m_NodeArray.push_back(dn);


	//親ノード生成
	Node parent;
	AABBS aabb;
	int i = 0;
	for (auto e : m_SphereArray){
		aabb.InSphere(e, m_SR);
		parent._Sphere.push_back(i++);
	}
	parent._ID = 1;
	parent._AABBs = aabb;
	parent._ID = 1;
	for (int i = 0; i < 8; i++){
		parent._Children[i] = i + 1;
	}
	m_NodeArray.push_back(parent);
	

	PCreateOctree(1, 1);
	CreateTree(1);

	return true;

}

OcTree::OcTree(){
}

OcTree::~OcTree(){
}