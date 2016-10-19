#include"OBVH.h"

using namespace PTUtility;

void  OBVH::Divide(Node& node, int GID, int Mode){

	int NElm = node._facets.size();

	//分割軸
	int Axis = -1;
	//分割場所leftが含むポリゴン数
	int Dpos = -1;
	//分割コスト
	float Cost = T_TRI * NElm;

	std::vector<Facet*> FA[3];
	FA[0] = node._facets;
	FA[1] = node._facets;
	FA[2] = node._facets;

	std::sort(FA[0].begin(), FA[0].end(), &Facet::compFacetXP);
	std::sort(FA[1].begin(), FA[1].end(), &Facet::compFacetYP);
	std::sort(FA[2].begin(), FA[2].end(), &Facet::compFacetZP);



	//三軸について
	for (int ax = 0; ax < 3; ax++){

		//表面積の控え
		//SufArray[i] には、i個のポリゴンが含まれた場合の表面積ね
		std::vector<float> SufArray(NElm + 2);

		AABBS Bright, Bleft;


		//rightの表面積を計算しとく
		SufArray[0] = 1000000;
		Bleft.Reset();
		for (int i = NElm - 1 ; i >= 0; i--){
			Bright.InTriangle(FA[ax][i]->GetPos(0), FA[ax][i]->GetPos(1), FA[ax][i]->GetPos(2));
			SufArray[NElm - i] = Bright.GetArea();
		}
		Bright.Reset();
		Bleft.Reset();

		//SAHの計算
		for (int i = 0; i < NElm; i++){
			Bleft.InTriangle(FA[ax][i]->GetPos(0), FA[ax][i]->GetPos(1), FA[ax][i]->GetPos(2));
			float tCost = CostAABB(node._AABBs[0], Bleft.GetArea(), i + 1, SufArray[NElm - i - 1], NElm - i - 1);
			if (Cost > tCost){
				Axis = ax;
				Dpos = i + 1;
				Cost = tCost;
			}
		}
	}

	
	if (Axis != -1 && Dpos > 0 && Dpos < NElm && NElm > MIN_ELEMENT){

		//分割先が見つかった

		Node Nright, Nleft;
		Nleft._facets.resize(Dpos);
		Nright._facets.resize(NElm - Dpos);
		AABBS Bright, Bleft;
		for (int i = 0; i < Dpos; i++){
			Bleft.InTriangle(FA[Axis][i]->GetPos(0), FA[Axis][i]->GetPos(1), FA[Axis][i]->GetPos(2));
			Nleft._facets[i] = FA[Axis][i];
		}
		for (int i = Dpos; i < NElm; i++){
			Bright.InTriangle(FA[Axis][i]->GetPos(0), FA[Axis][i]->GetPos(1), FA[Axis][i]->GetPos(2));
			Nright._facets[i - Dpos] = FA[Axis][i];
		}
		Nleft._AABBs.push_back(Bleft);
		Nright._AABBs.push_back(Bright);


		//OBVHだから、さらに二回分割をする
		if (Mode == 0){
			Divide(Nleft, GID, 1);
			Divide(Nright, GID, 1);
			std::vector<Facet*>().swap(node._facets);
		}
		if (Mode == 1){
			Divide(Nleft, GID, 2);
			Divide(Nright, GID, 2);
			std::vector<Facet*>().swap(node._facets);
			return;
		}
		if (Mode == 2){
			int cid = m_NArray.size();
			m_NArray.push_back(Nleft);
			m_NArray.push_back(Nright);
			
			m_NArray[GID]._Children[m_NArray[GID]._NumChildren] = cid;
			m_NArray[GID]._NumChildren++;

			m_NArray[GID]._Children[m_NArray[GID]._NumChildren] = cid + 1;
			m_NArray[GID]._NumChildren++;

			std::vector<Facet*>().swap(node._facets);

			return;
		}

		//そして時は動き出す
		for (int i = 0; i < 8; i++){
			int cd = m_NArray[GID]._Children[i];
			if (cd > 0){ Divide(m_NArray[cd], cd, 0); }
		}
	

	}
	else{
		//分割しない

		if (Mode == 1 || Mode == 2){
			int cid = m_NArray.size();
			m_NArray.push_back(node);
			m_NArray[GID]._Children[m_NArray[GID]._NumChildren] = cid;
			m_NArray[GID]._NumChildren++;
		}

	}

}


float OBVH::CostAABB(
	const PTUtility::AABBS& Parent,
	float sA, int NumTriangleA,
	float sB, int NumTriangleB) const{

	return 2.0f * T_AABB +
		NumTriangleA * T_TRI * (sA) / (Parent.GetArea()) +
		NumTriangleB * T_TRI * (sB) / (Parent.GetArea());

}

float OBVH::ZeroCost(const PTUtility::AABBS& Parent, int NumElements) const{
	return CostAABB(Parent, Parent.GetArea(), NumElements, 10000000000000, 0);
}

int OBVH::CreateAABBTree(int ID){
	//親ノードが子のAABB８個を一括で持つよーにする


	if (m_NArray[ID]._NumChildren == 0){
		return 1;
	}
	else{

		float AABBMax[3][8];
		float AABBMin[3][8];

		for (int i = 0; i < 8; i++){

			if (m_NArray[ID]._Children[i] > 0){
				const Node& cn = m_NArray[m_NArray[ID]._Children[i]];
				for (int j = 0; j < 3; j++){
					AABBMax[j][i] = cn._AABBs[0].GetMaxPos()[j];
					AABBMin[j][i] = cn._AABBs[0].GetMinPos()[j];
				}
			}
			else{
				for (int j = 0; j < 3; j++){
					AABBMax[j][i] = -1000000;
					AABBMin[j][i] = 1000000;
				}
			}
		}

		for (int i = 0; i < 3; i++){
			m_NArray[ID].CAABB[0][i] = _mm256_loadu_ps(AABBMin[i]);
			m_NArray[ID].CAABB[1][i] = _mm256_loadu_ps(AABBMax[i]);
		}
		m_NArray[ID]._AABBs.clear();

		for (int i = 0; i < 8; i++){
			if (m_NArray[ID]._Children[i] > 0){
				CreateAABBTree(m_NArray[ID]._Children[i]);
			}
		}

	}

	return 1;
}

void OBVH::CreateGraph(const std::vector<PTUtility::Facet>& FacetArray){

	m_NArray.clear();

	Node parent;
	AABBS pAB;

	parent._AABBs.push_back(pAB);
	m_NArray.push_back(parent);
	parent._AABBs.clear();

	m_data = FacetArray;

	for (int i = 0; i < m_data.size(); i++){
		pAB.InTriangle(m_data[i].GetPos(0), m_data[i].GetPos(1), m_data[i].GetPos(2));
	}
	m_Max = pAB.GetMaxPos();
	m_Min = pAB.GetMinPos();
	parent._AABBs.push_back(pAB);
	parent._facets.resize(m_data.size());
	for (int i = 0; i < m_data.size(); i++){
		parent._facets[i] = &m_data[i];
	}

	m_NArray.push_back(parent);

	Divide(m_NArray[1], 1, 0);

	CreateAABBTree(1);

	return;

}

void OBVH::ClearGraph(){
}


void OBVH::PWriteObj(std::ostream& ost, int PID) const{

	if (m_NArray[PID]._facets.size() > 0){

		AABBS A = m_NArray[PID]._AABBs[0];
		Vec3 c = A.GetCenter();
		float dx, dy, dz;
		dx = A.GetWidth()[0] * 0.5f;
		dy = A.GetWidth()[1] * 0.5f;
		dz = A.GetWidth()[2] * 0.5f;

		_VList.push_back(c + Vec3(-dx, -dy, dz));
		_VList.push_back(c + Vec3(-dx, dy, dz));
		_VList.push_back(c + Vec3(dx, dy, dz));
		_VList.push_back(c + Vec3(dx, -dy, dz));
		_VList.push_back(c + Vec3(-dx, -dy, -dz));
		_VList.push_back(c + Vec3(-dx, dy, -dz));
		_VList.push_back(c + Vec3(dx, dy, -dz));
		_VList.push_back(c + Vec3(dx, -dy, -dz));

		std::vector<int> tt;

		tt.push_back(1 + _NF * 8); tt.push_back(4 + _NF * 8); tt.push_back(3 + _NF * 8); tt.push_back(2 + _NF * 8);
		_FList.push_back(tt);
		tt.clear();

		tt.push_back(1 + _NF * 8); tt.push_back(5 + _NF * 8); tt.push_back(8 + _NF * 8); tt.push_back(4 + _NF * 8);
		_FList.push_back(tt);
		tt.clear();

		tt.push_back(2 + _NF * 8); tt.push_back(6 + _NF * 8); tt.push_back(5 + _NF * 8); tt.push_back(1 + _NF * 8);
		_FList.push_back(tt);
		tt.clear();

		tt.push_back(3 + _NF * 8); tt.push_back(7 + _NF * 8); tt.push_back(6 + _NF * 8); tt.push_back(2 + _NF * 8);
		_FList.push_back(tt);
		tt.clear();

		tt.push_back(4 + _NF * 8); tt.push_back(8 + _NF * 8); tt.push_back(7 + _NF * 8); tt.push_back(3 + _NF * 8);
		_FList.push_back(tt);
		tt.clear();

		tt.push_back(6 + _NF * 8); tt.push_back(7 + _NF * 8); tt.push_back(8 + _NF * 8); tt.push_back(5 + _NF * 8);
		_FList.push_back(tt);
		tt.clear();

		_NF++;
		_SSS += m_NArray[PID]._facets.size();
	}


	for (int i = 0; i < 8; i++){
		if (m_NArray[PID]._Children[i] > 0){ PWriteObj(ost, m_NArray[PID]._Children[i]); }
	}



}

void OBVH::PrintGraph(std::ostream& ost) const{ 
}


void OBVH::WriteObj(std::ostream& ost) const{

	_VList.clear();
	_FList.clear();
	_NF = 0;
	_SSS = 0;

	ost << " o BVH" << std::endl;

	PWriteObj(ost, 1);

	for (int v = 0; v < _VList.size(); v++){
		ost << "v " << _VList[v][0] << " " << _VList[v][1] << " " << _VList[v][2] << std::endl;
	}
	for (int f = 0; f < _FList.size(); f++){
		ost << "f ";
		for (int ff = 0; ff < _FList[0].size(); ff++){
			ost << _FList[f][ff] << " ";
		}
		ost << std::endl;
	}

}


int OBVH::PGetElements(const PTUtility::Ray& ray, std::vector<RF>& Facets, int ID) const{
	
	float u, v, t;
	u = v = t = 0.0f;
	if (m_NArray[ID]._NumChildren == 0){
		int nc = m_NArray[ID]._facets.size();
		
		for (int i = 0; i < nc; i++){
			if (Intersect::Ray_Triangle(ray, *m_NArray[ID]._facets[i], u, v, t)){
				if (t > FLT_MIN){
					Facets.push_back(AbsG::RF(m_NArray[ID]._facets[i], u, v, t));
				}
			}
		}
		return Facets.size();
	}

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


	if (Intersect::IsHit(m_NArray[ID].CAABB, org, idir, sign, mask)){
		for (int i = 0; i < 8; i++){
			int fg = mask & (1 << (i));
			if ((fg) > 0){
				int iid = m_NArray[ID]._Children[i];
				if (iid > 0){ PGetElements(ray, Facets, iid); }
			}
		}
	}

	return Facets.size();
}


int OBVH::GetElements(const PTUtility::Ray& ray, std::vector<RF>& Facets) const{
	
	Facets.clear();
	
	int n = PGetElements(ray, Facets, 1);
	std::sort(Facets.begin(), Facets.end());
	return n;
}

int OBVH::PGetElements(const Vec3& Center, float R, std::vector<RF>& Facets, int ID) const{

	float u, v, t;
	u = v = t = 0.0f;
	if (m_NArray[ID]._NumChildren == 0){
		int nc = m_NArray[ID]._facets.size();

		for (int i = 0; i < nc; i++){
			if (Intersect::Sphere_Triangle(Center,R, *m_NArray[ID]._facets[i], u, v, t)){
				Facets.push_back(AbsG::RF(m_NArray[ID]._facets[i], u, v, t));
			}
		}
		return Facets.size();
	}
	else{

		__m256 Pos[3];
		__m256 Rd[1];

		int mask = 0;

		for (int i = 0; i < 3; i++){
			Pos[i] = _mm256_set1_ps(Center[i]);
		}
		Rd[0] = _mm256_set1_ps(R);

		if (Intersect::Sphere_AABB_Fast_S(Pos, Rd, m_NArray[ID].CAABB, mask)){
			for (int i = 0; i < 8; i++){
				int fg = mask & (1 << (i));
				if ((fg) > 0){
					int iid = m_NArray[ID]._Children[i];
					if (iid > 0){ PGetElements(Center, R, Facets, iid); }
				}
			}
		}

		return Facets.size();

	}
	

}

int OBVH::GetElements(const Vec3& Center, float R, std::vector<RF>& Facets){
	Facets.clear();
	int n = PGetElements(Center, R, Facets, 1);
	std::sort(Facets.begin(), Facets.end());
	return n;
}

float OBVH::GetSceneWidth(){

	Vec3 w = m_Max - m_Min;
	return std::max(w[0], std::max(w[1], w[2]));
}

Vec3 OBVH::GetSceneCenter(){
	return 0.5 * (m_Max + m_Min);
}


OBVH::OBVH() : T_AABB(1.0f), T_TRI(1.0f), MAX_DEPTH(10), MIN_ELEMENT(1), ZEROP(FLT_MIN){
}

OBVH::~OBVH(){
}