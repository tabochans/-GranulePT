#include"VoxelGrid.h"

using namespace PTUtility;


VoxelGrid::VoxelData* VoxelGrid::GetVoxelFromPos(const Vec3& Pos) const{

	int ID[3];

	for (int i = 0; i < 3; i++){
		ID[i] = (Pos[i] - m_MinPos[i] + 0.50 * GetVoxelWidth()) / GetVoxelWidth();

		if (ID[i] < 0 || ID[i] >= m_NVoxel){
			return NULL;
		}
	}

	return GetVoxelFromID(ID[0], ID[1], ID[2]);

}

VoxelGrid::VoxelData* VoxelGrid::GetVoxelFromID(int IDX, int IDY, int IDZ) const{
	return &m_VoxelData[IDX][IDY][IDZ];
}



bool VoxelGrid::CreateVoxelGridFromMesh(const Scene& scene, int NVoxel){

	float MaxWidth = scene.GetSceneWidth();
	float dw = MaxWidth / (float)(NVoxel);
	Vec3 C = scene.GetSceneCenter();

	InitVoxelGrid(NVoxel + 4);
	m_VoxelWidth = dw;
	m_NVoxel = NVoxel + 4;

	m_MinPos = C - Vec3::Ones() * (dw * (m_NVoxel / 2.0f) - dw / 2.0f);

	for (int i = 0; i < (NVoxel + 4) * (NVoxel + 4); i++){

		RandomMT MT(i + time(NULL) * 1937);

		int w = i % (NVoxel + 4);
		int h = i / (NVoxel + 4);

		int IDX = -1;
		int IDY = -1;

		int Ssample = 10;
		for (int ssx = 0; ssx < Ssample; ssx++){
			for (int ssy = 0; ssy < Ssample; ssy++){
				std::vector<int> ta(NVoxel + 4, -1);

				Vec3 RV(MT.genrand64_real3() * 0.9 * dw + dw * 0.05, MT.genrand64_real3() * 0.9 * dw + dw * 0.05, 0);
				//Vec3 RV(0, 0, 0);

				if (ssx == 0){
					RV[0] = 0.01 * dw;
				}
				else if (ssx == Ssample - 1){
					RV[0] = 0.99 * dw;
				}
				if (ssy == 0){
					RV[1] = 0.01 * dw;
				}
				else if (ssy == Ssample - 1){
					RV[1] = 0.99 * dw;
				}


				Ray r = Ray(
					m_MinPos + Vec3(w, h, -0.0000001) * dw + Vec3::Ones() * dw * 0.5 - RV,
					Vec3(-0.00000001, 0.0000000001, 1));


				std::vector<AbsG::RF> ffa;
				int n = scene.GetFacetFromRay(r, ffa);


				if (n > 0){
					for (int v = 0; v < n; v++){
						Vec3 Pos = r.m_Org + r.m_Dir * ffa[v]._t;
						VoxelData* vv = GetVoxelFromPos(Pos);
						ta[vv->_ID[2]] = 5;
						if (v == 0){
							IDX = GetVoxelFromPos(Pos)->_ID[0];
							IDY = GetVoxelFromPos(Pos)->_ID[1];
						}
						else{
							if (IDX != GetVoxelFromPos(Pos)->_ID[0]){
								throw;
							}
							if (IDY != GetVoxelFromPos(Pos)->_ID[1]){
								throw;
							}
						}
						
					}
				}

				int InVoxel = -1;
				std::vector<int> ida;
				for (int v = 0; v < NVoxel + 4; v++){
					if (ta[v] == 5){
						ida.push_back(v);
					}
				}
				if (ida.size() == 1){
					int iid = ida[0] + 1;
					ta[iid] = 5;
				}
				else if (ida.size() % 2 == 1){
					int iid = ida[(ida.size() + 1) / 2 - 1];
					ta[iid] = -1;
				}
				for (int v = 0; v < NVoxel + 4; v++){
					if (ta[v] == 5){
						ta[v] = 1;
						InVoxel *= -1;
						continue;
					}
					else{
						if (InVoxel == 1){
							ta[v] = 1;
						}
						else{
							ta[v] = -1;
						}
					}

				}
				if (ta.size() > 0){

					for (int v = 0; v < NVoxel + 4; v++){
						if (ta[v] == 1){
							VoxelData* vox = GetVoxelFromID(IDX, IDY, v);
							vox->_Vstat = 1;
						}
					}
				}

			}
		}
	}


	//ExpandVoxel();
	CreateSurfaceVoxel();

	return true;
}


bool VoxelGrid::CreateVoxelGrid(const char* FileName){

	

	std::vector<Vec3> buf;


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
			buf.push_back(Vec3(x, y, z));
		}

	}

	file.close();

	Vec3 Min(100, 100, 100);
	Vec3 Max(-100, -100, -100);

	float ww = 100000;

	for (int i = 0; i < buf.size(); i++){

		float tw = abs(buf[0][0] - buf[i][0]);
		if (tw < ww && tw > 0.000000001){
			ww = tw;
		}
		
		for (int j = 0; j < 3; j++){
			if (buf[i][j] < Min[j]){
				Min[j] = buf[i][j];
			}
			if (buf[i][j] > Max[j]){
				Max[j] = buf[i][j];
			}
		}
	}

	m_MinPos = Min - Vec3::Ones() * ww * 2;

	float MM = std::max(Max[0], std::max(Max[1], Max[2]));
	float mm = std::min(Min[0], std::min(Min[1], Min[2]));

	m_NVoxel = 4 + (MM - mm + 0.00000001) / ww;
	m_VoxelWidth = ww;

	InitVoxelGrid(m_NVoxel);


	for (int i = 0; i < buf.size(); i++){

		VoxelData* pV = GetVoxelFromPos(buf[i]);
		if (pV){
			pV->_Vstat = 1;
		}
	}

	m_VoxelWidth *= 1.00001;

	return true;
}

bool VoxelGrid::CreateSurfaceVoxel(){
	
	for (int i = 2; i < m_NVoxel - 2; i++){
		for (int j = 2; j < m_NVoxel - 2; j++){
			for (int k = 2; k < m_NVoxel - 2; k++){


				VoxelData* vp = GetVoxelFromID(i, j, k);
				if (vp->_Vstat == 1){

					bool IfInter = true;
					for (int dx = -1; dx < 2; dx+=2){
						for (int dy = -1; dy < 2; dy+=2){
							for (int dz = -1; dz < 2; dz+=2){
								VoxelData* v = GetVoxelFromID(i + dx, j + dy, k + dz);
								if (v->_Vstat == -1){
									IfInter = false;
									goto LOOP_END;
								}
							}
						}
					}
					LOOP_END: 
					if (!IfInter){
						vp->_Vstat = 2;
					}
				}


			}
		}
	}

	return false;
}

bool VoxelGrid::ExpandVoxel(){



	for (int i = 2; i < m_NVoxel - 2; i++){
		for (int j = 2; j < m_NVoxel - 2; j++){
			for (int k = 2; k < m_NVoxel - 2; k++){

				VoxelData* vp = GetVoxelFromID(i, j, k);
				if (vp->_Vstat == 1){

					for (int dx = -1; dx < 2; dx++){
						for (int dy = -1; dy < 2; dy++){
							for (int dz = -1; dz < 2; dz++){
								if (dx == 0 && dy == 0 && dz == 0){ continue; }
								VoxelData* v = GetVoxelFromID(i + dx, j + dy, k + dz);
								if (v->_Vstat == -1){
									v->_Vstat = 2;
								}
							}
						}
					}
				}


			}
		}
	}
	for (int i = 2; i < m_NVoxel - 2; i++){
		for (int j = 2; j < m_NVoxel - 2; j++){
			for (int k = 2; k < m_NVoxel - 2; k++){

				VoxelData* vp = GetVoxelFromID(i, j, k);
				if (vp->_Vstat == 2){
					vp->_Vstat = 1;
				}


			}
		}
	}


	return false;
}




bool VoxelGrid::InitVoxelGrid(int w){
	m_VoxelData.clear();

	VoxelData VD;
	VD._Vstat = -1;

	for (int i = 0; i < w; i++){
		std::vector<std::vector<VoxelData>> vv;
		for (int j = 0; j < w; j++){

			std::vector<VoxelData> v;
			for (int k = 0; k < w; k++){
				VD._ID[0] = -1;
				VD._ID[1] = -1;
				VD._ID[2] = -1;
				v.push_back(VD);
			}
			vv.push_back(v);
		}
		m_VoxelData.push_back(vv);
	}


	for (int i = 0; i < w; i++){
		for (int j = 0; j < w; j++){
			for (int k = 0; k < w; k++){
				m_VoxelData[i][j][k]._ID[0] = i;
				m_VoxelData[i][j][k]._ID[1] = j;
				m_VoxelData[i][j][k]._ID[2] = k;
			}
		}
	}


	return true;
}


VoxelGrid::VoxelGrid(int Width) : m_NVoxel(Width){
}

VoxelGrid::~VoxelGrid(){

}


Vec3 VoxelGrid::GetCenterFromID(int IDX, int IDY, int IDZ)const{

	return Vec3((float)IDX, (float)IDY, (float)IDZ)*GetVoxelWidth() + m_MinPos;

}

Vec3 VoxelGrid::GetCenterPos(const VoxelData* v){

	return GetCenterFromID(v->_ID[0], v->_ID[1], v->_ID[2]);

}

int VoxelGrid::Intersect(const Ray& r, std::vector<const VPP>& result){

	std::vector<VPP> temp;
	VPP vpp;

	result.clear();

	if ((vpp._pp = GetVoxelFromPos(r.m_Org)) != NULL){
		if (vpp._pp->_Vstat > 0){

			float t = 1000;
			float f = -1000;

			VoxelData* pp = GetVoxelFromPos(GetCenterPos(vpp._pp));


			if (!Intersect::Voxel_Ray(GetCenterPos(vpp._pp), GetVoxelWidth()*1.0001, r, t, f)){
				return 0;
			}
			vpp._f = f;
			vpp._t = t;
			result.push_back(vpp);
			return 1;
		}
	}
	for (int i = 0; i < m_NVoxel; i++){
		for (int j = 0; j < m_NVoxel; j++){
			for (int k = 0; k < m_NVoxel; k++){

				if (m_VoxelData[i][j][k]._Vstat > 0){
					float t = 1000;
					float f = -1000;

					if (Intersect::Voxel_Ray(GetCenterFromID(i, j, k), GetVoxelWidth()*1.001, r, t, f)){
						vpp._t = t;
						vpp._f = f;
						vpp._pp = GetVoxelFromID(i, j, k);
						temp.push_back(vpp);
					}

				}

			}
		}
	}

	if (temp.size() > 0){
		std::sort(temp.begin(), temp.end());

		if (GetVoxelFromPos(r.m_Org)){
			bool fg = false;
			for (int i = 0; i < temp.size(); i++){
				
				if (!fg){
					if (temp[i]._f > 0.00001){
						fg = true;
					}
				}

				if (fg){ 
					result.push_back(temp[i]); 
					return 1;
				}
			}
		}
		else{
			result.push_back(temp[0]);
		}


		return result.size();
	}

	return 0;

}

void VoxelGrid::PrintVoxel(const Vec3& Center, std::ofstream& file){

	static int _NF = 0;

	float dx, dy, dz;
	dx = dy = dz = m_VoxelWidth * 0.5f;
	std::vector<Vec3> _VList;
	std::vector<std::vector<int>> _FList;

	_VList.push_back(Center + Vec3(-dx, -dy, dz));
	_VList.push_back(Center + Vec3(-dx, dy, dz));
	_VList.push_back(Center + Vec3(dx, dy, dz));
	_VList.push_back(Center + Vec3(dx, -dy, dz));
	_VList.push_back(Center + Vec3(-dx, -dy, -dz));
	_VList.push_back(Center + Vec3(-dx, dy, -dz));
	_VList.push_back(Center + Vec3(dx, dy, -dz));
	_VList.push_back(Center + Vec3(dx, -dy, -dz));

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


	for (int v = 0; v < _VList.size(); v++){
		file << "v " << _VList[v][0] << " " << _VList[v][1] << " " << _VList[v][2] << std::endl;
	}
	for (int f = 0; f < _FList.size(); f++){
		file << "f ";
		for (int ff = 0; ff < _FList[0].size(); ff++){
			file << _FList[f][ff] << " ";
		}
		file << std::endl;
	}

}


void VoxelGrid::WriteVoxelData(std::ofstream& file){

	

	for (int i = 0; i < m_NVoxel; i++){
		for (int j = 0; j < m_NVoxel; j++){
			for (int k = 0; k < m_NVoxel; k++){

				VoxelData* vvv = GetVoxelFromID(i, j, k);
				Vec3 p = GetCenterFromID(i, j, k);
				if (vvv->_Vstat == 2){
					PrintVoxel(p, file);
					//file << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
				}
			}
		}
	}


}