#pragma once
#include"Ray.h"
#include"Vec3.h"
#include"Intersect.h"
#include"Scene.h"
#include<sstream>
#include<algorithm>
#include<fstream>
#include<vector>

namespace PTUtility{

	class VoxelGrid{

	public:
		struct VoxelData{
			int _ID[3];
			
			int _Vstat;
			VoxelData(){
				_ID[0] = -1;
				_ID[1] = -1;
				_ID[2] = -1;
				_Vstat = -1;
			}
		};

		struct VPP{
			VoxelData* _pp;
			float _t;
			float _f;
			static bool Comp(const VPP& left, const VPP& right){
				return left._t < right._t;
			}


			bool operator<(const VPP& ref){
				return this->_t < ref._t;
			}
			bool operator>(const VPP& ref){
				return this->_t > ref._t;
			}
			bool operator<=(const VPP& ref){
				return this->_t <= ref._t;
			}
			bool operator>=(const VPP& ref){
				return this->_t >= ref._t;
			}



		};


	private:


		//こ、こんすと？　まあいいや。。。
		VoxelData* GetVoxelFromPos(const Vec3& Pos) const;
		VoxelData* GetVoxelFromID(int IDX, int IDY, int IDZ) const;
		Vec3 GetCenterFromID(int IDX, int IDY, int IDZ)const;

		bool InitVoxelGrid(int w);
		bool CreateSurfaceVoxel();
		bool ExpandVoxel();

	private:

		void PrintVoxel(const Vec3& Center, std::ofstream& file);

		//ボクセルの幅
		float m_VoxelWidth;

		//一軸にならぶボクセル数
		int m_NVoxel;

		//インデックス(0,0,0)のボクセルの中心座標であり、最も小さい座標
		Vec3 m_MinPos;

		mutable std::vector<std::vector<std::vector<VoxelData>>> m_VoxelData;


	public:
		explicit VoxelGrid(int Width);
		virtual ~VoxelGrid();

		//ボクセルの幅　いい？幅だからね？
		float GetVoxelWidth()const{ return m_VoxelWidth; }

		Vec3 GetCenterPos(const VoxelData* v);

		//Rayをクエリとして、交差するすべてのボクセルのポインタを返します
		int Intersect(const Ray& r, std::vector<const VPP>& result);
	

		//ポイントクラウドデータから、グリッドを作成する
		bool CreateVoxelGrid(const char* FileName);

		//メッシュをボクセライズして作成する
		bool CreateVoxelGridFromMesh(const Scene& scene, int NVoxel);


		void WriteVoxelData(std::ofstream& file);

	};

}