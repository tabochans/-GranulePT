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


		//���A���񂷂ƁH�@�܂�������B�B�B
		VoxelData* GetVoxelFromPos(const Vec3& Pos) const;
		VoxelData* GetVoxelFromID(int IDX, int IDY, int IDZ) const;
		Vec3 GetCenterFromID(int IDX, int IDY, int IDZ)const;

		bool InitVoxelGrid(int w);
		bool CreateSurfaceVoxel();
		bool ExpandVoxel();

	private:

		void PrintVoxel(const Vec3& Center, std::ofstream& file);

		//�{�N�Z���̕�
		float m_VoxelWidth;

		//�ꎲ�ɂȂ�ԃ{�N�Z����
		int m_NVoxel;

		//�C���f�b�N�X(0,0,0)�̃{�N�Z���̒��S���W�ł���A�ł����������W
		Vec3 m_MinPos;

		mutable std::vector<std::vector<std::vector<VoxelData>>> m_VoxelData;


	public:
		explicit VoxelGrid(int Width);
		virtual ~VoxelGrid();

		//�{�N�Z���̕��@�����H��������ˁH
		float GetVoxelWidth()const{ return m_VoxelWidth; }

		Vec3 GetCenterPos(const VoxelData* v);

		//Ray���N�G���Ƃ��āA�������邷�ׂẴ{�N�Z���̃|�C���^��Ԃ��܂�
		int Intersect(const Ray& r, std::vector<const VPP>& result);
	

		//�|�C���g�N���E�h�f�[�^����A�O���b�h���쐬����
		bool CreateVoxelGrid(const char* FileName);

		//���b�V�����{�N�Z���C�Y���č쐬����
		bool CreateVoxelGridFromMesh(const Scene& scene, int NVoxel);


		void WriteVoxelData(std::ofstream& file);

	};

}