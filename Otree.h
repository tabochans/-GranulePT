#pragma once

#include"Vec3.h"
#include"AABB_S.h"
#include"Intersect.h"
#include"Ray.h"
#include<vector>
#include <immintrin.h>
#include<algorithm>

namespace PTUtility{

	//�ق�Ƃ�AbsBVH����p�������������n
	//�I�N�c���[�@������
	class OcTree{

	public:

		struct SC{
			Vec3 _Pos;
			float _t;
			int _SID;

			bool operator<(const SC& ref){
				return this->_t < ref._t;
			}
			bool operator>(const SC& ref){
				return this->_t > ref._t;
			}
			bool operator<=(const SC& ref){
				return this->_t <= ref._t;
			}
			bool operator>=(const SC& ref){
				return this->_t >= ref._t;
			}

		};

	private:

		struct Node{

			int _ID;
			std::vector<int> _Children;
			std::vector<int> _Sphere;
			AABBS _AABBs;
			__m256 CAABB[2][3];
			bool IfLeaf;

			Node() : _Children(8,0), IfLeaf(false){
				
			}
		};

		//���̃f�[�^�z��
		std::vector<Vec3> m_SphereArray;
		float m_SR;

		//Node�z��
		std::vector<Node> m_NodeArray;
		
		//�ő�K�w
		int m_Depth;
		

		int PGetSphereFromRay(int ID, const Ray& ray, std::vector<SC>& result);
		bool PCreateOctree(int ParentID, int depth);

		int CreateTree(int PID);

	protected:


	public:

		int GetSphereFromRay(const Ray& ray, std::vector<SC>& result);
		bool CreateOctree(const std::vector<Vec3>& data, const float SphereRadius, int depth);

		OcTree();
		virtual ~OcTree();

	};



}