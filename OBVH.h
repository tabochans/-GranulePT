#pragma once


#include"Facet.h"
#include"AbsBVH.h"
#include"AABB_S.h"
#include"Vec3.h"
#include<algorithm>
#include<memory>
#include<ostream>
#include<deque>
#include <immintrin.h>


namespace PTUtility{

	//OBVH ‚æ‚­‚í‚©‚ñ‚È‚¢‚¯‚ÇSIMD‰»‚³‚ê‚Ä‚¢‚é‚ç‚µ‚¢
	class OBVH : public AbsG{

	private:

		const float T_AABB;
		const float T_TRI;
		const int MAX_DEPTH;
		const int MIN_ELEMENT;
		const float ZEROP;

		std::vector<PTUtility::Facet> m_data;
		Vec3 m_Max;
		Vec3 m_Min;

		enum DIVMODE{
			FIRST = 0,
			SECOND = 1,
			THIRD = 3
		};

		struct Node{

			std::vector<PTUtility::Facet*> _facets;
			std::vector<PTUtility::AABBS> _AABBs;
			int _Children[8];
			int _NumChildren;

			__m256 CAABB[2][3];

			Node(){
				_Children[0] = 0;
				_Children[1] = 0;
				_Children[2] = 0;
				_Children[3] = 0;
				_Children[4] = 0;
				_Children[5] = 0;
				_Children[6] = 0;
				_Children[7] = 0;
				_NumChildren = 0;
			}
			~Node(){
			}

		};
		std::deque<Node> m_NArray;

		void Divide(Node& node, int GID, int mode);



		float CostAABB(
			const PTUtility::AABBS& Parent,
			float sA, int NumTriangleA,
			float sB, int NumTriangleB) const;

		float ZeroCost(const PTUtility::AABBS& Parent, int NumElements) const;

		void PWriteObj(std::ostream& ost, int PID) const;


		mutable std::vector<std::vector<int>> _FList;
		mutable std::vector<Vec3> _VList;
		mutable int _NF;
		mutable int _SS;
		mutable int _SSS;

		int PGetElements(const PTUtility::Ray& ray, std::vector<RF>& Facets, int ID) const;
		int PGetElements(const Vec3& Center, float R, std::vector<RF>& Facets, int ID) const;
		int CreateAABBTree(int ID);

	protected:

	public:

		virtual const std::string Print() const{
			return std::string("OBVH");
		}
		virtual float GetSceneWidth();
		virtual Vec3 GetSceneCenter();
		virtual void CreateGraph(const std::vector<PTUtility::Facet>& FacetArray);
		virtual void ClearGraph();
		virtual void PrintGraph(std::ostream& ost) const;
		virtual void WriteObj(std::ostream& ost) const;
		virtual int GetElements(const PTUtility::Ray& ray, std::vector<RF>& Facets) const;
		virtual int GetElements(const Vec3& Center, float R, std::vector<RF>& Facets);

		OBVH();
		virtual ~OBVH();

	};



}