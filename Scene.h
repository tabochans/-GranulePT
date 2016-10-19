#pragma once

#define TINYOBJLOADER_IMPLEMENTATION
#include"tiny_obj_loader.h"
#include"Facet.h"
#include"AbsBVH.h"
#include"Texture2D.h"
#include"MT.h"
#include<vector>
#include<map>
#include<iostream>
#include<memory>
#include<deque>

namespace PTUtility{

	class Scene{

	public:

		static float DiffuseBRDF(const Vec3& in, const Vec3& out){
			return 1.0f / M_PI;
		}

	public:

		struct LightMesh{
			float _Area;
			Vec3 _O;
			Vec3 _V;
			Vec3 _U;
			Vec3 _Color;
			Vec3 _Normal;
			float _Power;
		};

	private:
		AbsG* m_SceneData;

		std::map<std::string, std::shared_ptr<Texture2D> > m_TextureMap;
		std::deque<LightMesh> m_LightMeshs;
		float m_SumLightArea;

	public:
		int IfExistLight()const{ return m_LightMeshs.size() > 0; }
		float GetSumLightArea() const;
	
	public:

		bool CreateScene(const char* ObjFileName);
		void WriteBVH(std::ostream& ost) const;
		float GetSceneWidth() const;
		Vec3 GetSceneCenter() const;

		int GetFacetFromRay(const PTUtility::Ray& ray, std::vector<AbsG::RF>& Facets) const;
		const LightMesh& GetLightMesh(RandomMT& MT, float& spdf, long int ID = -1) const;
		bool IfDeeperThanDepth(const Vec3& Pos, float Depth) const;

		const Vec3 GetTextureColor(const std::string& TextureName, const Vec2& UV) const;


		//îCà”ÇÃÉfÅ[É^ç\ë¢Ç(ÇΩÇ‘ÇÒ)ìnÇπÇÈäÛÉKÉX
		explicit Scene(AbsG* DataStructure);
		
		virtual ~Scene();

	};

}