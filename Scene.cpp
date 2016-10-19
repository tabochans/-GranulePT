#include"Scene.h"

using namespace PTUtility;




float Scene::GetSumLightArea() const{
	return m_SumLightArea;
}

const Scene::LightMesh& Scene::GetLightMesh(RandomMT& MT, float& spdf, long int ID) const{
	spdf = 0.0f;
	if (ID < 0){
		ID = 0;
		float sa = 0.0f;
		float ra = MT.genrand64_real2() * (m_SumLightArea - 0.0000001);
		for (int i = 0; i < m_LightMeshs.size(); i++){
			sa += m_LightMeshs[i]._Area;
			if (ra < sa){
				spdf = m_LightMeshs[i]._Area;
				return m_LightMeshs[i];
			}
		}
	}
	throw;
	return m_LightMeshs[ID];

}

bool Scene::IfDeeperThanDepth(const Vec3& Pos, float Depth) const{
	std::vector<AbsG::RF> Facets;
	int n = m_SceneData->GetElements(Pos, Depth, Facets);
	return n == 0;
}


bool Scene::CreateScene(const char* ObjFileName){

	std::cout << "Create Scene START" << std::endl;

	std::cout << "Load OBJ SATRT" << std::endl;
	time_t stt = time(NULL);

	std::vector<Facet> Facets;

	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	std::string err;
	bool ret = tinyobj::LoadObj(shapes, materials, err, ObjFileName, NULL, 1);

	if (!err.empty()) {
		std::cerr << err << std::endl;
	}
	if (!ret) {
		std::cerr << "Failed to load/parse .obj.\n" << std::endl;
		return false;
	}


	for (int m = 0; m < materials.size(); m++){
		if (!materials[m].diffuse_texname.empty()){
			m_TextureMap[materials[m].name] = std::shared_ptr<Texture2D>(new Texture2D(materials[m].diffuse_texname.c_str()));
		}
	}

	

	for (int sps = 0; sps < shapes.size(); sps++){

		for (int idx = 0, face = 0; idx + 2 < shapes[0].mesh.indices.size(); idx += 3){

			int mtid = shapes[sps].mesh.material_ids[face];

			if (materials.size() > 0){

				float r, g, b;
				r = materials[mtid].diffuse[0];
				g = materials[mtid].diffuse[1];
				b = materials[mtid].diffuse[2];

				Facet facet(
					Vertex(
					Vec3(
					shapes[sps].mesh.positions[0 + 3 * shapes[0].mesh.indices[idx + 0]],
					shapes[sps].mesh.positions[1 + 3 * shapes[0].mesh.indices[idx + 0]],
					shapes[sps].mesh.positions[2 + 3 * shapes[0].mesh.indices[idx + 0]]),

					Vec3(
					shapes[sps].mesh.normals[0 + 3 * shapes[0].mesh.indices[idx + 0]],
					shapes[sps].mesh.normals[1 + 3 * shapes[0].mesh.indices[idx + 0]],
					shapes[sps].mesh.normals[2 + 3 * shapes[0].mesh.indices[idx + 0]]),

					Vec2(
					shapes[sps].mesh.texcoords[0 + 2 * shapes[0].mesh.indices[idx + 0]],
					shapes[sps].mesh.texcoords[1 + 2 * shapes[0].mesh.indices[idx + 0]])
					),

					Vertex(
					Vec3(
					shapes[sps].mesh.positions[0 + 3 * shapes[0].mesh.indices[idx + 1]],
					shapes[sps].mesh.positions[1 + 3 * shapes[0].mesh.indices[idx + 1]],
					shapes[sps].mesh.positions[2 + 3 * shapes[0].mesh.indices[idx + 1]]),

					Vec3(
					shapes[sps].mesh.normals[0 + 3 * shapes[0].mesh.indices[idx + 1]],
					shapes[sps].mesh.normals[1 + 3 * shapes[0].mesh.indices[idx + 1]],
					shapes[sps].mesh.normals[2 + 3 * shapes[0].mesh.indices[idx + 1]]),

					Vec2(
					shapes[sps].mesh.texcoords[0 + 2 * shapes[0].mesh.indices[idx + 1]],
					shapes[sps].mesh.texcoords[1 + 2 * shapes[0].mesh.indices[idx + 1]])
					),

					Vertex(
					Vec3(
					shapes[sps].mesh.positions[0 + 3 * shapes[0].mesh.indices[idx + 2]],
					shapes[sps].mesh.positions[1 + 3 * shapes[0].mesh.indices[idx + 2]],
					shapes[sps].mesh.positions[2 + 3 * shapes[0].mesh.indices[idx + 2]]),

					Vec3(
					shapes[sps].mesh.normals[0 + 3 * shapes[0].mesh.indices[idx + 2]],
					shapes[sps].mesh.normals[1 + 3 * shapes[0].mesh.indices[idx + 2]],
					shapes[sps].mesh.normals[2 + 3 * shapes[0].mesh.indices[idx + 2]]),

					Vec2(
					shapes[sps].mesh.texcoords[0 + 2 * shapes[0].mesh.indices[idx + 2]],
					shapes[sps].mesh.texcoords[1 + 2 * shapes[0].mesh.indices[idx + 2]])
					),

					Vec3(r, g, b), Facet::ATT_MATERIAL::DIFFUSE, materials[mtid].name
					);

				
				//ambientをemissionとして用いる。。。
				if (materials[mtid].ambient[0] + materials[mtid].ambient[1] + materials[mtid].ambient[2] > 0.1){
					const Vec3 em(materials[mtid].ambient[0], materials[mtid].ambient[1], materials[mtid].ambient[2]);
					facet.SetLe(600.0f, em);

					LightMesh L;
					L._Color = em;
					L._U = Vec3(
						shapes[sps].mesh.positions[0 + 3 * shapes[0].mesh.indices[idx + 1]],
						shapes[sps].mesh.positions[1 + 3 * shapes[0].mesh.indices[idx + 1]],
						shapes[sps].mesh.positions[2 + 3 * shapes[0].mesh.indices[idx + 1]]) - 
						Vec3(
						shapes[sps].mesh.positions[0 + 3 * shapes[0].mesh.indices[idx + 0]],
						shapes[sps].mesh.positions[1 + 3 * shapes[0].mesh.indices[idx + 0]],
						shapes[sps].mesh.positions[2 + 3 * shapes[0].mesh.indices[idx + 0]]);
					L._V = Vec3(
						shapes[sps].mesh.positions[0 + 3 * shapes[0].mesh.indices[idx + 2]],
						shapes[sps].mesh.positions[1 + 3 * shapes[0].mesh.indices[idx + 2]],
						shapes[sps].mesh.positions[2 + 3 * shapes[0].mesh.indices[idx + 2]]) -
						Vec3(
						shapes[sps].mesh.positions[0 + 3 * shapes[0].mesh.indices[idx + 0]],
						shapes[sps].mesh.positions[1 + 3 * shapes[0].mesh.indices[idx + 0]],
						shapes[sps].mesh.positions[2 + 3 * shapes[0].mesh.indices[idx + 0]]);

					L._O = Vec3(
						shapes[sps].mesh.positions[0 + 3 * shapes[0].mesh.indices[idx + 0]],
						shapes[sps].mesh.positions[1 + 3 * shapes[0].mesh.indices[idx + 0]],
						shapes[sps].mesh.positions[2 + 3 * shapes[0].mesh.indices[idx + 0]]);

					L._Area = L._U.cross(L._V).norm() * 0.5f;
					m_SumLightArea += L._Area;
					L._U.normalize();
					L._V.normalize();
					L._Normal = 0.333333 * (facet.GetNormal(0) + facet.GetNormal(1) + facet.GetNormal(2));
					L._Normal.normalize();
					L._Power = 600.0f;
					m_LightMeshs.push_back(L);

				}

				Facets.push_back(facet);

			}
			else{

				Facet facet(
					Vertex(
					Vec3(
					shapes[sps].mesh.positions[0 + 3 * shapes[0].mesh.indices[idx + 0]],
					shapes[sps].mesh.positions[1 + 3 * shapes[0].mesh.indices[idx + 0]],
					shapes[sps].mesh.positions[2 + 3 * shapes[0].mesh.indices[idx + 0]]),

					Vec3(
					shapes[sps].mesh.normals[0 + 3 * shapes[0].mesh.indices[idx + 0]],
					shapes[sps].mesh.normals[1 + 3 * shapes[0].mesh.indices[idx + 0]],
					shapes[sps].mesh.normals[2 + 3 * shapes[0].mesh.indices[idx + 0]]),

					Vec2(0, 0)
					//shapes[sps].mesh.texcoords[0 + 2 * shapes[0].mesh.indices[idx + 0]],
					//shapes[sps].mesh.texcoords[1 + 2 * shapes[0].mesh.indices[idx + 0]])
					),

					Vertex(
					Vec3(
					shapes[sps].mesh.positions[0 + 3 * shapes[0].mesh.indices[idx + 1]],
					shapes[sps].mesh.positions[1 + 3 * shapes[0].mesh.indices[idx + 1]],
					shapes[sps].mesh.positions[2 + 3 * shapes[0].mesh.indices[idx + 1]]),

					Vec3(
					shapes[sps].mesh.normals[0 + 3 * shapes[0].mesh.indices[idx + 1]],
					shapes[sps].mesh.normals[1 + 3 * shapes[0].mesh.indices[idx + 1]],
					shapes[sps].mesh.normals[2 + 3 * shapes[0].mesh.indices[idx + 1]]),

					Vec2(0, 0)
					//shapes[sps].mesh.texcoords[0 + 2 * shapes[0].mesh.indices[idx + 1]],
					//shapes[sps].mesh.texcoords[1 + 2 * shapes[0].mesh.indices[idx + 1]])
					),

					Vertex(
					Vec3(
					shapes[sps].mesh.positions[0 + 3 * shapes[0].mesh.indices[idx + 2]],
					shapes[sps].mesh.positions[1 + 3 * shapes[0].mesh.indices[idx + 2]],
					shapes[sps].mesh.positions[2 + 3 * shapes[0].mesh.indices[idx + 2]]),

					Vec3(
					shapes[sps].mesh.normals[0 + 3 * shapes[0].mesh.indices[idx + 2]],
					shapes[sps].mesh.normals[1 + 3 * shapes[0].mesh.indices[idx + 2]],
					shapes[sps].mesh.normals[2 + 3 * shapes[0].mesh.indices[idx + 2]]),

					Vec2(0, 0)
					//shapes[sps].mesh.texcoords[0 + 2 * shapes[0].mesh.indices[idx + 2]],
					//shapes[sps].mesh.texcoords[1 + 2 * shapes[0].mesh.indices[idx + 2]])
					),

					Vec3(1, 1, 1), Facet::ATT_MATERIAL::DIFFUSE
					);

				Facets.push_back(facet);
			}
			
			face++;
		}
	}

	
	time_t ed = time(NULL);
	std::cout << "Load OBJ END " << (ed - stt) << "sec" << std::endl;
	stt = time(NULL);

	std::cout << "Create " << m_SceneData->Print() << " START" << std::endl;
	if (m_SceneData){
		m_SceneData->CreateGraph(Facets);
	}
	ed = time(NULL);
	std::cout << "Create " << m_SceneData->Print() << " END " << (ed - stt) << "sec" << std::endl;

	std::cout << "Create Scene END" << std::endl;

	std::vector<Facet>().swap(Facets);

	return true;

}

float Scene::GetSceneWidth() const{
	return m_SceneData->GetSceneWidth();
}
Vec3 Scene::GetSceneCenter() const{
	return m_SceneData->GetSceneCenter();
}

void Scene::WriteBVH(std::ostream& ost) const{
	std::cout << "Write OBJ SATRT" << std::endl;
	m_SceneData->WriteObj(ost);
	std::cout << "Write OBJ END" << std::endl;
	m_SceneData->PrintGraph(std::cout);
}


int Scene::GetFacetFromRay(const PTUtility::Ray& ray, std::vector<AbsG::RF>& Facets) const{
	return m_SceneData->GetElements(ray, Facets);
}

const Vec3 Scene::GetTextureColor(const std::string& TextureName, const Vec2& UV) const{
	
	
	if (m_TextureMap.find(TextureName) == m_TextureMap.end()){
		return Vec3(-1, -1, -1);
	}
	else{
		

		std::shared_ptr<Texture2D> tex = m_TextureMap.at(TextureName);
		return tex->Sample(UV.x(), UV.y());
	}

}


Scene::Scene(AbsG* DataStructure) : m_SceneData(DataStructure){
	m_SumLightArea = 0.0f;
}
Scene::~Scene(){
	delete m_SceneData;
}