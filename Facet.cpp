#include"Facet.h"


using namespace PTUtility;

unsigned long int Facet::id = 0;

//三角形ポリゴン
class Facet::Triangle{

private:

	const void Fc(){
		return;
	}

public:

	const Vec3& GetCenter() const{
		return m_FCenter;
	}
	const Vec3& GetFaceNormal() const{
		return m_FNormal;
	}
	const Vec3& GetPos(int idx) const{
		return m_vertex[idx].m_Pos;
	}
	const Vec3& GetNormal(int idx) const{
		return m_vertex[idx].m_Normal;
	}
	const Vec2& GetUV(int idx) const{
		return m_vertex[idx].m_UV;
	}


	Vec3 m_FNormal;
	Vec3 m_FCenter;
	Vertex m_vertex[3];

	Triangle(){
		m_FNormal = Vec3(0, 1, 0);
		m_FCenter = Vec3(0, 0, 0);
		Fc();
	}

	//うーん
	Triangle(const Vertex& v1, const Vertex& v2, const Vertex& v3){
		m_vertex[0] = v1;
		m_vertex[1] = v2;
		m_vertex[2] = v3;

		m_FNormal = 0.33333 * (v1.m_Normal + v2.m_Normal + v3.m_Normal);
		m_FCenter = 0.33333 * (v1.m_Pos + v2.m_Pos + v3.m_Pos);

		Fc();
	}

	virtual ~Triangle(){
	}


};


//真っ赤に染めます
Facet::Facet(
	const Vertex& v1,
	const Vertex& v2,
	const Vertex& v3)
	: m_triangle(new Triangle(v1, v2, v3)), m_attribute(new AttributeTemple(Vec3(1, 0, 0), ATT_MATERIAL::DIFFUSE)), m_ID(id){
	id++;
}

Facet::Facet(
	const Vertex& v1,
	const Vertex& v2,
	const Vertex& v3,
	const Vec3& color)
	: m_triangle(new Triangle(v1, v2, v3)), m_attribute(new AttributeTemple(color, ATT_MATERIAL::DIFFUSE)), m_ID(id){
	id++;
}

Facet::Facet(
	const Vertex& v1,
	const Vertex& v2,
	const Vertex& v3,
	const Vec3& color,
	AttributeBSDF::BF BRDF,
	AttributeBSDF::BF BTDF,
	const std::string& TextureName)
	: m_triangle(new Triangle(v1, v2, v3)), m_attribute(new AttributeBSDF(color, BRDF, BTDF, TextureName)), m_ID(id){
	id++;
}
Facet::Facet(
	const Vertex& v1,
	const Vertex& v2,
	const Vertex& v3,
	const Vec3& color,
	float Roughness,
	float Metal,
	float Trans,
	const std::string& TextureName)
	: m_triangle(new Triangle(v1, v2, v3)), m_attribute(new AttributeParam(color, Roughness, Metal, Trans, TextureName)), m_ID(id){
	id++;
}
Facet::Facet(
	const Vertex& v1,
	const Vertex& v2,
	const Vertex& v3,
	const Vec3& color,
	ATT_MATERIAL mat,
	const std::string& TextureName)
	: m_triangle(new Triangle(v1, v2, v3)), m_attribute(new AttributeTemple(color, mat, TextureName)), m_ID(id){
	id++;
}


Facet::Facet(
	const Vertex& v1,
	const Vertex& v2,
	const Vertex& v3,
	const Vec3& color,
	AttributeBSDF::BF BRDF,
	AttributeBSDF::BF BTDF)
	: m_triangle(new Triangle(v1, v2, v3)), m_attribute(new AttributeBSDF(color, BRDF, BTDF)), m_ID(id){
	id++;
}
Facet::Facet(
	const Vertex& v1,
	const Vertex& v2,
	const Vertex& v3,
	const Vec3& color,
	float Roughness,
	float Metal,
	float Trans)
	: m_triangle(new Triangle(v1, v2, v3)), m_attribute(new AttributeParam(color, Roughness, Metal, Trans)), m_ID(id){
	id++;
}
Facet::Facet(
	const Vertex& v1,
	const Vertex& v2,
	const Vertex& v3,
	const Vec3& color,
	ATT_MATERIAL mat)
	: m_triangle(new Triangle(v1, v2, v3)), m_attribute(new AttributeTemple(color, mat)), m_ID(id){
	id++;
}



Facet::~Facet(){
	//たぶんなにもしなくておk
}

const Vec3& Facet::GetCenter()const{
	return m_triangle->GetCenter();
}
float Facet::GetP(const Vec3& in, const Vec3& out)const{
	return m_attribute->GetP(in, out);
}
const Vec3& Facet::GetPos(int idx)const{
	return m_triangle->GetPos(idx);
}
const Vec3& Facet::GetNormal(int idx)const{
	return m_triangle->GetNormal(idx);
}
const Vec2& Facet::GetUV(int idx)const{
	return m_triangle->GetUV(idx);
}

Vec3 Facet::GetEmission(const Vec3& out)const{ return m_attribute->m_Light(out); }
Vec3 Facet::GetEmission()const{ return m_attribute->m_emit * m_attribute->m_Color; }









