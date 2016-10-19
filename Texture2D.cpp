#include"Texture2D.h"


using namespace PTUtility;

Texture2D::Texture2D(){
	m_Color = Vec3(1, 1, 1);
	m_data = NULL;
}
Texture2D::Texture2D(const char* FileName){
	m_data = NULL;
	Load(FileName);
}
Texture2D::Texture2D(const Vec3& color){
	m_data = NULL;
	m_Color = color;
}
Texture2D::~Texture2D(){
	stbi_image_free(m_data);
}

bool Texture2D::Load(const char* FileName){
	if (FileName == ""){
		return false;
	}
	if (FileName == NULL){
		return false;
	}

	if (m_data){
		stbi_image_free(m_data);
	}
	m_data = NULL;

	m_data = stbi_load(FileName, &m_X, &m_Y, &m_Comp, 0);

	if (m_data == NULL){
		return false;
	}
	m_Color = Vec3(1, 1, 1);
	
	return true;

}

Vec3 Texture2D::Sample(float u, float v)const{
	int x, y;
	x = (int(abs(u) * m_X)) % m_X;
	y = (int(abs(v) * m_Y)) % m_Y;
	return GetColor(x, y);
}
Vec3 Texture2D::GetColor(int x, int y)const{
	
	float r, g, b;
	if (m_data == NULL){
		r = g = b = 0.0f;
	}
	else{
		r = m_data[m_Comp*(x + y * m_X) + 0] / 255.01f;
		g = m_data[m_Comp*(x + y * m_X) + 1] / 255.01f;
		b = m_data[m_Comp*(x + y * m_X) + 2] / 255.01f;
	}
	return Vec3(r, g, b);
}

