#include"IBL.h"

using namespace PTUtility;


IBLTexture::IBLTexture(const char* FileName) : m_texture(FileName){

}

IBLTexture::IBLTexture(const Vec3& Color) : m_texture(Vec3(Color)){

}

IBLTexture::IBLTexture() : m_texture(Vec3(0, 0, 0)){

}


IBLTexture::~IBLTexture(){

}

Vec3 IBLTexture::GetColor1(const Vec3& dir)const{

	//とりあえずテキトーに
	float u, v;
	float x, y, z;
	float phc, phs, r;
	x = dir.x();
	y = dir.y();
	z = dir.z();

	r = 0.5f - 0.4 * abs(z);

	if (x*x + z*z > 0.01){
		phc = x / sqrt(x*x + z*z);
		phs = y / sqrt(x*x + z*z);
		return m_texture.Sample(r*phc + 0.5f, r*phs + 0.5f);
	}
	else{
		return m_texture.Sample(0, 0);
	}

	
}


Vec3 IBLTexture::GetColor2(const Vec3& dir) const{

	//とりあえずテキトーに
	float u, v;
	float x, y, z;
	float ph, r;
	x = dir.x();
	y = dir.y();
	z = dir.z();

	if (x*x + z* z > 0.01){
		ph = (acos(x / sqrt(x*x + z*z)));
	}
	else{
		ph = 0.0f;
	}

	return 1.5f * m_texture.Sample(ph / (3.1416), 0.6 + (-y + 1.0f));
}









IBLDLight::IBLDLight(const Vec3& Direction, float S, const Vec3& Color) : _Dir(-Direction), _S(S), _Color(Color){
}


IBLDLight::~IBLDLight(){

}

Vec3 IBLDLight::GetColor1(const Vec3& dir)const{
	
	float x, y, z;
	x = dir.x();
	y = dir.y();
	z = dir.z();

	float LL = std::max(0.0f, std::min(dir.dot(_Dir), 0.999f));

	return 10.0f * _Color * pow(LL, _S);
	return 5.0f * _Color * pow(LL, _S);

}


Vec3 IBLDLight::GetColor2(const Vec3& dir) const{

	float x, y, z;
	x = dir.x();
	y = dir.y();
	z = dir.z();

	float LL = std::max(0.0f, std::min(dir.dot(_Dir), 0.999f));
	float RL = pow(LL, _S);

	if (RL > 0.875){
		return 10.0f * _Color;
	}
	else{
		return Vec3::Zero();
	}
}


