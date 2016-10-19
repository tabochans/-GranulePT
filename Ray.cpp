#include"Ray.h"

using namespace PTUtility;

Vec3 Ray::ReflectRay(const Ray& in, const Vec3& normal){

	return (in.m_Dir - 2.0f * (normal.dot(in.m_Dir)) * normal).normalized();

}

//in �́A�����Ă�������@�s�����ގ��Ȃ�A�@���Ƌt�����Ŏw�肷��
Vec3 Ray::RefractRay(const Ray& in, const Vec3& normal, float eta_from, float eta_to){
	return normal;
}

//in �́A�����Ă�������@�s�����ގ��Ȃ�A�@���Ƌt�����Ŏw�肷��
float Ray::FresnelReflectance(const Ray& in, const Vec3& normal){
	return 1.0f;
}

//in �́A�����Ă�������@�s�����ގ��Ȃ�A�@���Ƌt�����Ŏw�肷��
float Ray::FresnelTransmittance(const Ray& in, const Vec3& normal, float eta_from, float eta_to){
	return 1.0f;
}


Ray::Ray(const Vec3& Origin, const Vec3& Direction) : m_Org(Origin), m_Dir(Direction.normalized()), m_Index(1.0f){
}

Ray::Ray(const Vec3& Origin, const Vec3& Direction, float Index) : m_Org(Origin), m_Dir(Direction.normalized()), m_Index(Index){

}

Ray::~Ray(){

}