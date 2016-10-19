#include"ImportanceSample.h"
#define M_PI 3.14159265358979323846

using namespace PTUtility;

float ISample::NextDirPhong(const Ray& in, const Vec3& Normal, Vec3& result, RandomMT& MT, float s){

	Vec3 e1, e2, e3;
	
	e1 = Ray::ReflectRay(in, Normal);
	if (abs(e1.x()) < 0.00001){
		e2 = Vec3(1, 0, 0).cross(e1).normalized();
	}
	else{
		e2 = Vec3(0, 1, 0).cross(e1).normalized();
	}
	e3 = e1.cross(e2).normalized();

	float ph = MT.genrand64_real1() * 2.0f * M_PI;
	float th = MT.genrand64_real1();
	float cth2 = pow(1-th,1.0f / (s + 1.0f));

	float cph = cos(ph);
	float sph = sin(ph);
	float sth2 = sqrt(1.0001f - cth2*cth2);


	result = Vec3(e3* sth2 * cph + e2*sth2 * sph + e1*cth2);
	result.normalize();

	return 1.0f;
}


float ISample::NextDirSphere(const Ray& in, const Vec3& Normal, Vec3& result, RandomMT& MT){

	float ph = MT.genrand64_real1() * 2.0f * M_PI;
	float th = MT.genrand64_real1() * M_PI;
	th = acos(MT.genrand64_real1() * 2.0f - 1.0f);

	result = Vec3(sinf(th)*cosf(ph), sinf(th)*sinf(ph), cosf(th));
	result.normalize();
	
	//‚È‚ñ‚Æ‚È‚­
	return 1.0f;

}



float ISample::NextDirDiffuse(
	const Ray& in, 
	const Vec3& Normal, 
	Vec3& result,
	RandomMT& MT){

	Vec3 e1, e2, e3;
	e1 = Normal;

	if (abs(e1.x()) > 0.000001){
		e2 = Vec3(0, 1, 0).cross(e1).normalized();
	}
	else{
		e2 = Vec3(1, 0, 0).cross(e1).normalized();
	}
	e3 = e1.cross(e2).normalized();

	float ph = MT.genrand64_real1() * 2.0f * M_PI;
	float th = MT.genrand64_real1();
	float th2 = sqrt(th);

	result = e1*sqrt(1.000000001 - th) + e2*cos(ph)*th2 + e3*sin(ph)*th2;
	result.normalize();

	//‚È‚ñ‚Æ‚È‚­
	return 1.0f;

}

float ISample::NextDirHemSphere(const Ray& in, const Vec3& Normal, Vec3& result, RandomMT& MT){

	Vec3 e1, e2, e3;
	e1 = Normal;
	if (abs(e1.x()) < 0.00001){
		e2 = Vec3(1, 0, 0).cross(e1).normalized();
	}
	else{
		e2 = Vec3(0, 1, 0).cross(e1).normalized();
	}
	e3 = e1.cross(e2).normalized();

	float ph = MT.genrand64_real1() * 2.0f * M_PI;
	float th = MT.genrand64_real1() * M_PI * 0.5f;
	th = acos(MT.genrand64_real1());


	result = Vec3(e1* cosf(th) + e2*cos(ph)*sinf(th) + e3*sin(ph)*sinf(th));
	result.normalize();

	//‚È‚ñ‚Æ‚È‚­
	return 1.0f;

}



float ISample::NextDirSpecular(const Ray& in, const Vec3& Normal, Vec3& result, RandomMT& MT){

	result = Ray::ReflectRay(in, Normal);
	return 1.0f;
}


float ISample::NextDirRefraction(const Ray& in, const Vec3& Normal, float& n_from, float& n_to, Vec3& result, RandomMT& MT){

	float p = 0.0f;
	float sn = (n_from) / (n_to);
	float ns = sn;
	if (sn > 1){ ns = (n_to) / (n_from); }
	bool into = in.m_Dir.dot(Normal) > 0;
	Vec3 NN = Normal;
	if (into){ NN = -Normal; }
	float fls = 0.0f;
	float tr = 0.0f;

	if(abs(n_from-n_to) < 0.01){
		result = in.m_Dir;
		tr = 1.0f;
		fls = 0.0f;
		return 1.0f;
	}
	if (n_from < n_to){
		result = -(1 + sn)*NN + sn*in.m_Dir;
		result.normalize();
		float R0 = pow((1.0f / ns - 1.0f) / (1.0f / ns + 1.0f), 2);
		float cs = into ? -in.m_Dir.dot(NN) : -result.dot(Normal);
		fls = R0 + (1.0f - R0)*pow(1.01 - cs, 5);
		float nn2 = pow(sn, 2);
		tr = (1 - fls)*nn2;
	}
	else if (n_from > n_to){

		if (in.m_Dir.dot(-NN) < 1.0f - pow(sn, 2) + 0.01){
			result = Ray::ReflectRay(in, Normal);
			n_to = n_from;
			return 1.0f;
		}
		result = -(1 + sn)*NN + sn*in.m_Dir;
		result.normalize();
		float R0 = pow((1.0f / ns - 1.0f) / (1.0f / ns + 1.0f), 2);
		float cs = into ? -in.m_Dir.dot(NN) : -result.dot(Normal);
		fls = R0 + (1.0f - R0)*pow(1.01 - cs, 5);
		float nn2 = pow(sn, 2);
		tr = (1 - fls)*nn2;

	}

	float tt = fls + tr;
	float pr = fls / tt;
	float pt = tr / tt;

	if (MT.genrand64_real1() < pr){
		result = Ray::ReflectRay(in, Normal);
		return 1.0f;
		return fls / 0.15;
	}
	else{
		return 1.0f;
		return tr / 0.85f;
	}
	
	return p;
}