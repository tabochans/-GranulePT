#include"Intersect.h"
#include"Eigen\Core"
#include"Eigen\Geometry"

using namespace PTUtility;

bool Intersect::Sphere_AABB_Fast_S(
	const __m256* CPos,
	const __m256* R,
	const __m256 value[2][3],
	int& mask
	){

	__m256 fg[3];
	__m256 zero = _mm256_set1_ps(0.0f);

	fg[0] = _mm256_min_ps(
		_mm256_sub_ps(_mm256_add_ps(CPos[0], R[0]), value[0][0]),
		_mm256_sub_ps(value[1][0], _mm256_sub_ps(CPos[0], R[0]))
		);
	fg[1] = _mm256_min_ps(
		_mm256_sub_ps(_mm256_add_ps(CPos[1], R[0]), value[0][1]),
		_mm256_sub_ps(value[1][1], _mm256_sub_ps(CPos[1], R[0]))
		);
	fg[2] = _mm256_min_ps(
		_mm256_sub_ps(_mm256_add_ps(CPos[2], R[0]), value[0][2]),
		_mm256_sub_ps(value[1][2], _mm256_sub_ps(CPos[2], R[0]))
		);

	__m256 mfg = _mm256_min_ps(fg[0], _mm256_min_ps(fg[1], fg[2]));


	mask = _mm256_movemask_ps(_mm256_cmp_ps(mfg, zero, _CMP_GT_OS));
	return (mask > 0);

}


bool Intersect::Sphere_Triangle(const Vec3& CPos, float R, const Facet& facet, float& u, float& v, float& t){
	
	Vec3 N = (0.3333 * (facet.GetNormal(0) + facet.GetNormal(1) + facet.GetNormal(2))).normalized();
	Vec3 p = facet.GetPos(1) - facet.GetPos(0);
	Vec3 q = facet.GetPos(2) - facet.GetPos(0);

	Eigen::Matrix3f mat;

	mat << -N.x(), -N.y(), -N.z(), p.x(), p.y(), p.z(), q.x(), q.y(), q.z();

	mat.inverse();
	Eigen::Vector3f pq = Eigen::Vector3f(
		CPos[0] - facet.GetPos(0)[0],
		CPos[1] - facet.GetPos(0)[1],
		CPos[2] - facet.GetPos(0)[2]);

	Eigen::Vector3f x = mat * pq;

	if (x[0] <= R && x[1] >= 0.0f && x[2] >= 0.0f && x[1] + x[2] <= 1.0f){
		return true;
	}
	else if ((CPos - facet.GetPos(2)).norm() <= R){
		return true;
	}
	else if ((CPos - facet.GetPos(1)).norm() <= R){
		return true;
	}
	else if ((CPos - facet.GetPos(0)).norm() <= R){
		return true;
	}
	else if (x[1] < 0.0f){
		if(Intersect::Ray_Sphere(Ray(facet.GetPos(0), p), CPos, R, t)){
			return true;
		}
	}
	else if (x[2] < 0.0f){
		if (Intersect::Ray_Sphere(Ray(facet.GetPos(0), q), CPos, R, t)){
			return true;
		}
	}
	else if (x[1] + x[2] < 1.0f){
		if (Intersect::Ray_Sphere(Ray(facet.GetPos(1), facet.GetPos(2) - facet.GetPos(1)), CPos, R, t)){
			return true;
		}
	}

	/*else if (x[1] <= 0.0f && x[2] >= 1.0f && abs(x[1]) + abs(x[2]) >= 1.0f && (CPos - facet.GetPos(2)).norm() <= R){
		return true;
	}
	else if (x[1] >= 1.0f && x[2] <= 0.0f && abs(x[1]) + abs(x[2]) >= 1.0f && (CPos - facet.GetPos(1)).norm() <= R){
		return true;
	}
	else if (x[1] <= 0.0f && x[2] <= 0.0f && abs(x[1]) + abs(x[2]) <= 1.0f && (CPos - facet.GetPos(0)).norm() <= R){
		return true;
	}*/


	return false;
}

bool Intersect::Voxel_Ray(const Vec3& Vpos, float width, const Ray& ray, float& tm, float& TM){

	tm = 100000;
	TM = -10000;

	float tmax, tmin;
	const Vec3 max = Vpos + Vec3::Ones() * 0.5 * width;
	const Vec3 min = Vpos - Vec3::Ones() * 0.5 * width;

	if (abs(ray.m_Dir.x()) < 0.0001){
		tmax = (2.0f * (ray.m_Dir.x() > 0) - 1) * 10000.0f * (max.x() - ray.m_Org.x());
		tmin = (2.0f * (ray.m_Dir.x() > 0) - 1) * 10000.0f * (min.x() - ray.m_Org.x());
	}
	else{
		tmax = (max.x() - ray.m_Org.x()) / ray.m_Dir.x();
		tmin = (min.x() - ray.m_Org.x()) / ray.m_Dir.x();
	}
	if (tmax < tmin){
		swap(tmax, tmin);
	}

	float tymin, tymax;
	if (abs(ray.m_Dir.y()) < 0.0001){
		tymax = (2.0f * (ray.m_Dir.y() > 0) - 1) * 10000.0f * (max.y() - ray.m_Org.y());
		tymin = (2.0f * (ray.m_Dir.y() > 0) - 1) * 10000.0f * (min.y() - ray.m_Org.y());
	}
	else{
		tymax = (max.y() - ray.m_Org.y()) / ray.m_Dir.y();
		tymin = (min.y() - ray.m_Org.y()) / ray.m_Dir.y();
	}
	if (tymax < tymin){
		swap(tymax, tymin);
	}

	if (tmin > tymax || tymin > tmax){
		return false;
	}
	tmin = tmin > tymin ? tmin : tymin;
	tmax = tmax < tymax ? tmax : tymax;

	if (abs(ray.m_Dir.z()) < 0.0001){
		tymax = (2.0f * (ray.m_Dir.z() > 0) - 1) * 10000.0f * (max.z() - ray.m_Org.z());
		tymin = (2.0f * (ray.m_Dir.z() > 0) - 1) * 10000.0f * (min.z() - ray.m_Org.z());
	}
	else{
		tymax = (max.z() - ray.m_Org.z()) / ray.m_Dir.z();
		tymin = (min.z() - ray.m_Org.z()) / ray.m_Dir.z();
	}
	if (tymax < tymin){
		swap(tymax, tymin);
	}

	if (tmin > tymax || tymin > tmax){
		return false;
	}
	tmin = tmin > tymin ? tmin : tymin;
	tmax = tmax < tymax ? tmax : tymax;

	tm = tmin;
	TM = tmax;

	//if (tmin < FLT_MIN){ return false; }

	return true;




}


bool Intersect::IsHit(
	const __m256 value[2][3], 
	const __m256* org,        
	const __m256* idir,       
	const int* sign,          
	int& mask
	)
{

	__m256 tmin = _mm256_set1_ps(-100000.0f);
	__m256 tmax = _mm256_set1_ps( 100000.0f);


	int idx0, idx1;


	// YŽ².
	idx0 = sign[1];
	idx1 = 1 - idx0;
	tmin = _mm256_max_ps(tmin, _mm256_mul_ps(_mm256_sub_ps(value[idx0][1], org[1]), idir[1]));
	tmax = _mm256_min_ps(tmax, _mm256_mul_ps(_mm256_sub_ps(value[idx1][1], org[1]), idir[1]));


	// ZŽ².
	idx0 = sign[2];
	idx1 = 1 - idx0;
	tmin = _mm256_max_ps(tmin, _mm256_mul_ps(_mm256_sub_ps(value[idx0][2], org[2]), idir[2]));
	tmax = _mm256_min_ps(tmax, _mm256_mul_ps(_mm256_sub_ps(value[idx1][2], org[2]), idir[2]));


	// XŽ².
	idx0 = sign[0];
	idx1 = 1 - idx0;
	tmin = _mm256_max_ps(tmin, _mm256_mul_ps(_mm256_sub_ps(value[idx0][0], org[0]), idir[0]));
	tmax = _mm256_min_ps(tmax, _mm256_mul_ps(_mm256_sub_ps(value[idx1][0], org[0]), idir[0]));


	mask = _mm256_movemask_ps(_mm256_cmp_ps(tmax, tmin, _CMP_GT_OS));
	return (mask > 0);

}

void Intersect::swap(float& right, float& left){
	float buf = right;
	right = left;
	left = buf;
}

bool Intersect::AABB_Ray(const AABB& aabb, const Ray& ray){

	float tmax, tmin;
	const Vec3& max = aabb.GetMaxPos();
	const Vec3& min = aabb.GetMinPos();

	if (abs(ray.m_Dir.x()) < 0.0001){
		tmax = (2.0f * (ray.m_Dir.x() > 0) - 1) * 10000.0f * (max.x() - ray.m_Org.x());
		tmin = (2.0f * (ray.m_Dir.x() > 0) - 1) * 10000.0f * (min.x() - ray.m_Org.x());
	}
	else{
		tmax = (max.x() - ray.m_Org.x()) / ray.m_Dir.x();
		tmin = (min.x() - ray.m_Org.x()) / ray.m_Dir.x();
	}
	if (tmax < tmin){
		swap(tmax, tmin);
	}

	float tymin, tymax;
	if (abs(ray.m_Dir.y()) < 0.0001){
		tymax = (2.0f * (ray.m_Dir.y() > 0) - 1) * 10000.0f * (max.y() - ray.m_Org.y());
		tymin = (2.0f * (ray.m_Dir.y() > 0) - 1) * 10000.0f * (min.y() - ray.m_Org.y());
	}
	else{
		tymax = (max.y() - ray.m_Org.y()) / ray.m_Dir.y();
		tymin = (min.y() - ray.m_Org.y()) / ray.m_Dir.y();
	}
	if (tymax < tymin){
		swap(tymax, tymin);
	}

	if (tmin > tymax || tymin > tmax){
		return false;
	}
	tmin = tmin > tymin ? tmin : tymin;
	tmax = tmax < tymax ? tmax : tymax;

	if (abs(ray.m_Dir.z()) < 0.0001){
		tymax = (2.0f * (ray.m_Dir.z() > 0) - 1) * 10000.0f * (max.z() - ray.m_Org.z());
		tymin = (2.0f * (ray.m_Dir.z() > 0) - 1) * 10000.0f * (min.z() - ray.m_Org.z());
	}
	else{
		tymax = (max.z() - ray.m_Org.z()) / ray.m_Dir.z();
		tymin = (min.z() - ray.m_Org.z()) / ray.m_Dir.z();
	}
	if (tymax < tymin){
		swap(tymax, tymin);
	}

	if (tmin > tymax || tymin > tmax){
		return false;
	}
	tmin = tmin > tymin ? tmin : tymin;
	tmax = tmax < tymax ? tmax : tymax;

	//if (tmin < FLT_MIN){ return false; }

	return true;

}
bool Intersect::Ray_Triangle(const Ray& ray, const Facet& facet, float& u, float& v ,float& t){

	const Vec3& d = ray.m_Dir;
	const Vec3& e1 = (facet.GetPos(1) - facet.GetPos(0));
	const Vec3& e2 = (facet.GetPos(2) - facet.GetPos(0));
	const Vec3& p = d.cross(e2);

	float a = e1.dot(p);
	if (a > -FLT_MIN && a < FLT_MIN) { u = v = t = 0.0; return false; }
	float f = 1.0f / a;
	Vec3 s = ray.m_Org - facet.GetPos(0);
	u = f * s.dot(p);
	if (u < -0.01 || u > 1.0 + 0.01) { u = v = t = 0.0; return false; }

	Vec3 q = s.cross(e1);
	v = f * d.dot(q);
	if (v < -0.01 || u + v > 1.0 + 0.01) { u = v = t = 0.0; return false; }

	t = f * e2.dot(q);

	return true;

}
bool Intersect::Ray_Sphere(
	const Ray& ray, const Vec3& CPos, float R, float& t){
	t = 10000;


	if ((ray.m_Org - CPos).norm2() < R * R){
		//‚»‚à‚»‚àRay‚ª“à•”

		return false;
	}
	else{

		Vec3 CO = (CPos - ray.m_Org);

		float Len = CO.norm();


		float cos = CO.dot(ray.m_Dir) / Len;
		if (cos < 0.0f){
			return false;
		}
		if (Len * Len * (1.0f - cos*cos) < R * R){
			t = Len * cos - sqrt(R*R - Len*Len *(1.0f - cos*cos));
			return true;
		}

	}

	
	


	return false;
}


bool Intersect::Ray_Sphere_In(
	const Ray& ray, const Vec3& CPos, float R, float& t){
	t = 10000;

	Vec3 CO = (CPos - ray.m_Org);

	if (CO.norm2() < R * R){
		//‚»‚à‚»‚àRay‚ª“à•”

		float Len = CO.norm();

		float cos = ray.m_Dir.dot(CO) / Len;

		t = 0.5* (-Len*cos + sqrt(Len*Len*cos*cos - 4.0f*(Len*Len - R*R)));


		return true;
	}
	else{

		float Len = CO.norm();

		float cos = CO.dot(ray.m_Dir) / Len;

		if (cos < 0.0f){
			return false;
		}
		if (Len * Len * (1.0f - cos*cos) < R * R){
			t = Len * cos - sqrt(R*R - Len*Len *(1.0f - cos*cos));
			return true;
		}

	}





	return false;
}