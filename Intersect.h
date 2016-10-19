#pragma once
#include"Vec3.h"
#include"AABB.h"
#include"AABB_S.h"
#include"Ray.h"
#include"Facet.h"
#include <immintrin.h>

//åç∑îªíËÇÇ™ÇÒÇŒÇÈêl

namespace PTUtility{

	namespace Intersect{
		void swap(float& right, float& left);
		bool Voxel_Ray(const Vec3& Vpos, float width, const Ray& ray, float& tm, float& TM);
		bool AABB_Ray(const AABB& aabb, const Ray& ray);
		bool Ray_Triangle(const Ray& ray, const Facet& facet, float& u, float& v, float& t);
		bool Ray_Sphere(
			const Ray& ray, const Vec3& CPos, float R, float& t);

		bool Ray_Sphere_In(
			const Ray& ray, const Vec3& CPos, float R, float& t);

		bool IsHit(
			const __m256 value[2][3],  
			const __m256* org,        
			const __m256* idir,      
			const int* sign, 
			int& mask
			);


		bool Sphere_AABB_Fast_S(
			const __m256* CPos,
			const __m256* R,
			const __m256 value[2][3],
			int& mask);

		bool Sphere_Triangle(const Vec3& CPos, float R, const Facet& facet, float& u, float& v, float& t);
	}

}