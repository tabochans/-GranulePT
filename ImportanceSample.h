#pragma once

#include"Vec3.h"
#include"Ray.h"
#include"Facet.h"
#include"MT.h"

//インポータンスサンプルをしてくれる？
namespace PTUtility{

	namespace ISample{

		float NextDirDiffuse(const Ray& in, const Vec3& Normal, Vec3& result, RandomMT& MT);
		float NextDirSpecular(const Ray& in, const Vec3& Normal, Vec3& result, RandomMT& MT);
		float NextDirRefraction(const Ray& in, const Vec3& Normal, float& n_from, float& n_to, Vec3& result, RandomMT& MT);
		float NextDirSphere(const Ray& in, const Vec3& Normal, Vec3& result, RandomMT& MT);
		float NextDirHemSphere(const Ray& in, const Vec3& Normal, Vec3& result, RandomMT& MT);
		float NextDirPhong(const Ray& in, const Vec3& Normal, Vec3& result, RandomMT& MT, float s);

	}


}