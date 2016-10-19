#pragma once
#include"Vec3.h"
#include<algorithm>

namespace PTUtility{

	class iIBL{

	private:

	public:

		virtual Vec3 GetColor1(const Vec3& Dir) const = 0;
		virtual Vec3 GetColor2(const Vec3& Dir) const = 0;

		iIBL(){}
		virtual ~iIBL(){}


	};


}