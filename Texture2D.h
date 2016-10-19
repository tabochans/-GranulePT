#pragma once
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.hpp"
#include"Vec3.h"

namespace PTUtility{

	class Texture2D{
	private:
		int m_X;
		int m_Y;
		int m_Comp;
		Vec3 m_Color;

		unsigned char* m_data;

		const Texture2D& operator=(const Texture2D&);
		Texture2D(const Texture2D&);

		bool Load(const char* FileName);

	protected:

	public:
		Texture2D();
		explicit Texture2D(const char* FileName);
		explicit Texture2D(const Vec3& color);
		virtual ~Texture2D();

		Vec3 Sample(float u, float v)const;
		Vec3 GetColor(int x, int y)const;

	};

}