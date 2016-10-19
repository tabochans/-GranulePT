#pragma once
#include"iIBL.h"
#include"Texture2D.h"

namespace PTUtility{

	class IBLTexture : public iIBL{

	private:
		Texture2D m_texture;

	public:

		virtual Vec3 GetColor1(const Vec3& Dir) const;
		virtual Vec3 GetColor2(const Vec3& Dir) const;

		IBLTexture(const char* FileName);
		IBLTexture(const Vec3& Color);
		IBLTexture();

		virtual ~IBLTexture();


	};



	class IBLDLight : public iIBL{

	private:

		Vec3 _Dir;
		float _S;
		Vec3 _Color;

	public:

		virtual Vec3 GetColor1(const Vec3& Dir) const;
		virtual Vec3 GetColor2(const Vec3& Dir) const;

		IBLDLight(const Vec3& Direction, float S, const Vec3& Color);

		virtual ~IBLDLight();
	};


}