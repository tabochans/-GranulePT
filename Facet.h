#pragma once


#include"Vec3.h"
#include"Vec2.h"
#include<string>
#include <memory>

/////////////////////////////
/*
Vertex class
Facet class
の二つをご利用いただけます
*/
/////////////////////////////
namespace PTUtility{


	struct Vertex;
	class Facet;


	//頂点かなー
	struct Vertex{

	public:

		Vec3 m_Pos;
		Vec3 m_Normal;
		Vec2 m_UV;

		Vertex() : m_Pos(0, 0, 0), m_Normal(0, 1, 0), m_UV(0, 0){
		}
		Vertex(
			const Vec3& pos,
			const Vec3& normal,
			const Vec2& uv
			) : m_Pos(pos), m_Normal(normal), m_UV(uv){
		}

		virtual ~Vertex(){
		}

	};


	//ポリゴンのこと？
	class Facet{
	private:
		static unsigned long int id;

	public:

		static bool compXX(const Facet& left, const Facet& right){
			if (left.m_ID == right.m_ID){ return false; }
			return (10000000000.0f * left.GetCenter().x() - 10000000000.0f * right.GetCenter().x()) + ((int)left.m_ID - (int)right.m_ID) < 0;
		}
		static bool compYY(const Facet& left, const Facet& right){
			if (left.m_ID == right.m_ID){ return false; }
			return (10000000000.0f * left.GetCenter().y() - 10000000000.0f * right.GetCenter().y()) + ((int)left.m_ID - (int)right.m_ID) < 0;
		}
		static bool compZZ(const Facet& left, const Facet& right){
			if (left.m_ID == right.m_ID){ return false; }
			return (10000000000.0f * left.GetCenter().z() - 10000000000.0f * right.GetCenter().z()) + ((int)left.m_ID - (int)right.m_ID) < 0;
		}


		static bool compFacetXP(const Facet* left, const Facet* right){
			return (left->GetCenter().x() - right->GetCenter().x()) < 0;
		}
		static bool compFacetYP(const Facet* left, const Facet* right){
			return (left->GetCenter().y() - right->GetCenter().y()) < 0;
		}
		static bool compFacetZP(const Facet* left, const Facet* right){
			return (left->GetCenter().z() - right->GetCenter().z()) < 0;
		}


		static bool compX(const Facet& left, const Facet& right){
			return left.GetCenter().x() < right.GetCenter().x();
		}
		static bool compY(const Facet& left, const Facet& right){
			return left.GetCenter().y() < right.GetCenter().y();
		}
		static bool compZ(const Facet& left, const Facet& right){
			return left.GetCenter().z() < right.GetCenter().z();
		}
		//マテリアルのテンプレ
		enum ATT_MATERIAL{
			DIFFUSE = 0,
			METAL = 1,
			GLASS = 2
		};

	private:

		class Triangle;


		//マテリアルを表現する感じのクラス
		class Attribute{

		protected:
			int m_Mode;

		public:

			static Vec3 NE(const Vec3& out){
				return Vec3::Zero();
			}

			Vec3 m_Color;
			std::string m_TextureName;
			float m_emit;

			typedef Vec3(*LF)(const Vec3& out);
			LF m_Light;

			void SetLightFunction(LF LFF){
				m_Light = (LF)LFF;
			}
			void SetEmit(float e, const Vec3& Color){
				m_emit = e;
				m_Color = Color;
			}
			
			//どんだけ反射とかをするかを返す
			virtual float GetP(const Vec3& in, const Vec3& out) = 0;

			Attribute(const Vec3& color, const std::string& tname) : m_Color(color),m_TextureName(tname),m_Light(NE),m_emit(0){
				m_Mode = -1;
			}
			Attribute(const Vec3& color, LF EmissionF) : m_Color(color), m_TextureName(""), m_Light(EmissionF), m_emit(0){
				m_Mode = -1;
			}
			Attribute(const Vec3& color) : m_Color(color), m_TextureName(""), m_Light(NE), m_emit(0){
				m_Mode = -1;
			}
			virtual ~Attribute(){

			}

		};

		////////////////////////
		//以下子孫たち　ここで宣言するのどうなんだろうね？
		////////////////////////

		//BSDFでマテリアルを表現する
		class AttributeBSDF : public Attribute{

		public:
			typedef float(*BF)(const Vec3& in, const Vec3& out);

		private:
			
			BF m_BRDF;
			BF m_BTDF;

		public:

			virtual float GetP(const Vec3& in, const Vec3& out){
				if (in.dot(out) < 0){
					return m_BTDF(in, out);
				}
				else{
					return m_BRDF(in, out);
				}
			}

			AttributeBSDF(
				const Vec3& color,
				BF BRDF,
				BF BTDF) : m_BRDF(BRDF), m_BTDF(BTDF), Attribute(color){
				m_Mode = 1;
			}
			AttributeBSDF(
				const Vec3& color,
				BF BRDF,
				BF BTDF,
				const std::string& TextureName) : m_BRDF(BRDF), m_BTDF(BTDF), Attribute(color,TextureName){
				m_Mode = 1;
			}
			virtual ~AttributeBSDF(){

			}

		};

		class AttributeTemple : public Attribute{

		public:

		private:

			ATT_MATERIAL m_Mat;

		public:

			virtual float GetP(const Vec3& in, const Vec3& out){
				switch (m_Mat)
				{
				case DIFFUSE:
					return 0.9f;
				case METAL:
					return 0.5f;
				case GLASS:
					return 0.3f;
				default:
					return 0.0f;
				}
			}

			AttributeTemple(
				const Vec3& color, ATT_MATERIAL mat) : m_Mat(mat), Attribute(color){
				m_Mode = 3;
			}
			AttributeTemple(
				const Vec3& color, ATT_MATERIAL mat,const std::string& TextureName) : m_Mat(mat), Attribute(color,TextureName){
				m_Mode = 3;
			}
			virtual ~AttributeTemple(){

			}

		};

		class AttributeParam : public Attribute{

		private:

			float m_Roughness;
			float m_metal;
			float m_trans;

		public:

			virtual float GetP(const Vec3& in, const Vec3& out){
				return 0.5f;
			}

			AttributeParam(
				const Vec3& color,
				float roughness,
				float metal,
				float trans) : m_Roughness(roughness), m_metal(metal), m_trans(trans), Attribute(color){
				m_Mode = 2;
			}
			AttributeParam(
				const Vec3& color,
				float roughness,
				float metal,
				float trans,
				const std::string& TextureName) : m_Roughness(roughness), m_metal(metal), m_trans(trans), Attribute(color,TextureName){
				m_Mode = 2;
			}
			virtual ~AttributeParam(){

			}


		};
		////////////////////////
		//以上子孫たち　ここで宣言するのどうなんだろうね？
		////////////////////////



		std::shared_ptr<const Triangle> m_triangle;
		std::shared_ptr<Attribute> m_attribute;


	protected:


	public:

		const std::string& GetTextureName()const{ return m_attribute->m_TextureName; }
		const Vec3& GetDiffuseColor()const{ return m_attribute->m_Color; }
		const Vec3& GetCenter() const;
		float GetP(const Vec3& in, const Vec3& out) const;
		const Vec3& GetPos(int idx) const;
		const Vec3& GetNormal(int idx) const;
		const Vec2& GetUV(int idx) const;

		Vec3 GetEmission(const Vec3& out)const;
		Vec3 GetEmission()const;
		unsigned long int m_ID;

		void SetLe(float power, const Vec3& color, Attribute::LF Le = Attribute::NE){
			m_attribute->SetEmit(power, color);
			m_attribute->SetLightFunction(Le);
		}

		Facet(){}
		Facet(
			const Vertex& v1, 
			const Vertex& v2, 
			const Vertex& v3, 
			const Vec3& color, 
			AttributeBSDF::BF BRDF, 
			AttributeBSDF::BF BTDF,
			const std::string& TextureName);

		Facet(
			const Vertex& v1, 
			const Vertex& v2, 
			const Vertex& v3,
			const Vec3& color,
			float Roughness, 
			float Metal, 
			float Trans,
			const std::string& TextureName);

		Facet(const Vertex& v1, 
			const Vertex& v2, 
			const Vertex& v3, 
			const Vec3& color,
			ATT_MATERIAL mat,
			const std::string& TextureName);

		Facet(
			const Vertex& v1,
			const Vertex& v2,
			const Vertex& v3,
			const Vec3& color,
			AttributeBSDF::BF BRDF,
			AttributeBSDF::BF BTDF);

		Facet(
			const Vertex& v1,
			const Vertex& v2,
			const Vertex& v3,
			const Vec3& color,
			float Roughness,
			float Metal,
			float Trans);

		Facet(const Vertex& v1,
			const Vertex& v2,
			const Vertex& v3,
			const Vec3& color,
			ATT_MATERIAL mat);

		Facet(const Vertex& v1, const Vertex& v2, const Vertex& v3);
		Facet(const Vertex& v1, const Vertex& v2, const Vertex& v3,const Vec3& color);

		virtual ~Facet();

	};


	

}


