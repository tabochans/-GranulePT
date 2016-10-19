#pragma once
#include"Ray.h"
#include"Vec3.h"
#include"Intersect.h"
#include"ImportanceSample.h"
#include"PhaseCal.h"
#include"OBVH.h"


namespace PTUtility{

	class Scene;

	//ó±àÍÇ¬ÇÃè⁄ç◊Ç»èÓïÒ
	class iGrain{

		friend PhaseEstimater;

	public:
		//PhaseFunctionÇÀ
		typedef float(*PhaseF)(const Vec3& in, const Vec3& out);

		//PhaseFunctionÇªÇÃÇQ
		

	private:

	protected:

		PhaseEstimater* m_PPHE;

		Vec3 m_Color;

		Vec3 m_Lamda_s;
		float m_Lamda_dt;
		Vec3 m_Alpha_s;
		Vec3 m_Sigma_t;
		float m_Beta;
		PhaseF m_PF;

		//GrainÇéÊÇËàÕÇﬁBounding Sphere ÇÃîºåa
		float m_R;

		

	public:

		virtual void InitGrain(float VoxelWidth, const char* PhaseFileName);

		float GetR()const{ return m_R; }
		float GetBeta()const{
			return m_Beta;
		}

		const Vec3& GetColor()const{
			return m_Color;
		}

		const Vec3& GetLamdaS()const{
			return m_Lamda_s;
		}

		float GetLamdaDt()const{
			return m_Lamda_dt;
		}

		const Vec3& GetAlphaS()const{
			return m_Alpha_s;
		}

		const Vec3& GetSigmaT()const{
			return m_Sigma_t;
		}

		float PPhase(float cos)const{
			return m_PPHE->GetFromCos(cos);
		}

		virtual float Phase(const Vec3& in, const Vec3& out)const{
			if (m_PPHE == NULL){
				return m_PF(in, out);
			}
			else{
				return PPhase(in.dot(out));
			}
		}


		iGrain(float Radi, const char* FileName) : m_Lamda_s(1, 1, 1), m_Lamda_dt(0), m_Alpha_s(1, 1, 1), m_Sigma_t(0, 0, 0), m_R(Radi){
			m_PPHE = NULL;
		}
		explicit iGrain(float Radi) : m_Lamda_s(1, 1, 1), m_Lamda_dt(0), m_Alpha_s(1, 1, 1), m_Sigma_t(0, 0, 0), m_R(Radi){
			m_PPHE = NULL;
		}
		virtual ~iGrain(){
			delete m_PPHE;
		}

		virtual bool IfIntersect(const Ray& r, const Vec3& RCenter, int RandomSeed)const = 0;
		
		virtual bool NextRay(
			const Ray& r,
			const Vec3& Center,
			int RandomSeed,
			Ray& result,
			float& PDF,
			RandomMT& MT)const = 0;

		virtual Vec3 GetGrainCenter() const = 0;
		virtual float GetGrainRadius() const = 0;

	};

	class Grain_Sphere : public iGrain{

	private:

	protected:

	public:

		virtual Vec3 GetGrainCenter()const{
			return Vec3::Zero();
		}
		virtual float GetGrainRadius()const{
			return m_R * 0.9;
		}

		virtual bool IfIntersect(const Ray& r, const Vec3& RCenter, int RandomSeed)const;

		virtual bool NextRay(
			const Ray& r,
			const Vec3& Center,
			int RandomSeed,
			Ray& result,
			float& PDF,
			RandomMT& MT)const = 0;


		explicit Grain_Sphere(float Radi) : iGrain(Radi){

		}
		virtual ~Grain_Sphere(){

		}

	};


	class Grain_TranslucentSphere : public Grain_Sphere{

	private:

		//ÇøÇÂÇ¡Ç∆ÇÀ
		static float Phase(const Vec3& in, const Vec3& out){
			return 0.25 / M_PI;
		}

		float m_SSS;
		bool m_Diffuse;

	protected:

	public:

		using Grain_Sphere::IfIntersect;


		virtual bool NextRay(
			const Ray& r,
			const Vec3& Center,
			int RandomSeed,
			Ray& result,
			float& PDF,
			RandomMT& MT)const;


		explicit Grain_TranslucentSphere(
			const Vec3& Color, float Radi, float SSS, bool Diffuse) : Grain_Sphere(Radi){
			m_Alpha_s = m_Color = Color;
			m_Lamda_s = Vec3(1.44, 1.44, 1.44);
			m_Lamda_dt = 0.0f;
			m_Sigma_t = Vec3(1, 1, 1);
			m_Beta = 1.0f;
			m_PF = Grain_TranslucentSphere::Phase;
			m_Diffuse = Diffuse;
			m_SSS = SSS * Radi;
		}
		virtual ~Grain_TranslucentSphere(){

		}

	};


	class Grain_GlassSphere : public Grain_Sphere{

	private:

		//ÇøÇÂÇ¡Ç∆ÇÀ
		static float Phase(const Vec3& in, const Vec3& out){
			float c = in.dot(out);
			c = std::max(-0.99999999f, std::min(0.999999f, c));
			return 25.0f * ( pow(c, 4)*(2 + c) + 0.01f) / (41.0f * M_PI);
		}

	protected:

	public:

		using Grain_Sphere::IfIntersect;
		

		virtual bool NextRay(
			const Ray& r,
			const Vec3& Center,
			int RandomSeed,
			Ray& result,
			float& PDF,
			RandomMT& MT)const;


		explicit Grain_GlassSphere(const Vec3& Color, float Radi) : Grain_Sphere(Radi){
			m_Alpha_s = m_Color = Color;
			m_Lamda_s = Vec3(1.44, 1.44, 1.44);
			m_Lamda_dt = 0.0f;
			m_Sigma_t = Vec3(1, 1, 1);
			m_Beta = 1.0f;
			m_PF = Grain_GlassSphere::Phase;
		}
		virtual ~Grain_GlassSphere(){

		}

	};



	class Grain_MetalSphere : public Grain_Sphere{

	private:

		static float Phase(const Vec3& in, const Vec3& out){
			return (0.999999) / (4.0f*M_PI);
		}


	protected:

	public:

		using Grain_Sphere::IfIntersect;

		virtual bool NextRay(
			const Ray& r,
			const Vec3& Center,
			int RandomSeed,
			Ray& result,
			float& PDF,
			RandomMT& MT)const;


		explicit Grain_MetalSphere(const Vec3& Color, float Radi) : Grain_Sphere(Radi){
			m_Alpha_s = m_Color = Color;
			m_Lamda_s = Vec3(0.01, 0.01, 0.01);
			m_Lamda_dt = 0.0f;
			m_Sigma_t = Vec3(1, 1, 1);
			m_Beta = 0.99f;
			m_PF = Grain_MetalSphere::Phase;
		}
		virtual ~Grain_MetalSphere(){

		}

	};



	class Grain_DiffuseSphere : public Grain_Sphere{

	private:

		static float Phase(const Vec3& in, const Vec3& out){
			return pow(-in.dot(out) + 1, 2) / (2.658f * 2.0f * M_PI);
		}


	protected:


	public:


		using Grain_Sphere::IfIntersect;

		virtual bool NextRay(
			const Ray& r,
			const Vec3& Center,
			int RandomSeed,
			Ray& result,
			float& PDF,
			RandomMT& MT)const;


		explicit Grain_DiffuseSphere(const Vec3& Color, float Radi) : Grain_Sphere(Radi){
			m_Alpha_s = m_Color = Color;
			m_Lamda_s = Vec3(0.01, 0.01, 0.01);
			m_Lamda_dt = 0.0f;
			m_Sigma_t = Vec3(1, 1, 1);
			m_Beta = 0.99f;
			m_PF = Grain_DiffuseSphere::Phase;


		}
		virtual ~Grain_DiffuseSphere(){

		}

	};

	//å¥ì_Ç™íÜêSÇ∂Ç·Ç»Ç¢Ç∆Ç‚ÇŒÇ¢
	class Grain_FromMesh : public iGrain{

	public :
		enum MATERIAL_TYPE{
			Diffuse,
			METAL,
			GLASS
		};


	private:

		static float Phase(const Vec3& in, const Vec3& out){
			return 0;
		}

		Scene* m_scene;
		Vec3 m_CenterPos;
		MATERIAL_TYPE m_mat;

		bool PNextRay_DIFFUSE(
			const Ray& r,
			const Vec3& Center,
			int RandomSeed,
			Ray& result,
			float& PDF,
			RandomMT& MT)const;
		bool PNextRay_METAL(
			const Ray& r,
			const Vec3& Center,
			int RandomSeed,
			Ray& result,
			float& PDF,
			RandomMT& MT)const;
		bool PNextRay_GLASS(
			const Ray& r,
			const Vec3& Center,
			int RandomSeed,
			Ray& result,
			float& PDF,
			RandomMT& MT)const;


	protected:

	public:

		virtual void InitGrain(float VoxelWidth, const char* PhaseFileName);

		virtual Vec3 GetGrainCenter()const;
		virtual float GetGrainRadius()const;

		virtual bool IfIntersect(const Ray& r, const Vec3& RCenter, int RandomSeed)const;

		virtual bool NextRay(
			const Ray& r,
			const Vec3& Center,
			int RandomSeed,
			Ray& result,
			float& PDF,
			RandomMT& MT)const;

		Grain_FromMesh(MATERIAL_TYPE type, const Vec3& Color, float Radi, const char* FileName);
		virtual ~Grain_FromMesh();

	};






}



















