#pragma once

#include"VoxelGrid.h"
#include"Grain.h"
#include"PhaseCal.h"
#include"Vec3.h"
#include"Otree.h"
#include<algorithm>
#include<string>
#include<sstream>
#include<vector>

namespace PTUtility{

	//つぶつぶ集合の管理をします
	class Aggrigate{

	public :
		using SC = OcTree::SC;

	private:

		std::vector<std::shared_ptr<iGrain>> m_Grains;
		std::vector<float> m_raito;
		float m_Traito;

		struct Pack{
			std::vector<Vec3> _Pos;
			float _R;
		};

		//ツリーが管理するつぶつぶと、Aggrigateが管理するつぶつぶはスケールが異なる　
		OcTree m_Otree;

	private:

		Pack m_Pack;
		VoxelGrid* m_VoxelGrid;
		float GetSphereR()const{
			return m_Pack._R;
		}

		bool LoadPack(const char* FileName);

		float GetSigma_b()const{
			static float r = -1;
			if (r < 0){
				r = (3.0f * GetPackingRaito()) / (4.0f * GetSphereR() * (1.00001 - GetPackingRaito()));
			}
			return r;
		}
		Vec3 GetLamda_s(int i)const{
			return GetSphereR() * m_Grains[i]->GetLamdaS();
		}
		float GetLamda_dt(int i)const{
			return GetSphereR() * m_Grains[i]->GetLamdaDt();
		}
		float GetBeta(int i)const{
			return m_Grains[i]->GetBeta();
		}
		float GetLamda_Beta(int i)const{
			float re = (GetLamda_dt(i) + (1.0f / GetSigma_b())) * ((1.00001 - GetBeta(i)) / (GetBeta(i))) + (1.0f / GetSigma_b());
			return re;
		}

		void InitCore(){

			//お好みのつぶつぶを追加できます
			//AddGrain(new Grain_TranslucentSphere(Vec3(0.99, 0.5, 0.3), m_Pack._R * m_VoxelGrid->GetVoxelWidth(), 0.1f, true), 0.7, "PHASE\\PTRANS.csv");
			//AddGrain(new Grain_GlassSphere(Vec3(0.99, 0.99, 0.99), m_Pack._R * m_VoxelGrid->GetVoxelWidth()), 1.1, "PHASE\\PGLASS.csv");
			//AddGrain(new Grain_MetalSphere(Vec3(0.66, 0.88, 0.5), m_Pack._R * m_VoxelGrid->GetVoxelWidth()), 0.3, "PHASE\\PMETAL.csv");
			AddGrain(new Grain_DiffuseSphere(Vec3(0.99, 0.5, 0.3), m_Pack._R * m_VoxelGrid->GetVoxelWidth()), 1.1, "PHASE\\PDIFFUSE.csv");
			//AddGrain(new Grain_DiffuseSphere(Vec3(0.99, 0.89, 0.89), m_Pack._R * m_VoxelGrid->GetVoxelWidth()), 1.0, "PHASE\\PDIFFUSE.csv");
			//AddGrain(new Grain_FromMesh(Grain_FromMesh::MATERIAL_TYPE::GLASS, Vec3(0.98, 0.99, 0.96), m_Pack._R * m_VoxelGrid->GetVoxelWidth(), "GRAIN\\GR1.obj"), 1.3, "PHASE\\PGRAIN2.csv");
			//AddGrain(new Grain_FromMesh(Grain_FromMesh::MATERIAL_TYPE::GLASS, Vec3(0.93, 0.99, 0.96), m_Pack._R * m_VoxelGrid->GetVoxelWidth(), "GRAIN\\GrainFB.obj"), 1.3, "PHASE\\PGLASS.csv");
			
			
			
			InitParam();
			InitSphereGraph();
		}

		void InitSphereGraph(){

			std::vector<Vec3> ddd = m_Pack._Pos;

			for (auto& c : ddd){
				c = (c - Vec3::Ones() * 0.5f) * m_VoxelGrid->GetVoxelWidth();
			}
			m_Otree.CreateOctree(ddd, m_Pack._R * m_VoxelGrid->GetVoxelWidth(), 4);
			ddd.clear();

		}


	protected:
		//謎の関数　マルチスレッドの闇を払うのかもしれないと感じた。
		void InitParam(){
			GetSigma_b();
			GetPackingRaito();
			GetPhase_g();
			GetAlpha_s();
			GetSigma_t();
		}

	public:

		




		//Widthには、ボクセルの幅を渡しておけばよいとのこと
		explicit Aggrigate(float AggrigateWidth);

		virtual ~Aggrigate();

		//ポイントクラウドデータから初期化　通常はこっちかな？
		bool InitAggrigate(const char* FileName);

		//ボクセライズからの初期化
		bool InitAggrigate(const Scene& scene, int NVoxel);


		int Intersect(const Ray& r, std::vector<const VoxelGrid::VPP>& result);
		Vec3 GetCenter(const VoxelGrid::VoxelData* v);


		int GetSphereFromRay(
			const Vec3& VoxelCenter, 
			const Ray& ray, 
			std::vector<OcTree::SC>& result);

		void AddGrain(
			iGrain* pgrain, 
			float PackRaito, 
			const std::string& PhaseFileName){
			m_raito.push_back(PackRaito);
			m_Traito += PackRaito;
			pgrain->InitGrain(m_VoxelGrid->GetVoxelWidth(), PhaseFileName.c_str());

			m_Grains.push_back(std::shared_ptr<iGrain>(pgrain));
		}

		bool NextRay(
			int GrainID,
			const Ray& r,
			const Vec3& Center,
			int RandomSeed,
			Ray& result,
			float& PDF,
			RandomMT& MT);

		int RandomSelectGrain(int Seed){
			RandomMT MT(Seed);

			int N = m_Grains.size();
			float tt = FLT_MIN;
			float ttt = MT.genrand64_real2() * m_Traito;

			for (int i = 0; i < N; i++){
				tt += m_raito[i];
				if (tt > ttt){
					return i;
				}
			}
			return 0;
		}

		const iGrain* GetGrain(int i)const{
			return m_Grains[i].get();
		}

		float GetVoxelWidth()const{
			return m_VoxelGrid->GetVoxelWidth();
		}

		float GetPackingRaito()const{
			static float r = -1.0f;
			if (r < 0){
				return (r = (m_Pack._Pos.size() * (4.0f / 3.0f) * M_PI * pow(GetSphereR(), 3)));
			}
			else{
				return r;
			}
		}

		


		float PhaseF(const Vec3& in, const Vec3& out){
			float result = 0.0f;

			for (int i = 0; i < m_Grains.size(); i++){
				result += (m_raito[i]) * m_Grains[i]->Phase(in, out);
			}
			return result / m_Traito;
		}


		float GetPhase_g(){
			static float r = -10;

			if (r < -9){
				float result = 0;
				//PhaseFunctionの積分を行う

				for (int th = 0; th < 314; th++){
					for (int ph = 0; ph < 628; ph++){
						float s_th = sinf((float)th*0.01);
						float c_th = cosf((float)th*0.01);
						float s_ph = sinf((float)ph*0.01);
						float c_ph = cosf((float)ph*0.01);

						//入射方向(1,0,0)
						float dot = s_th * c_ph;

						result += s_th *
							dot *
							PhaseF(Vec3(1, 0, 0), Vec3(s_th*c_ph, s_th*s_ph, c_th));

					}
				}
				r = result * 0.0000999;

			}
			return r;
		}
		
		Vec3 GetAlpha_s(){
			static Vec3 r(-1, -1, -1);
			if (r.x() < 0){
				Vec3 result(0, 0, 0);
				for (int i = 0; i < m_Grains.size(); i++){
					result += (m_raito[i]) * m_Grains[i]->GetAlphaS();
				}
				r = result / m_Traito;
			}
			
			return r;
		}
		Vec3 GetSigma_t(){

			static Vec3 r(-2, -3, -4);
			if (r.x() < 0){

				Vec3 result(0, 0, 0);
				for (int i = 0; i < m_Grains.size(); i++){
					for (int c = 0; c < 3; c++){
						result[c] += m_raito[i] * (GetLamda_Beta(i) + m_Grains[i]->GetAlphaS()[c] * GetLamda_s(i)[c]);
					}
				}
				r[0] = 1.0f / (m_VoxelGrid->GetVoxelWidth() * result[0] / m_Traito);
				r[1] = 1.0f / (m_VoxelGrid->GetVoxelWidth() * result[1] / m_Traito);
				r[2] = 1.0f / (m_VoxelGrid->GetVoxelWidth() * result[2] / m_Traito);

			}
			return r;
		}




	};


}