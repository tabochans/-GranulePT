#include"GPathTracing.h"


using namespace PTUtility;




class GPathTrace::GPTCore{

public:

	enum RENDERMODE{

		EPT_INTER,
		EPT_INTRA,
		VPT,
		DA
	};


	struct PTDATA{
		struct VERTEX{
			Vec3 Normal;
			Vec3 Pos;
			Vec2 UV;
			std::string TextureName;
			Vec3 Color;
			Vec3 Emit;
			Vec3 in_Color;
			Vec3 d_Color;
			Vec3 ibl_Color;
			Vec3 weight;
			Vec3 NN;
			int EPT_Depth;

			VERTEX() :
				Normal(0, 1, 0),
				Pos(0, 0, 0),
				UV(0, 0),
				TextureName(""),
				Color(0, 0, 0),
				Emit(0, 0, 0),
				in_Color(0, 0, 0),
				d_Color(0, 0, 0),
				ibl_Color(0, 0, 0),
				weight(1, 1, 1),
				NN(0, 1, 0),
				EPT_Depth(0){

			}

		};
		struct PDF{
			float Pdf;
			PDF() : Pdf(1.0f){

			}
		};

		struct FLAG{
			bool NextEvent;
			bool IFInObject;
			FLAG() : NextEvent(false), IFInObject(false){

			}
		};

		struct RAY{
			Ray NextRay;
			Ray CRay;
			float Index_From;
			float Index_To;
			RAY() : NextRay(Vec3::Zero(), Vec3(0, 1, 0)), CRay(Vec3::Zero(), Vec3(0, 1, 0)), Index_From(1.0f), Index_To(1.0f){
			}
		};

		struct EPT_INTRA{
			int PID;
			int GrainIDX;
			Vec3 In_Pos;
			Vec3 In_Dir;
			Vec3 GCenter;
			EPT_INTRA() : PID(0), GrainIDX(0), In_Pos(0, 0, 1), In_Dir(0, 0, -1), GCenter(0, 0, 0){

			}

		};

		struct DA{
			Vec3 Pos;
			Vec3 EN;
			int DA_Count;
			Vec3 DA_Weight;
			bool IfNextEvent;
			DA() : Pos(0, 0, 0), EN(0, 1, 0), DA_Count(0), DA_Weight(1, 1, 1), IfNextEvent(false){

			}
		};

		int Depth;
		VERTEX Vertex;
		PDF Pdf;
		FLAG Flag;
		RAY Ray;
		EPT_INTRA EPTintra;
		DA Diffuse;
		Vec3 Result;
		bool IfLoop;

		PTDATA() : Vertex(), Pdf(), Depth(0), Flag(), Ray(), EPTintra(), Result(0, 0, 0), IfLoop(true){

		}

	};



	float NextEventIBL(
		int NShadowRay,
		const Scene& scene,
		RandomMT& MT,
		const Vec3& Pos,
		const Vec3& Normal,
		const Ray& ray,
		Vec3& result){

		if (NShadowRay <= 0){
			return 0.0f;
		}

		int SpL = NShadowRay;
		result = Vec3::Zero();

		for (int l = 0; l < SpL; l++){

			Vec3 ndir;
			ISample::NextDirDiffuse(ray, Normal, ndir, MT);

			//光源から明示的にサンプル///////
			std::vector<AbsG::RF> ffa;
			scene.GetFacetFromRay(Ray(Pos + 0.0001 * Normal, ndir), ffa);

			if (ffa.size() == 0){
				result += m_IBL->GetColor1(ndir);
			}


		}
		result /= ((float)SpL);

		return 0.0f;
	}


	RENDERMODE m_Mode;

	iIBL const *m_IBL;

	Aggrigate* m_agg;

	typedef float(*BRDF)(const Vec3& in, const Vec3& out);

	const int kDepth;

	char* img;

	Vec3 m_st_d;


	GPTCore() : kDepth(50000){
	}
	~GPTCore(){
		delete[] img;
	}

	Vec3 GetRdt(
		PTDATA& Data,
		float R,
		const Vec3& Pos_i,
		const Vec3& Normal_i,
		const Vec3& Pos_a,
		const Vec3& Normal_a
		){


		////////////////////////////////////////////////
		static bool init = false;
		const static float C11 = 0.5 * (-9.23372 + 22.2272 - 20.9292 + 10.2291 - 2.54396 + 0.254913);
		const static float C22 = 0.333333 * (-1641.1 + 135.926 - 656.175 + 1376.53 + 1213.67 -
			568.556 + 164.798 - 27.0181 + 1.91826);
		const static float g = m_agg->GetPhase_g();
		const static Vec3 Sigma_s = Vec3(
			m_agg->GetSigma_t()[0] * m_agg->GetAlpha_s()[0],
			m_agg->GetSigma_t()[1] * m_agg->GetAlpha_s()[1],
			m_agg->GetSigma_t()[2] * m_agg->GetAlpha_s()[2]
			);
		const static Vec3 Sigma_s_d = (1.0f - g) *
			Vec3(
			m_agg->GetAlpha_s()[0] * m_agg->GetSigma_t()[0],
			m_agg->GetAlpha_s()[1] * m_agg->GetSigma_t()[1],
			m_agg->GetAlpha_s()[2] * m_agg->GetSigma_t()[2]
			);

		const static Vec3 Sigma_t_d = Sigma_s_d + (m_agg->GetSigma_t() - Sigma_s);
		const static Vec3 Alpha_s(Sigma_s_d[0] / Sigma_t_d[0], Sigma_s_d[1] / Sigma_t_d[1], Sigma_s_d[2] / Sigma_t_d[2]);
		const static Vec3 Sigma_a = (m_agg->GetSigma_t() - Sigma_s);
		m_st_d = Sigma_t_d;
		const static Vec3 L(1.0f / Sigma_t_d[0], 1.0f / Sigma_t_d[1], 1.0f / Sigma_t_d[2]);
		const static float A = (1.0f + 3.0f * C22) / (1.0f - 2.0f * C11);
		const static Vec3 D(
			(2.0f * Sigma_a[0] + Sigma_s_d[0]) / (3.0f * powf(Sigma_a[0] + Sigma_s_d[0], 2.0f)),
			(2.0f * Sigma_a[1] + Sigma_s_d[1]) / (3.0f * powf(Sigma_a[1] + Sigma_s_d[1], 2.0f)),
			(2.0f * Sigma_a[2] + Sigma_s_d[2]) / (3.0f * powf(Sigma_a[2] + Sigma_s_d[2], 2.0f))
			);
		const static Vec3 Sigma_tr(
			sqrt(Sigma_a[0] / D[0]),
			sqrt(Sigma_a[1] / D[1]),
			sqrt(Sigma_a[2] / D[2]));
		const static Vec3 Zb = 2.0f * A * D;
		const static float Cf = 0.25 * (1.0f - 2.0f * C11);
		const static float Ce = 0.5 * (1.0f - 3.0f * C22);
		////////////////////////////////////////////////


		float zz = abs(Normal_i.dot(Pos_i - Data.Diffuse.Pos));
		float Dist_i = (Data.Diffuse.Pos - Pos_i).norm();
		float Dist_a = (Data.Diffuse.Pos - Pos_a).norm();
		float Za = abs(Normal_a.dot(Pos_a - Data.Diffuse.Pos));
		float d = zz + Za;
		Vec3 Zr[3];
		Zr[0] = L;
		Zr[1] = 2.0f * 1.0f * (Vec3::Ones() * d + 2.0f * Zb) + L;
		Zr[2] = (-2.0f * 1.0f * (Vec3::Ones() * d + 2.0f * Zb) + L);
		float r = R;

		Vec3 Zv[3];
		Zv[0] = -L - 2.0f * Zb;
		Zv[1] = 2.0f * 1.0f * (Vec3::Ones() * d + 2.0f * Zb) - L - 2.0f * Zb;
		Zv[2] = (-2.0f * 1.0f * (Vec3::Ones() * d + 2.0f * Zb) - L - 2.0f * Zb);

		Vec3 dr[3];
		dr[0] = Vec3(sqrt(r*r + Zr[0][0] * Zr[0][0]), sqrt(r*r + Zr[0][1] * Zr[0][1]), sqrt(r*r + Zr[0][2] * Zr[0][2]));
		dr[1] = Vec3(sqrt(r*r + Zr[1][0] * Zr[1][0]), sqrt(r*r + Zr[1][1] * Zr[1][1]), sqrt(r*r + Zr[1][2] * Zr[1][2]));
		dr[2] = Vec3(sqrt(r*r + Zr[2][0] * Zr[2][0]), sqrt(r*r + Zr[2][1] * Zr[2][1]), sqrt(r*r + Zr[2][2] * Zr[2][2]));

		Vec3 dv[3];
		dv[0] = Vec3(sqrt(Zv[0][0] * Zv[0][0] + r*r), sqrt(Zv[0][1] * Zv[0][1] + r*r), sqrt(Zv[0][2] * Zv[0][2] + r*r));
		dv[1] = Vec3(sqrt(Zv[1][0] * Zv[1][0] + r*r), sqrt(Zv[1][1] * Zv[1][1] + r*r), sqrt(Zv[1][2] * Zv[1][2] + r*r));
		dv[2] = Vec3(sqrt(Zv[2][0] * Zv[2][0] + r*r), sqrt(Zv[2][1] * Zv[2][1] + r*r), sqrt(Zv[2][2] * Zv[2][2] + r*r));

		Vec3 temp(0, 0, 0);


		for (int i = 0; i < 3; i++){
			for (int c = 0; c < 3; c++){

				temp[c] +=

					Alpha_s[c] * Alpha_s[c] / (4.0f * M_PI) *
					(
					(
					Ce * ((Zr[i][c] * (Sigma_tr[c] * dr[i][c] + 1.0f)) / (dr[i][c] * dr[i][c])) + Cf
					) * exp(-Sigma_tr[c] * dr[i][c]) / (dr[i][c] * dr[i][c])
					-
					(
					Ce * ((Zv[i][c] * (Sigma_tr[c] * dv[i][c] + 1.0f)) / (dv[i][c] * dv[i][c])) + Cf
					) * exp(-Sigma_tr[c] * dv[i][c]) / (dv[i][c] * dv[i][c])
					);

			}
		}


		return temp;
	}


	Vec3 GetRd_DA(
		PTDATA& Data,
		const Vec3& Pos_i,
		const Vec3& Normal_i,
		const Vec3& Pos_a,
		const Vec3& Normal_a
		){

		////////////////////////////////////////////////
		static bool init = false;
		const static float C11 = 0.5 * (-9.23372 + 22.2272 - 20.9292 + 10.2291 - 2.54396 + 0.254913);
		const static float C22 = 0.333333 * (-1641.1 + 135.926 - 656.175 + 1376.53 + 1213.67 -
			568.556 + 164.798 - 27.0181 + 1.91826);
		const static float g = m_agg->GetPhase_g();
		const static Vec3 Sigma_s = Vec3(
			m_agg->GetSigma_t()[0] * m_agg->GetAlpha_s()[0],
			m_agg->GetSigma_t()[1] * m_agg->GetAlpha_s()[1],
			m_agg->GetSigma_t()[2] * m_agg->GetAlpha_s()[2]
			);
		const static Vec3 Sigma_s_d = (1.0f - g) *
			Vec3(
			m_agg->GetAlpha_s()[0] * m_agg->GetSigma_t()[0],
			m_agg->GetAlpha_s()[1] * m_agg->GetSigma_t()[1],
			m_agg->GetAlpha_s()[2] * m_agg->GetSigma_t()[2]
			);

		const static Vec3 Sigma_t_d = Sigma_s_d + (m_agg->GetSigma_t() - Sigma_s);
		const static Vec3 Alpha_s(Sigma_s_d[0] / Sigma_t_d[0], Sigma_s_d[1] / Sigma_t_d[1], Sigma_s_d[2] / Sigma_t_d[2]);
		const static Vec3 Sigma_a = (m_agg->GetSigma_t() - Sigma_s);
		m_st_d = Sigma_t_d;
		const static Vec3 L(1.0f / Sigma_t_d[0], 1.0f / Sigma_t_d[1], 1.0f / Sigma_t_d[2]);
		const static float A = (1.0f + 3.0f * C22) / (1.0f - 2.0f * C11);
		const static Vec3 D(
			(2.0f * Sigma_a[0] + Sigma_s_d[0]) / (3.0f * powf(Sigma_a[0] + Sigma_s_d[0], 2.0f)),
			(2.0f * Sigma_a[1] + Sigma_s_d[1]) / (3.0f * powf(Sigma_a[1] + Sigma_s_d[1], 2.0f)),
			(2.0f * Sigma_a[2] + Sigma_s_d[2]) / (3.0f * powf(Sigma_a[2] + Sigma_s_d[2], 2.0f))
			);
		const static Vec3 Sigma_tr(
			sqrt(Sigma_a[0] / D[0]),
			sqrt(Sigma_a[1] / D[1]),
			sqrt(Sigma_a[2] / D[2]));
		const static Vec3 Zb = 2.0f * A * D;
		const static float Cf = 0.25 * (1.0f - 2.0f * C11);
		const static float Ce = 0.5 * (1.0f - 3.0f * C22);

		if (init == false){
			init = true;
			return Vec3::Zero();
		}
		////////////////////////////////////////////////


		float zz = abs(Normal_i.dot(Pos_i - Data.Diffuse.Pos));
		float Dist_i = (Data.Diffuse.Pos - Pos_i).norm();
		float Dist_a = (Data.Diffuse.Pos - Pos_a).norm();
		float Za = abs(Normal_a.dot(Pos_a - Data.Diffuse.Pos));
		float d = zz + Za;
		Vec3 Zr[3];
		Zr[0] = L;
		Zr[1] = 2.0f * 1.0f * (Vec3::Ones() * d + 2.0f * Zb) + L;
		Zr[2] = (-2.0f * 1.0f * (Vec3::Ones() * d + 2.0f * Zb) + L);
		float r = sqrt(Dist_i*Dist_i + FLT_MIN - zz * zz);

		Vec3 Zv[3];
		Zv[0] = -L - 2.0f * Zb;
		Zv[1] = 2.0f * 1.0f * (Vec3::Ones() * d + 2.0f * Zb) - L - 2.0f * Zb;
		Zv[2] = (-2.0f * 1.0f * (Vec3::Ones() * d + 2.0f * Zb) - L - 2.0f * Zb);

		Vec3 dr[3];
		dr[0] = Vec3(sqrt(r*r + Zr[0][0] * Zr[0][0]), sqrt(r*r + Zr[0][1] * Zr[0][1]), sqrt(r*r + Zr[0][2] * Zr[0][2]));
		dr[1] = Vec3(sqrt(r*r + Zr[1][0] * Zr[1][0]), sqrt(r*r + Zr[1][1] * Zr[1][1]), sqrt(r*r + Zr[1][2] * Zr[1][2]));
		dr[2] = Vec3(sqrt(r*r + Zr[2][0] * Zr[2][0]), sqrt(r*r + Zr[2][1] * Zr[2][1]), sqrt(r*r + Zr[2][2] * Zr[2][2]));

		Vec3 dv[3];
		dv[0] = Vec3(sqrt(Zv[0][0] * Zv[0][0] + r*r), sqrt(Zv[0][1] * Zv[0][1] + r*r), sqrt(Zv[0][2] * Zv[0][2] + r*r));
		dv[1] = Vec3(sqrt(Zv[1][0] * Zv[1][0] + r*r), sqrt(Zv[1][1] * Zv[1][1] + r*r), sqrt(Zv[1][2] * Zv[1][2] + r*r));
		dv[2] = Vec3(sqrt(Zv[2][0] * Zv[2][0] + r*r), sqrt(Zv[2][1] * Zv[2][1] + r*r), sqrt(Zv[2][2] * Zv[2][2] + r*r));

		Vec3 temp(0, 0, 0);

		//for (int i = 0; i < 3; i++){
		//	for (int c = 0; c < 3; c++){
		//		float C1 = (Zr[i][c]) * (1.0f + Sigma_tr[c] * dr[i][c]) * exp(-Sigma_tr[c] * dr[i][c]) / pow(dr[i][c], 3);
		//		float C2 = (Zv[i][c]) * (1.0f + Sigma_tr[c] * dv[i][c]) * exp(-Sigma_tr[c] * dv[i][c]) / pow(dv[i][c], 3);
		//		temp[c] += C1;
		//		temp[c] -= C2;
		//	}
		//}
		//MultVec(Alpha_s / (4.0f * M_PI * 100.0f), temp);


		for (int i = 0; i < 3; i++){
			for (int c = 0; c < 3; c++){

				temp[c] +=

					Alpha_s[c] * Alpha_s[c] / (4.0f * M_PI) *
					(
					(
					Ce * ((Zr[i][c] * (Sigma_tr[c] * dr[i][c] + 1.0f)) / (dr[i][c] * dr[i][c])) + Cf
					) * exp(-Sigma_tr[c] * dr[i][c]) / (dr[i][c])
					-
					(
					Ce * ((Zv[i][c] * (Sigma_tr[c] * dv[i][c] + 1.0f)) / (dv[i][c] * dv[i][c])) + Cf
					) * exp(-Sigma_tr[c] * dv[i][c]) / (dv[i][c])
					);

			}
		}





		static bool ff = false;
		if (!ff){
			std::ofstream file;
			file.open("dd.csv");
			float L = 1.0f / m_st_d[0];
			for (int f = 0; f < 2000; f++){

				file << f / 3000.0f;
				file << ",";
				Vec3 rd = GetRdt(Data, f / 3000.0f, Pos_i, Normal_i, Pos_a, Normal_a);
				file << rd[0] << ", ";
				file << rd[1] << ", ";
				file << rd[2];
				file << std::endl;
			}
			file.close();
			ff = true;
		}



		return temp / (1.0f);
	}



	//ちょっといみふ
	Vec3 GetRd_DiffusionDipole(
		PTDATA& Data,
		const Vec3& exitP
		){
		Vec3 result(0, 0, 0);

		float g = m_agg->GetPhase_g();
		const Vec3 sigma_s_d = m_agg->GetSigma_t()*(1.0f - g);
		const Vec3 sigma_a = Vec3(0.1, 0.1, 0.1);
		const Vec3 sigma_t_d = sigma_s_d + sigma_a;
		const Vec3 alpha_d = Vec3(
			sigma_s_d[0] / sigma_t_d[0],
			sigma_s_d[1] / sigma_t_d[1],
			sigma_s_d[2] / sigma_t_d[2]);

		const Vec3 Zr = Vec3(1.0f / sigma_t_d[0], 1.0f / sigma_t_d[1], 1.0f / sigma_t_d[2]);
		const float A = 1.0f;
		const Vec3 Zv = -Zr - Zr * 4.0f * A / 3.0f;

		const float s = (Data.Diffuse.Pos - exitP).norm();
		const float sn = abs((Data.Diffuse.Pos - exitP).dot(Data.Diffuse.EN));
		const float r = sqrt(s*s - sn*sn) + 0.01;


		const Vec3 dv = 0.001 * Vec3::Ones() + Vec3(
			sqrt(r*r + Zv[0] * Zv[0]), sqrt(r*r + Zv[1] * Zv[1]), sqrt(r*r + Zv[2] * Zv[2]));

		const Vec3 dr = 0.001 * Vec3::Ones() + Vec3(
			sqrt(r*r + Zr[0] * Zr[0]), sqrt(r*r + Zr[1] * Zr[1]), sqrt(r*r + Zr[2] * Zr[2]));

		const Vec3 sigma_tr = 1.732 * Vec3(
			sqrt(sigma_a[0] * sigma_t_d[0]),
			sqrt(sigma_a[1] * sigma_t_d[1]),
			sqrt(sigma_a[2] * sigma_t_d[2]));

		const Vec3 C1(
			Zr[0] * (sigma_tr[0] + 1.0f / dr[0]),
			Zr[1] * (sigma_tr[1] + 1.0f / dr[1]),
			Zr[2] * (sigma_tr[2] + 1.0f / dr[2]));

		const Vec3 C2(
			Zv[0] * (sigma_tr[0] + 1.0f / dv[0]),
			Zv[1] * (sigma_tr[1] + 1.0f / dv[1]),
			Zv[2] * (sigma_tr[2] + 1.0f / dv[2]));

		for (int i = 0; i < 3; i++){
			result[i] = alpha_d[i] *
				(
				C1[i] * exp(-sigma_tr[i] * dr[i]) / (dr[i] * dr[i]) +
				C2[i] * exp(-sigma_tr[i] * dv[i]) / (dv[i] * dv[i])
				) /
				(4.0f*M_PI);
		}

		return result;
	}

	void DiffusionAproximate_Dipole(
		PTDATA& Data,
		const Scene& scene,
		const Camera& camera,
		RandomMT& MT
		){

		Data.Vertex.Emit = Vec3::Zero();

		Vec3 NNormal(0, 1, 0);
		Vec3 Wt(1, 1, 1);

		int Ns = 1;
		for (int ss = 0; ss < Ns; ss++){

			std::vector<AbsG::RF> fa;
			Vec3 RD;
			ISample::NextDirSphere(Data.Ray.CRay, Data.Vertex.Normal, RD, MT);
			int nf = scene.GetFacetFromRay(Ray(Data.Diffuse.Pos, RD), fa);
			if (nf > 0){
				Vec3 Ni =
					fa[0]._u * (fa[0]._pF->GetNormal(1) - fa[0]._pF->GetNormal(0)) +
					fa[0]._v * (fa[0]._pF->GetNormal(2) - fa[0]._pF->GetNormal(0)) +
					fa[0]._pF->GetNormal(0);
				Ni.normalize();

				if (Ni.dot(RD) > FLT_MIN){
					
					Vec3 Pi = (0.000001 + fa[0]._t) * RD + Data.Diffuse.Pos;

					Data.Ray.NextRay.m_Org = Pi;
					NNormal = Ni;
					Data.Diffuse.EN = Ni;
					int nf = scene.GetFacetFromRay(Ray(Data.Diffuse.Pos, -RD), fa);
					Vec3 Na =
						fa[0]._u * (fa[0]._pF->GetNormal(1) - fa[0]._pF->GetNormal(0)) +
						fa[0]._v * (fa[0]._pF->GetNormal(2) - fa[0]._pF->GetNormal(0)) +
						fa[0]._pF->GetNormal(0);
					Na.normalize();
					Vec3 Pa = (0.000001 + fa[0]._t) * RD + Data.Diffuse.Pos;

					Wt = GetRd_DA(Data, Pi, Ni, Pa, Na);


					//Vec3 ir = 10.0f * std::max(-Ni.y(), 0.0f) * Vec3::Ones();
					//MultVec(Data.Vertex.weight, ir);
					//MultVec(Wt, ir);
					//Data.Result += ir;

				}
				else{
				}
			}
			else{
			}
		}

		//よーわからん4.0fが正しい気がするんだがー
		//Wt *= (2.0f / (float)Ns);
		Wt *= (4.0f / (float)Ns);

		MultVec(Wt, Data.Vertex.weight);

		Data.Diffuse.DA_Count++;
		//Data.Diffuse.DA_Count = 10;


		//3回もカウントしたらもういいだろという風潮
		//biased 
		if (Data.Diffuse.DA_Count > 4){
			Data.IfLoop = false;
			return;
		}
		else{
			m_Mode = VPT;
			ISample::NextDirDiffuse(Data.Ray.CRay, NNormal, Data.Ray.NextRay.m_Dir, MT);
			Data.Ray.NextRay.m_Org += 0.001 * Data.Ray.NextRay.m_Dir;
		}

		return;

	}

	void VolumePathTrace(
		PTDATA& Data,
		const Scene& scene,
		const Camera& camera,
		RandomMT& MT
		){


		std::vector<AbsG::RF> fa;
		int nf = scene.GetFacetFromRay(Data.Ray.CRay, fa);

		//一番近いポリゴンとの交差を見る
		int fi = 0;

		if (nf == 0){
			//交差が取れない
			Data.Vertex.Emit = m_IBL->GetColor1(Data.Ray.CRay.m_Dir);

			MultVec(Data.Vertex.weight, Data.Vertex.Emit);

			if (!Data.Diffuse.IfNextEvent){
				Data.Result += Data.Vertex.Emit;
			}

			Data.IfLoop = false;
			return;
		}
		Data.Diffuse.IfNextEvent = false;

		Data.Vertex.Normal = fa[fi]._u * (fa[fi]._pF->GetNormal(1) - fa[fi]._pF->GetNormal(0)) +
			fa[fi]._v * (fa[fi]._pF->GetNormal(2) - fa[fi]._pF->GetNormal(0)) +
			fa[fi]._pF->GetNormal(0);
		Data.Vertex.Normal.normalize();

		Data.Vertex.Pos = fa[fi]._t * Data.Ray.CRay.m_Dir + Data.Ray.CRay.m_Org;

		Data.Vertex.UV =
			fa[fi]._u * (fa[fi]._pF->GetUV(1) - fa[fi]._pF->GetUV(0)) +
			fa[fi]._v * (fa[fi]._pF->GetUV(2) - fa[fi]._pF->GetUV(0)) +
			fa[fi]._pF->GetUV(0);

		Data.Vertex.TextureName = fa[fi]._pF->GetTextureName();
		Data.Vertex.Color = scene.GetTextureColor(Data.Vertex.TextureName, Data.Vertex.UV);
		if (Data.Vertex.Color.x() < 0){
			Data.Vertex.Color = fa[fi]._pF->GetDiffuseColor();
		}

		Data.Vertex.Emit = (fa[fi]._pF->GetEmission());


		float N_from = Data.Ray.CRay.m_Index;
		float N_to = 1.0f;


		float russian_p = 1.0f;
		if (Data.Depth > kDepth){
			russian_p = std::max(
				std::max(Data.Vertex.Color[0], Data.Vertex.Color[1]), Data.Vertex.Color[2]) *
				pow(0.4, Data.Depth - kDepth);

			if (russian_p < MT.genrand64_real1()){
				Data.IfLoop = false;
				return;
			}
		}

		//VPLで用いるパラメータの計算・取得
		/////
		float g = m_agg->GetPhase_g();
		Vec3 Sigma_t = m_agg->GetSigma_t();
		float at = 0.3333333 * (Sigma_t[0] + Sigma_t[1] + Sigma_t[2]);
		Vec3 Alpha_s = m_agg->GetAlpha_s();
		Vec3 Sigma_a = Vec3::Zero();
		/////


		//直接指定するという選択
		/////////////////////////////////////////////////////////////////////////////
		//float sc = 10.0f;
		//Alpha_s = Vec3(0.83, 0.79, 0.75);
		//Vec3 sigma_s = Vec3(2.19, 2.62, 3.00) * sc;
		//Sigma_a = Vec3(0.0021, 0.0041, 0.0071) * sc;
		//Sigma_t = Sigma_a + sigma_s;
		/////////////////////////////////////////////////////////////////////////////



		if (Data.Vertex.Normal.dot(Data.Ray.CRay.m_Dir) > 0){
			//Rayが物体内部

			float FL = fa[fi]._t;
			float uu = MT.genrand64_real3();
			float PSL = -(1.0f / at)*log(uu + 0.000001);

			if (FL > PSL){

				Data.Vertex.Pos = Data.Ray.CRay.m_Org + PSL * Data.Ray.CRay.m_Dir;
				float p = ISample::NextDirSphere(Data.Ray.CRay, Data.Vertex.Normal, Data.Ray.NextRay.m_Dir, MT);
				p = 1.0f / (4.0f * M_PI);

				Vec3 ww =
					Alpha_s *
					m_agg->PhaseF(Data.Ray.CRay.m_Dir, Data.Ray.NextRay.m_Dir) / (russian_p * p);

				MultVec(Vec3(exp(-Sigma_a[0] * PSL), exp(-Sigma_a[1] * PSL), exp(-Sigma_a[2] * PSL)), Data.Vertex.weight);

				Data.Ray.NextRay.m_Org = Data.Vertex.Pos;
				//MultVec(ww, Data.Vertex.weight);
				MultVec(Alpha_s, Data.Vertex.weight);


				////DAへの切り替え処理
				if (scene.IfDeeperThanDepth(Data.Vertex.Pos, 1.0f / at)){
					m_Mode = DA;
					Data.Diffuse.Pos = Data.Vertex.Pos;
					return;
				}

			}
			else{
				Vec3 ww = Vec3::Ones() / (russian_p);
				MultVec(ww, Data.Vertex.weight);
				MultVec(Vec3(exp(-Sigma_a[0] * FL), exp(-Sigma_a[1] * FL), exp(-Sigma_a[2] * FL)), Data.Vertex.weight);
				Data.Vertex.Pos = (fa[fi]._t + 0.000001) * Data.Ray.CRay.m_Dir + Data.Ray.CRay.m_Org;
				Data.Ray.NextRay.m_Dir = Data.Ray.CRay.m_Dir;
				Data.Ray.NextRay.m_Org = Data.Vertex.Pos;
			}
		}
		else{

			Ray rrr = Data.Ray.CRay;
			float R = 0.027;

			if (MT.genrand64_real1() < R){
			}
			else{
				//透過
				Data.Vertex.Normal = -Data.Vertex.Normal;
			}
			ISample::NextDirDiffuse(Data.Ray.CRay, Data.Vertex.Normal, rrr.m_Dir, MT);

			Data.Vertex.Pos = (fa[fi]._t) * Data.Ray.CRay.m_Dir + Data.Ray.CRay.m_Org;
			Data.Ray.NextRay.m_Dir = rrr.m_Dir;
			//Data.Ray.NextRay.m_Dir = Data.Ray.CRay.m_Dir;
			Data.Ray.NextRay.m_Org = Data.Vertex.Pos + Data.Ray.NextRay.m_Dir * 0.0001;

			MultVec(Alpha_s / (russian_p), Data.Vertex.weight);

		}



	}


	//Explicit Path Tracing　つぶつぶ描画編
	void EPT_Intra(
		PTDATA& Data,
		const Scene& scene,
		const Camera& camera,
		RandomMT& MT){

		float russian_p = 1.0f;
		if (Data.Depth > kDepth){
			russian_p = std::max(std::max(Data.Vertex.Color[0], Data.Vertex.Color[1]), Data.Vertex.Color[2]) * pow(0.4, Data.Depth - kDepth);
			if (russian_p < MT.genrand64_real1()){
				Data.IfLoop = false;
				return;
			}
		}

		if (m_agg->NextRay(Data.EPTintra.GrainIDX, Data.Ray.CRay, Data.EPTintra.GCenter, Data.EPTintra.PID, Data.Ray.NextRay, Data.Pdf.Pdf, MT)){
			Data.Vertex.EPT_Depth++;
			Data.Flag.IFInObject = false;
			MultVec(Data.Vertex.weight, Data.Vertex.Emit);
			Data.Result += Data.Vertex.Emit / russian_p;

			MultVec(m_agg->GetGrain(Data.EPTintra.GrainIDX)->GetColor() / (Data.Pdf.Pdf), Data.Vertex.weight);
		}
		else{
			m_Mode = EPT_INTER;
		}

	}


	//Explicit Path Tracing　つぶつぶ間のトラバーサル
	void EPT_Inter(
		PTDATA& Data,
		const Scene& scene,
		const Camera& camera,
		RandomMT& MT){

		//////とりあえず
		//if (Data.Vertex.EPT_Depth > 30){
		//	m_Mode = VPT;
		//	Data.Depth--;
		//	//Data.IfLoop = false;
		//	return;
		//}

		std::vector<const VoxelGrid::VPP> vv;
		if (m_agg->Intersect(Data.Ray.CRay, vv)){

			if (vv[0]._pp->_Vstat > 0){

				//ロシアのしきたりにより、追跡を打ち切る
				float russian_p = 1.0f;
				if (Data.Depth > kDepth){
					russian_p = std::max(std::max(Data.Vertex.Color[0], Data.Vertex.Color[1]), Data.Vertex.Color[2]) * pow(0.4, Data.Depth - kDepth);
					if (russian_p < MT.genrand64_real1()){
						Data.IfLoop = false;
						return;
					}
				}

				Data.Vertex.TextureName = "Agg";
				if (Data.Vertex.TextureName == "Agg"){
					//つぶつぶ材質

					Data.Vertex.Color = Vec3::Ones();

					Vec3 Center = m_agg->GetCenter(vv[0]._pp);

					Data.Vertex.Pos = (vv[0]._f) * Data.Ray.CRay.m_Dir + Data.Ray.CRay.m_Org;
					Data.Ray.NextRay.m_Dir = Data.Ray.CRay.m_Dir;
					Data.Ray.NextRay.m_Org = Data.Vertex.Pos;

					std::vector<Aggrigate::SC> rs;
					if (m_agg->GetSphereFromRay(Center, Data.Ray.CRay, rs)){


						int CurrS = -1;
						if (vv[0]._pp->_Vstat != 100){
							//Partially Inside
							for (int ss = 0; ss < rs.size(); ss++){

								std::vector<AbsG::RF> fa;
								Vec3 RD;
								ISample::NextDirSphere(Data.Ray.CRay, Data.Vertex.Normal, RD, MT);
								int nf = scene.GetFacetFromRay(Ray(rs[ss]._Pos, RD), fa);
								if (nf > 0){
									if ((fa[0]._pF->GetNormal(0) +
										fa[0]._pF->GetNormal(1) +
										fa[0]._pF->GetNormal(2)).normalized().dot(RD) > 0.001){

										CurrS = ss;
										break;
									}
								}



							}
						}
						else{
							//Fully Inside
							CurrS = 0;
						}

						if (CurrS != -1){

							Data.Ray.NextRay.m_Org = Data.Ray.CRay.m_Org + (rs[CurrS]._t) * Data.Ray.CRay.m_Dir;
							Data.Ray.NextRay.m_Dir = Data.Ray.CRay.m_Dir;

							Data.Ray.NextRay.m_Org += 0.000001 * Data.Ray.NextRay.m_Dir;
							Data.EPTintra.PID =
								(vv[0]._pp->_ID[0] + vv[0]._pp->_ID[1] * 15 + vv[0]._pp->_ID[2] * 225 + (rs[CurrS]._SID + 112) * 4321);
							Data.EPTintra.GCenter = rs[CurrS]._Pos;

							Data.EPTintra.GrainIDX = m_agg->RandomSelectGrain(Data.EPTintra.PID);

							m_Mode = EPT_INTRA;



						}
					}
					//Data.Depth--;

					/*Data.Flag.IFInObject = false;
					MultVec(Data.Vertex.weight, Data.Vertex.Emit);
					Data.Result += Data.Vertex.Emit;

					MultVec(Data.Vertex.Color / (russian_p * Data.Pdf.Pdf), Data.Vertex.weight);*/


				}
				else{
					//拡散面

					//d_Color = Vec3::Zero();


					///*
					//NextEvent = true;
					//NextEventEstimate(
					//std::max(0, 0 - depth),
					//scene,
					//MT,
					//Pos,
					//Normal,
					//r,
					//Scene::DiffuseBRDF,
					//d_Color);

					//if (m_IBL != NULL){
					//NextEventIBL(
					//std::max(0, 0 - depth),
					//scene,
					//MT,
					//Pos,
					//Normal,
					//r,
					//Scene::DiffuseBRDF,
					//ibl_Color);
					//}
					//*/

					//pdf = ISample::NextDirDiffuse(CRay, Normal, NextRay.m_Dir, MT);
					//NextRay.m_Org = Pos + 0.00000001 * NextRay.m_Dir;

					//IFInObject = false;
					//MultVec(weight, Emit);
					//result += Emit;

					//MultVec(Color / (russian_p * pdf), weight);
				}

			}
			else{
				Data.Vertex.Pos = (vv[0]._f) * Data.Ray.CRay.m_Dir + Data.Ray.CRay.m_Org;
				Data.Ray.NextRay.m_Dir = Data.Ray.CRay.m_Dir;
				Data.Ray.NextRay.m_Org = Data.Vertex.Pos;
				//Data.Depth--;
			}

		}
		else{
			//シーンと交差なし
			//IBLを寄与に加えてターンエンド
			Data.Vertex.Emit = m_IBL->GetColor1(Data.Ray.CRay.m_Dir);

			MultVec(Data.Vertex.weight, Data.Vertex.Emit);
			Data.Result += Data.Vertex.Emit;

			Data.IfLoop = false;
			return;
		}
	}



	Vec3 RenderSceneR(
		const Ray& r,
		int depth,
		const Scene& scene,
		const Camera& camera,
		RandomMT& MT,
		bool InObject){


		m_Mode = EPT_INTER;
		//m_Mode = VPT;

		PTDATA Data;
		Data.Ray.CRay = r;


		while (Data.IfLoop){


			switch (m_Mode){

			case EPT_INTER:
				///////////////////////////////////////////////////////////////////////////
				///////////////////////////////////////////////////////////////////////////

				EPT_Inter(Data, scene, camera, MT);

				///////////////////////////////////////////////////////////////////////////
				///////////////////////////////////////////////////////////////////////////
				break;


			case EPT_INTRA:
				///////////////////////////////////////////////////////////////////////////
				///////////////////////////////////////////////////////////////////////////

				EPT_Intra(Data, scene, camera, MT);
				///////////////////////////////////////////////////////////////////////////
				///////////////////////////////////////////////////////////////////////////
				break;


			case VPT:
				///////////////////////////////////////////////////////////////////////////
				///////////////////////////////////////////////////////////////////////////
				VolumePathTrace(Data, scene, camera, MT);

				///////////////////////////////////////////////////////////////////////////
				///////////////////////////////////////////////////////////////////////////
				break;


			case DA:
				///////////////////////////////////////////////////////////////////////////
				///////////////////////////////////////////////////////////////////////////
				DiffusionAproximate_Dipole(Data, scene, camera, MT);

				///////////////////////////////////////////////////////////////////////////
				///////////////////////////////////////////////////////////////////////////
				break;



			}
			Data.Ray.CRay = Data.Ray.NextRay;
			Data.Ray.CRay.m_Dir.normalize();
			Data.Depth++;

		}


		return Data.Result;
	}



	Vec3 RenderScene(
		const Ray& r,
		int depth,
		const Scene& scene,
		const Camera& camera,
		RandomMT& MT,
		bool InObject,
		bool NE){


		int Depth = 0;
		bool NextEvent = false;
		bool IFInObject = false;
		float ToIndex = 1.0f;
		Ray NextRay = r;
		Ray CRay = r;

		Vec3 Normal(0, 1, 0);
		Vec3 Pos(0, 0, 0);
		Vec2 UV(0, 0);
		std::string TextureName;
		Vec3 Color(0, 1, 0);
		Vec3 Emit(0, 0, 0);
		Vec3 in_Color(0, 0, 0);
		Vec3 d_Color(0, 0, 0);
		Vec3 ibl_Color(0, 0, 0);
		Vec3 weight = Vec3::Ones();
		Vec3 NN = Normal;
		float pdf = 1.0f;

		//最終結果
		Vec3 result(0, 0, 0);

		while (true){

			std::vector<AbsG::RF> fa;
			int nf = scene.GetFacetFromRay(CRay, fa);

			if (nf > 0){
				//シーンと交差あり


				//有効な面のうち、最も交差が早い面を見つける
				int fi = -1;
				for (int ffi = 0; ffi < nf; ffi++){

					Normal = fa[ffi]._u * (fa[ffi]._pF->GetNormal(1) - fa[ffi]._pF->GetNormal(0)) +
						fa[ffi]._v * (fa[ffi]._pF->GetNormal(2) - fa[ffi]._pF->GetNormal(0)) +
						fa[ffi]._pF->GetNormal(0);
					Normal.normalize();

					if (fa[ffi]._t > 0.000001){
						if (IFInObject){
							fi = ffi;
							break;
						}
						else{
							if (Normal.dot(CRay.m_Dir) < -0.00000001){
								fi = ffi;
								break;
							}
						}
					}
				}
				if (fi == -1){ return Vec3::Zero(); }



				//必ず視線と反対方向を向く
				NN = Normal;
				if (Normal.dot(-CRay.m_Dir) < 0.0f){
					NN = -Normal;
				}

				//今いる位置
				Pos = fa[fi]._t * CRay.m_Dir + CRay.m_Org;

				//テクスチャ座標
				UV =
					fa[fi]._u * (fa[fi]._pF->GetUV(1) - fa[fi]._pF->GetUV(0)) +
					fa[fi]._v * (fa[fi]._pF->GetUV(2) - fa[fi]._pF->GetUV(0)) +
					fa[fi]._pF->GetUV(0);

				//テクスチャファイル名。。。と思いきや、マテリアルの名前
				TextureName = fa[fi]._pF->GetTextureName();

				//表面の色
				Color = scene.GetTextureColor(TextureName, UV);
				if (Color.x() < 0){
					//テクスチャが無い場合は、Diffuseの色を使う
					Color = fa[fi]._pF->GetDiffuseColor();
				}

				//光ってる度合い　光ってないならゼロ
				Emit = (fa[fi]._pF->GetEmission());


				//ロシアのしきたりにより、追跡を打ち切る
				float russian_p = 1.0f;
				if (Depth > kDepth){
					russian_p = std::max(std::max(Color[0], Color[1]), Color[2]) * pow(0.4, Depth - kDepth);
					if (russian_p < MT.genrand64_real1()){
						return result;
					}
				}


				if (TextureName == "Agg"){
					////つぶつぶ材質


					//std::vector<const VoxelGrid::VoxelData*> vv;
					//if (m_agg->Intersect(r, vv)){

					//	RandomMT mm(vv[0]->_ID[0] + vv[0]->_ID[1] * 100 + vv[0]->_ID[2] * 10000);
					//	float c1 = mm.genrand64_real1();
					//	float c2;
					//	float c3;

					//	return Vec3(c1, 1 - c1, c1*c1);
					//	//return Vec3(c1, c2, c3);
					//}
					//else{
					//	return Vec3::Zero();
					//}


				}
				else{
					//拡散面

					d_Color = Vec3::Zero();


					/*
					NextEvent = true;
					NextEventEstimate(
					std::max(0, 0 - depth),
					scene,
					MT,
					Pos,
					Normal,
					r,
					Scene::DiffuseBRDF,
					d_Color);

					if (m_IBL != NULL){
					NextEventIBL(
					std::max(0, 0 - depth),
					scene,
					MT,
					Pos,
					Normal,
					r,
					Scene::DiffuseBRDF,
					ibl_Color);
					}
					*/

					pdf = ISample::NextDirDiffuse(CRay, Normal, NextRay.m_Dir, MT);
					NextRay.m_Org = Pos + 0.00000001 * NextRay.m_Dir;

					IFInObject = false;
					MultVec(weight, Emit);
					result += Emit;

					MultVec(Color / (russian_p * pdf), weight);
				}

			}
			else{
				//シーンと交差なし
				//IBLを寄与に加えてターンエンド
				Emit = m_IBL->GetColor1(CRay.m_Dir);

				MultVec(weight, Emit);
				result += Emit;

				return result;
			}

			CRay = NextRay;
			Depth++;

		}


		return result;
	}

	float NextEventEstimate(
		int NShadowRay,
		const Scene& scene,
		RandomMT& MT,
		const Vec3& Pos,
		const Vec3& Normal,
		const Ray& ray,
		BRDF brdf,
		Vec3& result){

		result = Vec3::Zero();
		if (!scene.IfExistLight()){
			return 0.0f;
		}
		if (NShadowRay <= 0){
			return 0.0f;
		}

		int SpL = NShadowRay;

		for (int l = 0; l < SpL; l++){
			float pp;
			Scene::LightMesh LM = scene.GetLightMesh(MT, pp);
			float Lu, Lv;
			Lu = MT.genrand64_real1();
			Lv = MT.genrand64_real1();
			if (Lu + Lv > 1){
				Lu = 1.0f - Lu;
				Lv = 1.0f - Lv;
			}
			Vec3 LP = LM._U *Lu + LM._V * Lv + LM._O;
			if (LM._Normal.dot(LP - Pos) > 0){
				continue;
			}

			//光源から明示的にサンプル///////
			std::vector<AbsG::RF> ffa;
			scene.GetFacetFromRay(Ray(Pos + 0.0001 * Normal, LP - Pos), ffa);
			float NdL = std::max(0.0f, std::min(1.0f, (LP - Pos).normalized().dot(Normal)));
			float NdNL = std::max(0.0f, std::min(1.0f, (LP - Pos).normalized().dot(-LM._Normal)));

			if (ffa.size() > 0){
				if (ffa[0]._t + 0.005 > (LP - Pos).norm()){
					float G = NdL * NdNL / pow(ffa[0]._t, 2.0f);
					result +=
						brdf((Pos - LP).normalized(), -ray.m_Dir) *
						LM._Color * LM._Power *
						G *
						scene.GetSumLightArea() / (float)SpL;
				}
			}


		}

		return 0.0f;
	}


	float NextEventIBL(
		int NShadowRay,
		const Scene& scene,
		RandomMT& MT,
		const Vec3& Pos,
		const Vec3& Normal,
		const Ray& ray,
		BRDF brdf,
		Vec3& result){

		if (NShadowRay <= 0){
			return 0.0f;
		}

		int SpL = NShadowRay;
		result = Vec3::Zero();

		for (int l = 0; l < SpL; l++){

			Vec3 ndir;
			ISample::NextDirDiffuse(ray, Normal, ndir, MT);
			//ISample::NextDirHemSphere(ray, Normal, ndir, MT);

			//光源から明示的にサンプル///////
			std::vector<AbsG::RF> ffa;
			scene.GetFacetFromRay(Ray(Pos + 0.0001 * Normal, ndir), ffa);

			if (ffa.size() == 0){
				result += m_IBL->GetColor2(ndir);
			}


		}
		result /= ((float)SpL);

		return 0.0f;
	}


	void MultVec(const Vec3& a, Vec3& b){
		b = Vec3(b[0] * a[0], b[1] * a[1], b[2] * a[2]);
	}

};












GPathTrace::GPathTrace(int IMGwidth, int IMGheight) : m_Core(new GPTCore()), iRender(IMGwidth, IMGheight){

}

GPathTrace::~GPathTrace(){

}


void GPathTrace::RenderScene(const Scene& scene, const Camera& camera, const iIBL* ibl, int sample, int s_sample){

	m_Core->m_IBL = ibl;


	m_Core->m_agg = new Aggrigate(1.0f);
	//m_Core->m_agg->InitAggrigate("Voxel\\vvv.obj");
	m_Core->m_agg->InitAggrigate(scene, 3);


	m_Core->GetRd_DA(GPTCore::PTDATA(), Vec3::Zero(), Vec3(0, 1, 0), Vec3::Zero(), Vec3(0, 1, 0));



	delete[] m_Core->img;
	m_Core->img = new char[GetIMGWidth()*GetIMGHeight() * 4];
	for (int i = 0; i < 4 * GetIMGHeight()*GetIMGWidth(); i++){
		m_Core->img[i] = 0;
	}


	time_t st = time(NULL);

	Vec3 ssdx = 0.5f * camera._sx / (float)GetIMGWidth();
	Vec3 ssdy = 0.5f * camera._sy / (float)GetIMGHeight();








	/*std::ofstream file;
	file.open("aaa.obj");
	RandomMT mm(11);
	for (int i = 0; i < 50000; i++){

	float x, y, z;
	x = mm.genrand64_real1()*2.0f - 1.0f;
	y = mm.genrand64_real1()*2.0f - 1.0f;
	z = mm.genrand64_real1()*2.0f - 1.0f;

	bool Inter = true;
	for (int i = 0; i < 3; i++){
	std::vector<AbsG::RF> fa;
	Vec3 RD;
	ISample::NextDirSphere(
	Ray(Vec3(0, 1, 0), Vec3(x, y, z)),
	Vec3(0, 1, 0), RD, mm);

	int nf = scene.GetFacetFromRay(Ray(Vec3(x, y, z), RD), fa);

	if (nf > 0){
	if (
	(
	fa[0]._u * (fa[0]._pF->GetNormal(1) - fa[0]._pF->GetNormal(0)) +
	fa[0]._v * (fa[0]._pF->GetNormal(2) - fa[0]._pF->GetNormal(0)) +
	fa[0]._pF->GetNormal(0)
	).dot(RD) < 0.001){
	Inter = false;
	break;
	}
	}
	else{
	Inter = false;
	break;
	}
	}
	if (Inter){
	if (scene.IfDeeperThanDepth(Vec3(x, y, z), 0.1)){
	file << "v " << x << " " << y << " " << z << std::endl;
	}
	}

	}
	file.close();*/


#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
	for (int i = 0; i < (GetIMGWidth()) * (GetIMGHeight()); i++){

		Vec3 Color(0, 0, 0);

		RandomMT MT(i + time(NULL) * 1937);



		int w = i % GetIMGWidth();
		int h = i / GetIMGWidth();

		for (int ssx = 0; ssx < s_sample; ssx++){
			for (int ssy = 0; ssy < s_sample; ssy++){

				Vec3 RV(MT.genrand64_real1() * ssdx + MT.genrand64_real1() * ssdy);

				Ray r(
					camera._CameraPos,
					camera._Center + RV +
					camera._sx * (w - 0.5 * GetIMGWidth()) / (float)GetIMGWidth() +
					camera._sy * (0.5 * GetIMGHeight() - h) / (float)GetIMGHeight() - camera._CameraPos);

				//r = Ray(Vec3(2.0f*w / (float)GetIMGWidth() - 1.0f, 2.0f * h / (float)GetIMGHeight() - 1.0f, 5), Vec3(-0.01, 0.0, -1));

				for (int s = 0; s < sample; s++){
					//Color += m_Core->RenderScene(r, 0, scene, camera, MT, false, false);
					Color += m_Core->RenderSceneR(r, 0, scene, camera, MT, false);
				}

			}
		}
		Color /= (float)(s_sample*s_sample*sample);

		m_Core->img[4 * i + 0] = std::max(std::min(0.999f, Color.x()), 0.0f) * 255;
		m_Core->img[4 * i + 1] = std::max(std::min(0.999f, Color.y()), 0.0f) * 255;
		m_Core->img[4 * i + 2] = std::max(std::min(0.999f, Color.z()), 0.0f) * 255;
		m_Core->img[4 * i + 3] = 255;


		fprintf(stderr, "\r つぶつぶー %5.2f%%", 100.0*i / ((float)GetIMGWidth() * GetIMGHeight()));

	}
	printf("\nしゅーりょ \n");


	printf("\n");
	time_t ed = time(NULL);
	int rtm = int(ed) - int(st);
	int th, tm, ts;
	th = rtm / 3600;
	tm = (rtm - th * 3600) / 60;
	ts = rtm - 3600 * th - 60 * tm;

	printf("経過　%d時間 %d分　%d秒\n", th, tm, ts);












	delete m_Core->m_agg;

}

void GPathTrace::WriteImage(const char* FileName){

	char* img;
	img = new char[GetIMGWidth()*GetIMGHeight() * 4];

	// //レインボーで初期化　ふつくしい
	//for (int h = 0; h < GetIMGHeight(); h++){
	//	for (int w = 0; w < GetIMGWidth(); w++){
	//		img[4 * (w + h*GetIMGWidth()) + 0] = 255 * h/(float)GetIMGHeight();
	//		img[4 * (w + h*GetIMGWidth()) + 1] = 255 * w/(float)GetIMGWidth();
	//		img[4 * (w + h*GetIMGWidth()) + 2] = 255;
	//		img[4 * (w + h*GetIMGWidth()) + 3] = 255;
	//	}
	//}



	for (int h = 0; h < GetIMGHeight(); h++){
		for (int w = 0; w < GetIMGWidth(); w++){
			img[4 * (w + (GetIMGHeight() - h - 1)*GetIMGWidth()) + 0] = m_Core->img[4 * (w + h*GetIMGWidth()) + 0];
			img[4 * (w + (GetIMGHeight() - h - 1)*GetIMGWidth()) + 1] = m_Core->img[4 * (w + h*GetIMGWidth()) + 1];
			img[4 * (w + (GetIMGHeight() - h - 1)*GetIMGWidth()) + 2] = m_Core->img[4 * (w + h*GetIMGWidth()) + 2];
			img[4 * (w + (GetIMGHeight() - h - 1)*GetIMGWidth()) + 3] = m_Core->img[4 * (w + h*GetIMGWidth()) + 3];
		}
	}



	//for (int h = 0; h < GetIMGHeight(); h++){
	//	for (int w = 0; w < GetIMGWidth(); w++){
	//		img[4 * (w + h*GetIMGWidth()) + 0] = 255 * h/(float)GetIMGHeight();
	//		img[4 * (w + h*GetIMGWidth()) + 1] = 255 * w/(float)GetIMGWidth();
	//		img[4 * (w + h*GetIMGWidth()) + 2] = 255;
	//		img[4 * (w + h*GetIMGWidth()) + 3] = 255;
	//	}
	//}


	stbi_write_bmp(FileName, GetIMGWidth(), GetIMGHeight(), 4, (char*)img);

	delete[] img;

}