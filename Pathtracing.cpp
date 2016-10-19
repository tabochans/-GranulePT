#include"Pathtracing.h"

using namespace PTUtility;

class Pathtracing::PTCore{

public:

	iIBL const *m_IBL;

	typedef float(*BRDF)(const Vec3& in, const Vec3& out);

	const int kDepth;

	float* img;

	PTCore() : kDepth(30){
	}
	~PTCore(){
		delete[] img;
	}

	Vec3 RenderSceneR(
		const Ray& r,
		int depth,
		const Scene& scene, 
		const Camera& camera,
		RandomMT& MT,
		bool InObject){

		
		return Vec3::Zero();
	}



	Vec3 RenderScene(
		const Ray& r, 
		int depth,
		const Scene& scene, 
		const Camera& camera,
		RandomMT& MT,
		bool InObject,
		bool NE){

		std::vector<AbsG::RF> fa;
		int nf = scene.GetFacetFromRay(r, fa);

		//SSSのための何か
		float L = 0.05;
		float Ac = 2.5;
		Vec3 SSSc(1, 1, 1);
		bool NextEvent = true;
		bool InObj = InObject;
		Vec3 Normal(0, 1, 0);

		if (nf > 0){

			int fi = -1;
			for (int ffi = 0; ffi < nf; ffi++){

				Normal = fa[ffi]._u * (fa[ffi]._pF->GetNormal(1) - fa[ffi]._pF->GetNormal(0)) +
					fa[ffi]._v * (fa[ffi]._pF->GetNormal(2) - fa[ffi]._pF->GetNormal(0)) +
					fa[ffi]._pF->GetNormal(0);
				Normal.normalize();

				if (fa[ffi]._t > 0.001){
					if (r.m_Index < 1.05){
						if (Normal.dot(-r.m_Dir) > 0.00001f){
							fi = ffi;
							break;
						}
					}
					else if (r.m_Index > 1.1){
						fi = ffi;
						break;
					}
				}
			}
			if (fi == -1){ return Vec3::Zero(); }

			Vec3 Pos(0, 0, 0);
			Vec2 UV(0, 0);
			std::string TextureName;
			Vec3 Color(0, 1, 0);
			Vec3 Emit(0, 0, 0);
			Vec3 in_Color(0, 0, 0);
			Vec3 d_Color(0, 0, 0);
			Vec3 ibl_Color(0, 0, 0);
			Vec3 weight = Vec3::Ones();
			float N;
			float p = 1.0f;
			float fr = 1.0f;

			//必ず視線と反対方向を向く
			Vec3 NN = Normal;
			if (Normal.dot(-r.m_Dir) < 0.0f){
				NN = -Normal;
			}

			Pos = fa[fi]._t * r.m_Dir + r.m_Org;

			UV =
				fa[fi]._u * (fa[fi]._pF->GetUV(1) - fa[fi]._pF->GetUV(0)) +
				fa[fi]._v * (fa[fi]._pF->GetUV(2) - fa[fi]._pF->GetUV(0)) +
				fa[fi]._pF->GetUV(0);



			TextureName = fa[fi]._pF->GetTextureName();
			Color = scene.GetTextureColor(TextureName, UV);
			if (Color.x() < 0){
				Color = fa[fi]._pF->GetDiffuseColor();
			}

			Emit = (fa[fi]._pF->GetEmission());
			

			float N_from = r.m_Index;
			float N_to = 1.0f;


			float russian_p = 1.0f;
			if (depth > kDepth){
				russian_p = std::max(std::max(Color[0], Color[1]), Color[2]) * pow(0.4, depth - kDepth);
				if (russian_p < MT.genrand64_real1()){
					//return Vec3(0, 0, 1);
					return Vec3::Zero();
					return Emit;
				}
			}

			
			////////////////////////////////
					
			Vec3 nDir;



			//マテリアルで処理を場合分け～
			//結局レンダラの方にかいちゃってるじゃねーか
			if (fa[fi]._pF->GetTextureName() == "tile"){
				//Phong

				float p = ISample::NextDirPhong(r, Normal, nDir, MT, 100);
				if (nDir.dot(Normal) < 0.001){
					return Vec3::Zero();
				}

				in_Color += RenderScene(Ray(Pos + 0.0001 * nDir, nDir, N_to), depth + 1, scene, camera, MT, InObj, false);

				weight = Color / (russian_p * p);
				MultVec(weight, in_Color);


			}
			else if (fa[fi]._pF->GetTextureName() == "Default" || fa[fi]._pF->GetTextureName() == "default.001"){
				//ポット
				d_Color = Vec3(0, 0, 0);
				bool nexx = false;
				float p = 1.0f;
				if ((float)Color.x() + 0.05 < (float)Color.z()){
					p = ISample::NextDirSpecular(r, Normal, nDir, MT);
				}
				else{

					/*NextEventEstimate(
						std::max(0, 0 - depth),
						scene,
						MT,
						Pos,
						Normal,
						r,
						Scene::DiffuseBRDF,
						d_Color);*/

					p = ISample::NextDirDiffuse(r, Normal, nDir, MT);
					//nexx = true;
				}

				in_Color += RenderScene(Ray(Pos + 0.0001 * nDir, nDir, N_to), depth + 1, scene, camera, MT, InObj, nexx);

				weight = Color / (russian_p * p);
				in_Color += d_Color;
				MultVec(weight, in_Color);

			}
			else if (fa[fi]._pF->GetTextureName() == "LPB"){
				//Phong

				float p = ISample::NextDirPhong(r, Normal, nDir, MT, 10);
				if (nDir.dot(Normal) < 0.001){
					return Vec3::Zero();
				}

				in_Color += RenderScene(Ray(Pos + 0.0001 * nDir, nDir, N_to), depth + 1, scene, camera, MT, InObj, false);

				weight = Color / (russian_p * p);
				MultVec(weight, in_Color);


			}
			else if (fa[fi]._pF->GetTextureName() == "mateal"){
				///////////////////////
				//ファンタスティックな金属
				//よくわかんねーけどファンタスティック

				//ファンタスティックパラメータ
				//増やせば増やすほど、ファンタスティックな金属になる
				float rgp = 0.75f;

				
				if (MT.genrand64_real1() > rgp){
					//ノンファンタスティック
					float p = ISample::NextDirSpecular(r, Normal, nDir, MT);

					in_Color += RenderScene(Ray(Pos + 0.0001 * nDir, nDir, N_to), depth + 1, scene, camera, MT, InObj, false);

					weight = Color / (russian_p * p);
					MultVec(weight, in_Color);
				}
				else{
					//ジャストファンタスティック
					NextEventEstimate(
						std::max(0, 5 - depth),
						scene,
						MT,
						Pos,
						Normal,
						r,
						Scene::DiffuseBRDF,
						d_Color);
					d_Color *= M_PI;

					p = ISample::NextDirDiffuse(r, Normal, nDir, MT);
					NextEvent = true;

					in_Color = d_Color;
					in_Color += RenderScene(Ray(Pos + 0.0001 * nDir, nDir, N_to), depth + 1, scene, camera, MT, InObj, true);

					weight = Color / (russian_p * p);
					MultVec(weight, in_Color);

				}

				


			}
			else if (fa[fi]._pF->GetTextureName() == "LCY"){
				///////////////////////
				//ガラス
				Color = Vec3::Ones() * 0.95;
				if (N_from < 1.05){
					N_to = 1.77;
				}
				else{
					N_to = 1.0;
				}
				fr = ISample::NextDirRefraction(r, Normal, N_from, N_to, nDir, MT);

				in_Color += fr * RenderScene(Ray(Pos + 0.0001 * nDir, nDir, N_to), depth + 1, scene, camera, MT, InObj, false);

				weight = Color / (russian_p * p);
				MultVec(weight, in_Color);


			}
			else if (fa[fi]._pF->GetTextureName() == "sdb"){
				///////////////////////
				//金属

				float p = ISample::NextDirSpecular(r, Normal, nDir, MT);

				in_Color += RenderScene(Ray(Pos + 0.0001 * nDir, nDir, N_to), depth + 1, scene, camera, MT, InObj, false);

				weight = Color / (russian_p * p);
				MultVec(weight, in_Color);

			}
			else if (fa[fi]._pF->GetTextureName() == "mat"){
				//屈折付きはんとーめ
				//なんかびみょー　うまくいってない

				Vec3 alpha_s = Vec3(0.96, 0.86, 0.72);
				Vec3 alpha_t = Vec3(2.7, 3.1, 2.9) * 3.0;


				if (InObject){

					float FL = fa[fi]._t;
					float uu = MT.genrand64_real3() * 0.96 + 0.02;
					float PSL = 8.0f * L * pow(uu - 0.5, 3) + L;


					if (FL > PSL){

						Pos = r.m_Org + PSL * r.m_Dir;
						p = ISample::NextDirSphere(r, Normal, nDir, MT);
						InObj = true;
						N_to = 2.0f;

						in_Color += RenderScene(Ray(Pos + 0.0001 * nDir, nDir, N_to), depth + 1, scene, camera, MT, InObj, false);

						weight = Vec3::Ones() / (russian_p * p);
						MultVec(weight, in_Color);
						MultVec(Vec3(exp(-alpha_t[0] * PSL) * alpha_s[0], exp(-alpha_t[1] * PSL) * alpha_s[1], exp(-alpha_t[2] * PSL) * alpha_s[2]), in_Color);

					}
					else{

						N_to = 1.0f;
						N_from = 2.0f;

						fr = ISample::NextDirRefraction(r, Normal, N_from, N_to, nDir, MT);
						if (nDir.dot(Normal) > 0){
							InObj = false;
							N_to = 1.0f;
						}
						else{
							InObj = true;
							N_to = 2.0f;
						}

						in_Color += fr * RenderScene(Ray(Pos + 0.0001 * nDir, nDir, N_to), depth + 1, scene, camera, MT, InObj, false);

						weight = Color / (russian_p * p);
						MultVec(weight, in_Color);
						MultVec(Vec3(exp(-alpha_t[0] * FL) * alpha_s[0], exp(-alpha_t[1] * FL) * alpha_s[1], exp(-alpha_t[2] * FL) * alpha_s[2]), in_Color);

					}


				}
				else{

					N_to = 2.0f;
					N_from = 1.0f;

					fr = ISample::NextDirRefraction(r, Normal, N_from, N_to, nDir, MT);
					


					if (nDir.dot(Normal) < 0){

						InObj = true;
						N_to = 2.0f;
						
						in_Color += fr * RenderScene(Ray(Pos + 0.000001 * nDir, nDir, N_to), depth + 1, scene, camera, MT, InObj, false);
						
					}
					else{
						
						InObj = false;
						N_to = 1.0f;
						in_Color += fr * RenderScene(Ray(Pos + 0.000001 * nDir, nDir, N_to), depth + 1, scene, camera, MT, InObj, false);

					}

					weight = Color / (russian_p * p);
					MultVec(weight, in_Color);

				}

			}
			else if (fa[fi]._pF->GetTextureName() == "Maaat"){

				float sc = 12.0f;
				Vec3 alpha_s = Vec3(0.83, 0.79, 0.75);
				Vec3 sigma_s = Vec3(2.19, 2.62, 3.00) * sc;
				Vec3 sigma_a = Vec3(0.0021, 0.0041, 0.0071) * sc;

				Vec3 sigma_t = sigma_a + sigma_s;

				float at = 0.333 * (sigma_t[0] + sigma_t[1] + sigma_t[2]);

				//速攻魔法拡散する半透明物体

				if (InObject){

					float FL = fa[fi]._t;
					float uu = MT.genrand64_real3();
					float PSL = -(1.0f / at)*log(uu + 0.000001);

					if (FL > PSL){

						Pos = r.m_Org + PSL * r.m_Dir;
						p = ISample::NextDirSphere(r, Normal, nDir, MT);
						InObj = true;
						N_to = 2.0f;

						in_Color += 4.0f * M_PI * RenderScene(Ray(Pos + 0.000001 * nDir, nDir, N_to), depth + 1, scene, camera, MT, InObj, false);
						
						weight = Vec3::Ones() / (russian_p * p);
						MultVec(weight, in_Color);
						MultVec(alpha_s, in_Color);
						MultVec(Vec3(exp(-sigma_a[0] * PSL), exp(-sigma_a[1] * PSL), exp(-sigma_a[2] * PSL)), in_Color);

					}
					else{

						InObj = false;
						N_to = 1.0f;
						p = ISample::NextDirDiffuse(r, Normal, nDir, MT);
						
						in_Color += 4.0f * M_PI * RenderScene(Ray(Pos + 0.000001 * nDir, nDir, N_to), depth + 1, scene, camera, MT, InObj, false);

						weight = Color / (russian_p * p);
						MultVec(weight, in_Color);
						MultVec(Vec3(exp(-sigma_a[0] * FL), exp(-sigma_a[1] * FL), exp(-sigma_a[2] * FL)), in_Color);

					}


				}
				else{

					float RP = MT.genrand64_real1();
					

					
					if (RP < 0.8){
						p = ISample::NextDirDiffuse(r, -Normal, nDir, MT);
						InObj = true;
						N_to = 2.0f;
						in_Color += 1.0f * RenderScene(Ray(Pos + 0.000000001 * nDir, nDir, N_to), depth + 1, scene, camera, MT, InObj, false);
					}
					else{
						
						p = ISample::NextDirDiffuse(r, Normal, nDir, MT);
						InObj = false;
						N_to = 1.0f;
						in_Color += 1.0f * RenderScene(Ray(Pos + 0.000000001 * nDir, nDir, N_to), depth + 1, scene, camera, MT, InObj, false);

					}

					weight = 1.0f * Color / (russian_p * p);
					MultVec(weight, in_Color);

				}


			}
			else if(!InObj){
				///////////////////////
				//拡散面

				d_Color = Vec3::Zero();
				NextEventEstimate(
					std::max(0, 5 - depth),
					scene,
					MT,
					Pos,
					Normal,
					r,
					Scene::DiffuseBRDF,
					d_Color);

				/*if (m_IBL != NULL){
					NextEventIBL(
						std::max(0, 0 - depth),
						scene,
						MT,
						Pos,
						Normal,
						r,
						Scene::DiffuseBRDF,
						ibl_Color);
				}*/
				ibl_Color = Vec3(0, 0, 0);
				
				p = ISample::NextDirDiffuse(r, Normal, nDir, MT);
				NextEvent = true;

				in_Color += RenderScene(Ray(Pos + 0.0001 * nDir, nDir, N_to), depth + 1, scene, camera, MT, InObj, true);

				weight = Color / (russian_p * p);
				MultVec(weight, in_Color);

				d_Color += ibl_Color / (M_PI);
				weight = Color / (russian_p * p);
				MultVec(weight, d_Color);

				in_Color += d_Color;

			}

			///////////////////////////////////////////////////////////
			//return Emit + in_Color;
			if (NE){
				return in_Color;
			}
			else{
				return in_Color + Emit;
			}
			


			
		}
		else{
			if (m_IBL){
				return m_IBL->GetColor2(r.m_Dir);
			}
			else{
				return Vec3(0, 0, 0);
			}
		}

		return Vec3(1, 0, 0);

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
			float NdNL = std::max(0.0f, std::min(1.0f, ((LP - Pos).normalized().dot(-LM._Normal))));

			if (ffa.size() > 0){
				if (ffa[0]._t + 0.000005 > (LP - Pos).norm()){
					float G = NdL * NdNL / pow(ffa[0]._t, 2.0f);
					result += 
						brdf((Pos - LP).normalized(),-ray.m_Dir) * 
						LM._Color * LM._Power *  
						G * 
						scene.GetSumLightArea();
				}
			}
			

		}

		result /= SpL;

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


Pathtracing::Pathtracing(int IMGwidth, int IMGheight) : m_Core(new PTCore), iRender(IMGwidth, IMGheight){

}

Pathtracing::~Pathtracing(){

}


void Pathtracing::RenderScene(const Scene& scene, const Camera& camera, const iIBL* ibl, int sample, int s_sample){
	return;
}



void Pathtracing::RenderSceneT(
	const Scene& scene, 
	const Camera& camera, 
	const iIBL* ibl, 
	int sample, 
	int s_sample,
	int offset){

	m_Core->m_IBL = ibl;


	delete[] m_Core->img;
	m_Core->img = new float[GetIMGWidth()*GetIMGHeight() * 4];
	for (int i = 0; i < GetIMGWidth()*GetIMGHeight() * 4; i++){
		m_Core->img[i] = 0.0f;
	}



	time_t st = time(NULL);
	time_t stt = st;

	Vec3 ssdx = 0.5f * camera._sx / (float)GetIMGWidth();
	Vec3 ssdy = 0.5f * camera._sy / (float)GetIMGHeight();


	int ImgCount = 0;
	const int TITB = 30;
	const int EDTIME = 5 * 60 - 5;
	int NLOOP = 0;
	float AveTime = 0.0f;
	
	
	time_t Ltime;
	Ltime = time(NULL);

	printf("ぱすとれー\n");
	while (true){

		if (NLOOP % 1 == 0){
			std::stringstream st;
			int RA = 6 + 5.001 * (sin(NLOOP / 1.0f));
			int LA = 11 - RA;
			st << "\r";
			for (int aa = 0; aa < RA; aa++){
				st << "*";
			}
			st << " ☆彡 ";
			for (int aa = 0; aa < LA; aa++){
				st << "*";
			}
			fprintf(stderr, st.str().c_str());
		}


		time_t ed = time(NULL);
		int rtm = int(ed) - int(st);



		if (NLOOP > 0){
			time_t Pass = time(NULL) - Ltime;
			Ltime = time(NULL);

			float w = 1.0f / (float)(NLOOP);
			AveTime = AveTime * (1.0f - w) + (float)Pass*(w);
		}



		if (rtm > TITB){

			std::stringstream FileName("");
			FileName << "MMTIMG//IMG_";
			FileName << ImgCount;
			FileName << "_";
			FileName << NLOOP*s_sample*s_sample*sample;
			FileName << "N";
			FileName << ".bmp";
			WriteImage(FileName.str().c_str());
			st = ed;
			ImgCount++;
		}
		if (AveTime + (ed - stt) + offset > EDTIME){
			printf("\nしゅーりょ\n");
			break;
		}


#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
		for (int i = 0; i < GetIMGWidth() * GetIMGHeight(); i++){

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
						Color += m_Core->RenderScene(r, 0, scene, camera, MT, false, false);
					}

				}
			}
			Color /= (float)(s_sample*s_sample*sample);

			for (int ci = 0; ci < 3; ci++){
				float W = (1.0f / float(NLOOP + 1));
				float CL =
					(1.0f - W) * (float)m_Core->img[4 * (w + h*GetIMGWidth()) + ci] +
					W * (float)std::max(Color[ci], 0.0f);
				m_Core->img[4 * (w + h*GetIMGWidth()) + ci] = CL;
			}
			m_Core->img[4 * (w + h*GetIMGWidth()) + 3] = 1.0f;


		}



		NLOOP++;
		

	}

	printf("\n");
	time_t ed = time(NULL);
	int rtm = offset + int(ed) - int(stt);
	int th, tm, ts;
	th = rtm / 3600;
	tm = (rtm - th * 3600) / 60;
	ts = rtm - 3600 * th - 60 * tm;
	printf("経過　%d時間 %d分　%d秒\n", th, tm, ts);


}


void Pathtracing::WriteImage(const char* FileName){

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
			for (int cc = 0; cc < 4; cc++){
				img[4 * (w + (GetIMGHeight() - h - 1)*GetIMGWidth()) + cc] = std::max(0, std::min(255, int(m_Core->img[4 * (w + h*GetIMGWidth()) + cc] * 255)));
			}
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




