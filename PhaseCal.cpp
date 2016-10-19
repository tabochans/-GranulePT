#include"PhaseCal.h"
#include"Grain.h"

using namespace PTUtility;


float PhaseEstimater::GetFromCos(float c){
	float ac = acos(c);
	int ss = (ac / M_PI)*_D;
	int ID = std::max(0, std::min(int(_D - 1), ss));
	float a;
	if (ID > _D - 10){
		return m_data[_D -  10];
	}
	else if (ID < 10){
		return m_data[10];
	}
	a = (ac - M_PI * (ID / (float)_D));
	a /= (M_PI / (float)_D);
	a = std::max(0.0f, std::min(1.0f, a));
	
	float re = (1 - a) * m_data[ID] + a * m_data[ID + 1];
	return re;
}

void PhaseEstimater::EstimatePhaseFunction(iGrain* grain, float& Lamda_Dt, float& Lamda_S, float& Beta){

	RandomMT mm(88);
	int Dt = 0;
	int S = 0;
	int S2 = 0;
	std::vector<float> Fg(_D, 0.0f);
	int Nsample = 100;

	float LDT = 0.0f;
	float LS = 0.0f;

	int _DD = _D * 3;

	for (int i = 0; i < _DD; i++){
		for (int j = 0; j < _DD; j++){

			float ru = 2.0f * grain->GetR() * (i - _DD * 0.5f) / float(_DD);
			float rv = 2.0f * grain->GetR() * (j - _DD * 0.5f) / float(_DD);

			if (ru*ru + rv*rv > grain->GetR() * grain->GetR()){
				continue;
			}
			
			
			for (int s = 0; s < Nsample; s++){
				Vec3 InPos(ru, rv, -grain->GetR() * 10.0f);
				Vec3 InDir(0, 0, 1);
				Vec3 EnterPos(0, 0, 0);
				float pdf = 1.0f;
				Ray result(Vec3::Zero(), Vec3(0, 1, 0));

				int ct = 0;
				S2++;
				while (ct < 100 && grain->NextRay(Ray(InPos, InDir), Vec3::Zero(), 1, result, pdf, mm)){
					
					if (ct > 0){
						
					}
					else{
						EnterPos = result.m_Org;
					}

					ct++;
					InPos = result.m_Org;
					InDir = result.m_Dir.normalized();


				}
				if (ct == 0){
					//交差とれなかったとき
					LDT += 2.0f * sqrt(grain->GetR()*grain->GetR() * (1.0001) - ru*ru - rv*rv);
					continue;
				}
				
				//交差した
				LS += (InPos - EnterPos).norm();
				
				S++;
				
				float CC = InDir.dot(Vec3(0, 0, 1));
				int ID = _D * acos(CC) / (M_PI + 0.00001);
				float SS = sqrt(1.0 + FLT_MIN - CC*CC);

				Fg[ID] += 1.0f / SS;

			}
			
		}
	}

	if (S2 - S > 0){
		LDT /= (S2 - S);
	}
	LDT /= grain->GetR();

	LS /= (S);
	LS /= grain->GetR();

	m_data.resize(_D);

	Beta = S / (float)S2;
	float total = 0.0f;
	for (int i = 0; i < _D; i++){
		total += (Fg[std::max(10, std::min(_D - 10, i))] * sin(M_PI * i / (float)(_D)));
	}
	total *= ((M_PI / (float)(_D)));

	for (int i = 0; i < _D; i++){
		m_data[i] = 0.9999 * Fg[std::max(10, std::min(_D - 10, i))] / (total * 2.0f * M_PI);
	}


	Lamda_S = LS;
	Lamda_Dt = LDT;

	m_Beta = Beta;
	m_Dt = LDT;
	m_Ds = LS;

	grain->m_Beta = Beta;
	grain->m_Lamda_dt = LDT;
	grain->m_Lamda_s = Vec3::Ones() * LS;


}


void PhaseEstimater::EstimatePhaseFunction2(iGrain* grain, float& Lamda_Dt, float& Lamda_S, float& Beta){

	RandomMT mm(88);
	int Dt = 0;
	int S = 0;
	int S2 = 0;
	std::vector<float> Fg(_D, 0.0f);
	int Nsample = 1;

	float LDT = 0.0f;
	float LS = 0.0f;

	int _DD = 1 * _D;

	std::cout << "つぶつぶ計算。。。" << std::endl;
	for (int i = 0; i < _DD; i++){

		fprintf(stderr, "\r %5.2f%%", 100.0*i / _DD);

		for (int j = 0; j < _DD; j++){

			float ru = 2.0f * grain->GetR() * (i - _DD * 0.5f) / float(_DD);
			float rv = 2.0f * grain->GetR() * (j - _DD * 0.5f) / float(_DD);

			if (ru*ru + rv*rv > grain->GetR() * grain->GetR()){
				continue;
			}


			for (int t = 0; t < _DD; t++){
				for (int p = 0; p < _DD; p++){


					bool isHit = true;
					for (int s = 0; s < Nsample; s++){


						Vec3 InPos(ru, rv, -grain->GetR() * 10.0f);
						Vec3 InDir(0, 0, 1);
						Eigen::Vector3f tt;

						float Rot = 2.0f * M_PI * t / (float)_DD;
						Eigen::Vector3f jI(cosf(M_PI * p / (float)_DD), sinf(M_PI * p / (float)_DD), 0);
						Eigen::Quaternionf q(Eigen::AngleAxisf(Rot, jI));
						 
						tt = q * Eigen::Vector3f(InPos[0], InPos[1], InPos[2]);
						InPos = Vec3(tt[0], tt[1], tt[2]);

						tt = q * Eigen::Vector3f(InDir[0], InDir[1], InDir[2]);
						InDir = Vec3(tt[0], tt[1], tt[2]);

						InDir.normalize();
						Vec3 DInDir = InDir;



						float pdf = 1.0f;
						Ray result(Vec3::Zero(), Vec3(0, 1, 0));

						int ct = 0;
						while (ct < 100 && grain->NextRay(Ray(InPos, InDir), Vec3::Zero(), 1, result, pdf, mm)){

							if (ct > 0){
								LS += (InPos - result.m_Org).norm();
							}

							ct++;
							InPos = result.m_Org;
							InDir = result.m_Dir.normalized();


						}
						if (ct == 0){
							//交差とれなかったとき
							LDT += 2.0f * sqrt(grain->GetR()*grain->GetR() + FLT_MIN - ru*ru - rv*rv);
							s = Nsample;
							isHit = false;
							continue;
						}


						float CC = InDir.dot(DInDir)*0.9999;
						int ID = _D * acos(CC) / (M_PI + 0.00001);
						float SS = sqrt(1.0 + FLT_MIN - CC*CC);

						Fg[ID] += 1.0f / SS;

					}
					if (isHit){
						S2++;
						S++;
					}
					else{
						S2++;
					}



				}
			}
		}
	}

	if (S2 - S > 0){
		LDT /= (S2 - S);
	}
	LDT /= grain->GetR();

	LS /= (S);
	LS /= grain->GetR();

	m_data.resize(_D);

	Beta = S / (float)S2;
	float total = 0.0f;
	for (int i = 0; i < _D; i++){
		total += (Fg[std::max(10, std::min(_D - 10, i))] * sin(M_PI * i / (float)(_D)));
	}
	total *= ((M_PI / (float)(_D)));

	for (int i = 0; i < _D; i++){
		m_data[i] = 0.9999 * Fg[std::max(10, std::min(_D - 10, i))] / (total * 2.0f * M_PI);
	}


	Lamda_S = LS;
	Lamda_Dt = LDT;

	m_Beta = Beta;
	m_Dt = LDT;
	m_Ds = LS;

	std::cout << std::endl << "つぶつぶ計算。。。おわり" << std::endl;

	grain->m_Beta = Beta;
	grain->m_Lamda_dt = LDT;
	grain->m_Lamda_s = Vec3::Ones() * LS;

}

void PhaseEstimater::WritePhaseFunction(std::ofstream& file){

	for (int i = 0; i < _D; i++){
		int si = i;
		file << "P " << ", ";
		file << M_PI * i / (float)_D;
		file << ", ";
		file << m_data[si] << std::endl;
	}
	file << "Dt , " << m_Dt << std::endl;
	file << "Ds , " << m_Ds << std::endl;
	file << "Beta , " << m_Beta << std::endl;
}


void PhaseEstimater::LoadPhaseFunction(iGrain* grain, std::ifstream& file){
	
	m_data.clear();
	while (!file.eof()){
		std::string Line;
		float th, p;
		float ff;
		std::string CM("");
		std::string ID("");

		std::getline(file, Line);

		std::istringstream LS(Line);
		LS >> ID;
		if (ID == std::string("P")){
			LS >> CM;
			LS >> th;
			LS >> CM;
			LS >> p;
			m_data.push_back(p);
		}
		else if (ID == std::string("Dt")){
			LS >> CM;
			LS >> ff;
			m_Dt = ff;
		}
		else if (ID == std::string("Ds")){
			LS >> CM;
			LS >> ff;
			m_Ds = ff;
		}
		else if (ID == std::string("Beta")){
			LS >> CM;
			LS >> ff;
			m_Beta = ff;
		}

	}

	_D = m_data.size();

	grain->m_Beta = m_Beta;
	grain->m_Lamda_dt = m_Dt;
	grain->m_Lamda_s = Vec3::Ones() * m_Ds;
}









