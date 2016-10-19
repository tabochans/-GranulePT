#pragma once

#include"Ray.h"
#include"Eigen\Core"
#include"Eigen\Geometry"
#include<sstream>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<vector>

namespace PTUtility{

	class iGrain;

	class PhaseEstimater{


	private:
		int _D;
		std::vector<float> m_data;
		float m_Dt;
		float m_Ds;
		float m_Beta;

	public:

		float GetFromCos(float c);

		void EstimatePhaseFunction(iGrain* grain, float& Lamda_Dt, float& Lamda_S, float& Beta);
		void EstimatePhaseFunction2(iGrain* grain, float& Lamda_Dt, float& Lamda_S, float& Beta);

		void LoadPhaseFunction(iGrain* grain, std::ifstream& file);
		void WritePhaseFunction(std::ofstream& file);

		explicit PhaseEstimater(int D) : _D(std::max(50,D)){

		}
		virtual ~PhaseEstimater(){

		}

	};


}