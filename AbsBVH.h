#pragma once
#include"Vec3.h"
#include"Facet.h"
#include"Ray.h"
#include"Intersect.h"
#include<vector>

class AbsG{

private:


protected:


public:

	//れいとれデータの受け渡し用の構造体なの
	struct RF{
		const PTUtility::Facet* _pF;
		float _u, _v, _t;
		RF(const PTUtility::Facet* fp, float u, float v, float t) : _pF(fp), _u(u), _v(v), _t(t){
		}
		bool operator<(const RF& ref){
			return this->_t < ref._t;
		}
		bool operator>(const RF& ref){
			return this->_t > ref._t;
		}
		bool operator<=(const RF& ref){
			return this->_t <= ref._t;
		}
		bool operator>=(const RF& ref){
			return this->_t >= ref._t;
		}
		RF& operator=(const RF& ref){
			_pF = ref._pF;
			_u = ref._u;
			_v = ref._v;
			_t = ref._t;
			return *this;
		}
	};

	virtual void CreateGraph(const std::vector<PTUtility::Facet>& FacetArray) = 0;
	virtual void ClearGraph() = 0;
	virtual void PrintGraph(std::ostream& ost) const{ ; }
	virtual void WriteObj(std::ostream& ost) const{ ; }
	virtual const std::string Print() const{
		return std::string("iGraph");
	}

	//Rayと交差した三角形をすべて返してくれるの
	virtual int GetElements(const PTUtility::Ray& ray, std::vector<RF>& Facets) const = 0;
	
	//球と交差した三角形をすべて返してくれるの
	virtual int GetElements(const PTUtility::Vec3& Center, float R, std::vector<RF>& Facets){ return 0; }


	virtual float GetSceneWidth(){ return 1.0f; }
	virtual PTUtility::Vec3 GetSceneCenter(){ return PTUtility::Vec3::Zero(); }

	AbsG(){

	}
	virtual ~AbsG(){

	}

};