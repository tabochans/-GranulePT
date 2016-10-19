#pragma once
#include"Vec.h"
#include<cmath>


namespace PTUtility{

	//åpè≥Ç∑ÇÈÇ»ÇÊÅAê‚ëŒÇæÇºÅHÅI
	class Vec2{

	private:

		float a[2];

	public:

		static Vec2 Ones(){
			return Vec2(1, 1);
		}
		static Vec2 Zero(){
			return Vec2(0, 0);
		}

		static Vec2 MultVec(const Vec2& left, const Vec2& right){
			return Vec2(left.x() * right.x(), left.y() * right.y());
		}
		static Vec2 DivVec(const Vec2& left, const Vec2& right){
			return Vec2(left.x() / right.x(), left.y() / right.y());
		}

		float dot(const Vec2& ref)const{
			return a[0] * ref.x() + a[1] * ref.y();
		}

		float norm()const{
			return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
		}
		float norm2()const{
			return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
		}

		void normalize(){
			float d = norm();
			a[0] /= d;
			a[1] /= d;
			a[2] /= d;
		}
		Vec2 normalized()const{
			float d = norm();
			return Vec2(a[0] / d, a[1] / d);
		}

		inline const Vec2 operator+(const Vec2& ref)const{
			return Vec2(a[0] + ref.x(), a[1] + ref.y());
		}
		inline const Vec2 operator-(const Vec2& ref)const{
			return Vec2(a[0] - ref.x(), a[1] - ref.y());
		}
		inline const Vec2 operator-()const{
			return Vec2(-a[0], -a[1]);
		}
		inline const Vec2 operator+()const{
			return *this;
		}



		inline float x()const{
			return a[0];
		}
		inline float y()const{
			return a[1];
		}
		inline float& x(){
			return a[0];
		}
		inline float& y(){
			return a[1];
		}

		inline float operator[](size_t id)const{
			return a[id];
		}
		inline float& operator[](size_t id){
			return a[id];
		}

		inline float operator()(size_t id)const{
			return a[id];
		}
		inline float& operator()(size_t id){
			return a[id];
		}


		Vec2(const Vec2& ref){
			a[0] = ref.x();
			a[1] = ref.y();
		}
		Vec2& operator=(const Vec2& ref){
			a[0] = ref.x();
			a[1] = ref.y();
			return *this;
		}

		Vec2(){
			a[0] = a[1] = 0.0f;
		}
		Vec2(float x, float y){
			a[0] = x;
			a[1] = y;
		}


		inline Vec2& operator+=(const Vec2& ref){
			a[0] += ref.x();
			a[1] += ref.y();
			return *this;
		}
		inline Vec2& operator-=(const Vec2& ref){
			a[0] -= ref.x();
			a[1] -= ref.y();
			return *this;
		}


		~Vec2(){

		}

	};


	inline const Vec2 operator*(float left, const Vec2& right){
		return Vec2(left*right.x(), left*right.y());
	}
	inline const Vec2 operator*(const Vec2& left, float right){
		return Vec2(right*left.x(), right*left.y());
	}

	inline const Vec2 operator*=(Vec2& left, float right){
		left.x() *= right;
		left.y() *= right;
		return left;
	}
	inline const Vec2 operator/=(Vec2& left, float right){
		left.x() /= right;
		left.y() /= right;
		return left;
	}

	inline const Vec2 operator/(const Vec2& left, float right){
		return Vec2(left.x() / right, left.y() / right);
	}


}