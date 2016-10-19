#pragma once
#include"Vec3.h"


namespace PTUtility{

	struct Camera{

		//カメラ位置
		Vec3 _CameraPos;

		//カメラ方向
		Vec3 _CameraDir;

		//上方向
		Vec3 _CameraUp;

		//ワールド座標でのスクリーンの幅
		float _ScreenWidth;

		//ワールド座標でのスクリーンの高さ
		float _ScreenHeight;

		//射影面までの距離
		float _ScreenDist;

		//スクリーンを張るベクトル
		Vec3 _sx;
		Vec3 _sy;

		Vec3 _Center;

		Camera(
			const Vec3& pos, 
			const Vec3& dir,
			const Vec3& up,
			float width,
			float height,
			float A) : 
			_CameraPos(pos),_CameraDir(dir.normalized()),_CameraUp(up),
			_ScreenWidth(3.0f * width / height),_ScreenHeight(3.0f), _ScreenDist(3.0f * A){

			//カメラから見て右向き
			_sx = (_CameraDir.cross(_CameraUp).normalized()) * _ScreenWidth;

			//カメラから見て下向き
			_sy = (_CameraDir.cross(_sx).normalized())*_ScreenHeight;

			_Center = _CameraPos + _ScreenDist * _CameraDir;

		}

		~Camera(){

		}


	};

}