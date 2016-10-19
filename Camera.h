#pragma once
#include"Vec3.h"


namespace PTUtility{

	struct Camera{

		//�J�����ʒu
		Vec3 _CameraPos;

		//�J��������
		Vec3 _CameraDir;

		//�����
		Vec3 _CameraUp;

		//���[���h���W�ł̃X�N���[���̕�
		float _ScreenWidth;

		//���[���h���W�ł̃X�N���[���̍���
		float _ScreenHeight;

		//�ˉe�ʂ܂ł̋���
		float _ScreenDist;

		//�X�N���[���𒣂�x�N�g��
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

			//�J�������猩�ĉE����
			_sx = (_CameraDir.cross(_CameraUp).normalized()) * _ScreenWidth;

			//�J�������猩�ĉ�����
			_sy = (_CameraDir.cross(_sx).normalized())*_ScreenHeight;

			_Center = _CameraPos + _ScreenDist * _CameraDir;

		}

		~Camera(){

		}


	};

}