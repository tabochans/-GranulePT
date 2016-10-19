#pragma once
#include"Camera.h"
#include"Scene.h"
#include"iIBL.h"

namespace PTUtility{

	class iRender{

	private:

		int m_Width;
		int m_Height;

	protected:


	public:

		int GetIMGWidth(){ return m_Width; }
		int GetIMGHeight(){ return m_Height; }

		iRender(int IMGwidth, int IMGheight) : m_Width(IMGwidth > 0 ? IMGwidth : 10), m_Height(IMGheight > 0 ? IMGheight : 10){ ; }
		virtual ~iRender(){ ; }

		virtual void RenderScene(
			const Scene& scene, 
			const Camera& camera, 
			const iIBL* ibl,
			int sample, 
			int s_sample) = 0;

		virtual void WriteImage(const char* FileName) = 0;

	};

}

