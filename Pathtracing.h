#pragma once
#include"Render.h"
#include"stb_image_write.h"
#include"MT.h"
#include"ImportanceSample.h"
#include<sstream>
#include<memory>
#include<algorithm>
#include<time.h>

namespace PTUtility{

	class Pathtracing : public iRender{

		class PTCore;

	private:

		std::shared_ptr<PTCore> m_Core;

	protected:


	public:

		Pathtracing(int IMGwidth, int IMGheight);
		virtual ~Pathtracing();

		virtual void RenderScene(
			const Scene& scene,
			const Camera& camera,
			const iIBL* ibl,
			int sample,
			int s_sample);

		virtual void WriteImage(const char* FileName);


		void RenderSceneT(
			const Scene& scene,
			const Camera& camera,
			const iIBL* ibl,
			int sample,
			int s_sample,
			int offset);


	};


}