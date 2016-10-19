#pragma once

#include"Render.h"
#include"stb_image_write.h"
#include"MT.h"
#include"ImportanceSample.h"
#include"Aggrigate.h"
#include"Vec3.h"
#include"Vec2.h"
#include<sstream>
#include<memory>
#include<algorithm>
#include<time.h>

namespace PTUtility{

	class GPathTrace : public iRender{

	private:

		class GPTCore;
		std::shared_ptr<GPTCore> m_Core;


	protected:


	public:


		virtual void RenderScene(
			const Scene& scene,
			const Camera& camera,
			const iIBL* ibl,
			int sample,
			int s_sample);

		virtual void WriteImage(const char* FileName);

		GPathTrace(int IMGwidth, int IMGheight);
		virtual ~GPathTrace();

	};

}