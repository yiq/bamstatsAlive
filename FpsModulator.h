#ifndef YICPPLIB_FPSMODULATOR_H_
#define YICPPLIB_FPSMODULATOR_H_

#include <chrono>
#include <memory>

namespace YiCppLib {

	/** Frame-per-Second Modulation
	 *  
	 *  Often there is need to adjust one value so that some other periodical
	 *  updates can happen at a regular interval. A most straightforward case
	 *  would be redraw rate based on rendering complexity, but this concept
	 *  is equally applicable to stats update based on records processed in an
	 *  iobio stats-alive server.
	 *
	 *  Users of this class would instanciate an object with the target
	 *  variable to modulate, and call "redraw" everytime an update is
	 *  produced. The object will then directly affect the variable so that a
	 *  target FPS can be approximated
	 */

	template<class T>
	class FpsModulator {
		public:
			/** The constructor
			 * @param modulated The target variable to modulate
			 * @param targetFps The desired FPS in milliseconds
			 * @param tolerance The acceptable discrepencies between actual and desired FPS in milliseconds
			 */
			FpsModulator(T& modulated, size_t targetFps, size_t tolerance) 
				: m_modulated(modulated), m_targetFps(targetFps), m_tolerance(tolerance), m_lastRedraw(std::chrono::steady_clock::now())
			{}

			/** Notify the modulator that a redraw just occurred
			 *
			 *  It is the caller's responsibility to notify the modulator that a
			 *  redraw just occurred, so that it can then calculate the time
			 *  elapsed, and modulate target variable accordingly
			 */
			void redraw() {
				auto now = std::chrono::steady_clock::now();
				auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(now - m_lastRedraw).count();
				m_lastRedraw = std::chrono::steady_clock::now();

				auto lowerAcceptance = m_targetFps < m_tolerance ? 0 : m_targetFps - m_tolerance;
				auto upperAcceptance = m_targetFps + m_tolerance;

				if(dur < lowerAcceptance) {
					// too fast, increase target
					m_modulated *= 1.5;
				}
				else if (dur > upperAcceptance) {
					// too slow, decrease
					m_modulated /= 1.5;
				}
			}

		private:
			T& m_modulated;
			size_t m_targetFps;
			size_t m_tolerance;
			std::chrono::steady_clock::time_point m_lastRedraw;
	};
}

#endif
