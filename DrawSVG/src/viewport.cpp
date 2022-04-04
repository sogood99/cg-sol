#include "viewport.h"

#include "CMU462.h"

namespace CMU462 {

void ViewportImp::set_viewbox(float centerX, float centerY, float vspan) {
    // Task 5 (part 2):
    // Set svg coordinate to normalized device coordinate transformation. Your input
    // arguments are defined as normalized SVG canvas coordinates.
    this->centerX = centerX;
    this->centerY = centerY;
    this->vspan = vspan;

    double inv_doub_vs = 1 / (vspan * 2);
    double m_flat[9] = {inv_doub_vs, 0, -centerX * inv_doub_vs + 0.5, 0, inv_doub_vs, -centerY * inv_doub_vs + 0.5, 0, 0, 1};
    set_svg_2_norm(Matrix3x3(m_flat));
}

void ViewportImp::update_viewbox(float dx, float dy, float scale) {
    this->centerX -= dx;
    this->centerY -= dy;
    this->vspan *= scale;
    set_viewbox(centerX, centerY, vspan);
}

}  // namespace CMU462
