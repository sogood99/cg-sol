#include "software_renderer.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "triangulation.h"

using namespace std;

namespace CMU462 {

// Implements SoftwareRenderer //

void SoftwareRendererImp::draw_svg(SVG& svg) {
    // set top level transformation
    transformation = svg_2_screen;

    // draw all elements
    for (size_t i = 0; i < svg.elements.size(); ++i) {
        draw_element(svg.elements[i]);
    }

    // draw canvas outline
    Vector2D a = transform(Vector2D(0, 0));
    a.x--;
    a.y--;
    Vector2D b = transform(Vector2D(svg.width, 0));
    b.x++;
    b.y--;
    Vector2D c = transform(Vector2D(0, svg.height));
    c.x--;
    c.y++;
    Vector2D d = transform(Vector2D(svg.width, svg.height));
    d.x++;
    d.y++;

    rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
    rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
    rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
    rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

    // resolve and send to render target
    resolve();
}

void SoftwareRendererImp::set_sample_rate(size_t sample_rate) {
    // Task 4:
    // You may want to modify this for supersampling support
    this->sample_rate = sample_rate;
    init_supersample();
}

void SoftwareRendererImp::set_render_target(unsigned char* render_target,
                                            size_t width, size_t height) {
    // Task 4:
    // You may want to modify this for supersampling support
    this->render_target = render_target;
    this->target_w = width;
    this->target_h = height;

    init_supersample();
}

void SoftwareRendererImp::init_supersample() {
    // Initialize the supersampling buffer

    if (this->supersample_target != nullptr) {
        delete[] this->supersample_target;
    }
    this->supersample_target = new unsigned char[4 * this->target_w * this->target_h * sample_rate * sample_rate];
    memset(this->supersample_target, 255, 4 * this->target_w * this->target_h * sample_rate * sample_rate);
}

void SoftwareRendererImp::draw_element(SVGElement* element) {
    // Task 5 (part 1):
    // Modify this to implement the transformation stack

    // save the transformation for stack
    Matrix3x3 past_transformation = transformation;
    transformation = transformation * element->transform;

    switch (element->type) {
        case POINT:
            draw_point(static_cast<Point&>(*element));
            break;
        case LINE:
            draw_line(static_cast<Line&>(*element));
            break;
        case POLYLINE:
            draw_polyline(static_cast<Polyline&>(*element));
            break;
        case RECT:
            draw_rect(static_cast<Rect&>(*element));
            break;
        case POLYGON:
            draw_polygon(static_cast<Polygon&>(*element));
            break;
        case ELLIPSE:
            draw_ellipse(static_cast<Ellipse&>(*element));
            break;
        case IMAGE:
            draw_image(static_cast<Image&>(*element));
            break;
        case GROUP:
            draw_group(static_cast<Group&>(*element));
            break;
        default:
            break;
    }

    // restore the stack
    transformation = past_transformation;
}

// Primitive Drawing //

void SoftwareRendererImp::draw_point(Point& point) {
    Vector2D p = transform(point.position);
    rasterize_point(p.x, p.y, point.style.fillColor);
}

void SoftwareRendererImp::draw_line(Line& line) {
    Vector2D p0 = transform(line.from);
    Vector2D p1 = transform(line.to);
    rasterize_line(p0.x, p0.y, p1.x, p1.y, line.style.strokeColor);
}

void SoftwareRendererImp::draw_polyline(Polyline& polyline) {
    Color c = polyline.style.strokeColor;

    if (c.a != 0) {
        int nPoints = polyline.points.size();
        for (int i = 0; i < nPoints - 1; i++) {
            Vector2D p0 = transform(polyline.points[(i + 0) % nPoints]);
            Vector2D p1 = transform(polyline.points[(i + 1) % nPoints]);
            rasterize_line(p0.x, p0.y, p1.x, p1.y, c);
        }
    }
}

void SoftwareRendererImp::draw_rect(Rect& rect) {
    Color c;

    // draw as two triangles
    float x = rect.position.x;
    float y = rect.position.y;
    float w = rect.dimension.x;
    float h = rect.dimension.y;

    Vector2D p0 = transform(Vector2D(x, y));
    Vector2D p1 = transform(Vector2D(x + w, y));
    Vector2D p2 = transform(Vector2D(x, y + h));
    Vector2D p3 = transform(Vector2D(x + w, y + h));

    // draw fill
    c = rect.style.fillColor;
    if (c.a != 0) {
        rasterize_triangle(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c);
        rasterize_triangle(p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c);
    }

    // draw outline
    c = rect.style.strokeColor;
    if (c.a != 0) {
        rasterize_line(p0.x, p0.y, p1.x, p1.y, c);
        rasterize_line(p1.x, p1.y, p3.x, p3.y, c);
        rasterize_line(p3.x, p3.y, p2.x, p2.y, c);
        rasterize_line(p2.x, p2.y, p0.x, p0.y, c);
    }
}

void SoftwareRendererImp::draw_polygon(Polygon& polygon) {
    Color c;

    // draw fill
    c = polygon.style.fillColor;
    if (c.a != 0) {
        // triangulate
        vector<Vector2D> triangles;
        triangulate(polygon, triangles);

        // draw as triangles
        for (size_t i = 0; i < triangles.size(); i += 3) {
            Vector2D p0 = transform(triangles[i + 0]);
            Vector2D p1 = transform(triangles[i + 1]);
            Vector2D p2 = transform(triangles[i + 2]);
            rasterize_triangle(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c);
        }
    }

    // draw outline
    c = polygon.style.strokeColor;
    if (c.a != 0) {
        int nPoints = polygon.points.size();
        for (int i = 0; i < nPoints; i++) {
            Vector2D p0 = transform(polygon.points[(i + 0) % nPoints]);
            Vector2D p1 = transform(polygon.points[(i + 1) % nPoints]);
            rasterize_line(p0.x, p0.y, p1.x, p1.y, c);
        }
    }
}

void SoftwareRendererImp::draw_ellipse(Ellipse& ellipse) {
    // Extra credit
}

void SoftwareRendererImp::draw_image(Image& image) {
    Vector2D p0 = transform(image.position);
    Vector2D p1 = transform(image.position + image.dimension);

    rasterize_image(p0.x, p0.y, p1.x, p1.y, image.tex);
}

void SoftwareRendererImp::draw_group(Group& group) {
    for (size_t i = 0; i < group.elements.size(); ++i) {
        draw_element(group.elements[i]);
    }
}

// Rasterization //

// The input arguments in the rasterization functions
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= target_w) return;
    if (sy < 0 || sy >= target_h) return;

    // fill sample - NOT doing alpha blending!
    for (int k = 0; k < sample_rate; k++) {
        supersample_target[4 * (sx * sample_rate + sy * target_w * sample_rate + k) + 0] = color.r * 255;
        supersample_target[4 * (sx * sample_rate + sy * target_w * sample_rate + k) + 1] = color.g * 255;
        supersample_target[4 * (sx * sample_rate + sy * target_w * sample_rate + k) + 2] = color.b * 255;
        supersample_target[4 * (sx * sample_rate + sy * target_w * sample_rate + k) + 3] = color.a * 255;
    }
}

inline void swap(float* a, float* b) {
    float temp;
    temp = *a;
    *a = *b;
    *b = temp;
}

void SoftwareRendererImp::rasterize_line(float x0, float y0, float x1, float y1,
                                         Color color) {
    // Task 2:
    // Implement line rasterization
    // Bresenham algorithm

    if (x1 < x0) {
        swap(&x0, &x1);
        swap(&y0, &y1);
    }

    int eps = 0, dx = x1 - x0, dy = y1 - y0, x = x0, y = y0;

    if (dy >= 0) {
        if (dy <= dx) {
            for (int x = x0; x <= x1; x++) {
                rasterize_point(x, y, color);
                if (2 * (eps + dy) < dx) {
                    eps += dy;
                } else {
                    y++;
                    eps += dy - dx;
                }
            }
        } else {
            for (int y = y0; y <= y1; y++) {
                rasterize_point(x, y, color);
                if (2 * (eps + dx) < dy) {
                    eps += dx;
                } else {
                    x++;
                    eps += dx - dy;
                }
            }
        }
    } else {
        if (abs(dy) <= dx) {
            for (int x = x0; x <= x1; x++) {
                rasterize_point(x, y, color);
                if (abs(2 * (eps + dy)) < dx) {
                    eps += dy;
                } else {
                    y--;
                    eps += dy + dx;
                }
            }
        } else {
            // y1 < y0 so swap y0 and y1
            swap(&x0, &x1);
            swap(&y0, &y1);
            x = x0;
            dx = -dx;
            dy = -dy;

            for (int y = y0; y <= y1; y++) {
                rasterize_point(x, y, color);
                if (abs(2 * (eps + dx)) < dy) {
                    eps += dx;
                } else {
                    x--;
                    eps += dx + dy;
                }
            }
        }
    }
}

inline bool test_half_plane(float x, float y, float x0, float y0,
                            float x1, float y1, float x2, float y2) {
    // Test if point (x,y) is in same half plane as (x2, y2)
    // splitted by line through (x0,y0) -> (x1, y1)
    float a = y1 - y0, b = x0 - x1, c = y0 * (x1 - x0) - x0 * (y1 - y0);
    float eval_pt = a * x + b * y + c, eval_test = a * x2 + b * y2 + c;
    if ((eval_pt >= 0 && eval_test >= 0) || (eval_pt < 0 && eval_test < 0)) {
        return true;
    }
    return false;
}

void SoftwareRendererImp::rasterize_triangle(float x0, float y0, float x1,
                                             float y1, float x2, float y2,
                                             Color color) {
    // Task 3:
    // Implement triangle rasterization
    // Naive traingle rasterizer

    int min_x = (int)floor(min(x0, min(x1, x2))), max_x = (int)ceil(max(x0, max(x1, x2)));
    int min_y = (int)floor(min(y0, min(y1, y2))), max_y = (int)ceil(max(y0, max(y1, y2)));

    for (int i = max(min_y, 0); i < min(max_y, (int)target_h); i++) {
        for (int j = max(min_x, 0); j < min(max_x, (int)target_w); j++) {
            for (int k = 0; k < sample_rate; k++) {
                float x = j + (float)(2 * k + 1) / (2 * sample_rate), y = i + (float)(2 * k + 1) / (2 * sample_rate);
                if (test_half_plane(x, y, x0, y0, x1, y1, x2, y2) &&
                    test_half_plane(x, y, x1, y1, x2, y2, x0, y0) &&
                    test_half_plane(x, y, x0, y0, x2, y2, x1, y1)) {
                    supersample_target[4 * (j * sample_rate + i * target_w * sample_rate + k) + 0] = color.r * 255;
                    supersample_target[4 * (j * sample_rate + i * target_w * sample_rate + k) + 1] = color.g * 255;
                    supersample_target[4 * (j * sample_rate + i * target_w * sample_rate + k) + 2] = color.b * 255;
                    supersample_target[4 * (j * sample_rate + i * target_w * sample_rate + k) + 3] = color.a * 255;
                }
            }
        }
    }
}

void SoftwareRendererImp::rasterize_image(float x0, float y0, float x1,
                                          float y1, Texture& tex) {
    // Task 6:
    // Implement image rasterization
    int min_x = (int)floor(min(x0, x1)), max_x = (int)ceil(max(x0, x1));
    int min_y = (int)floor(min(y0, y1)), max_y = (int)ceil(max(y0, y1));

    for (int i = max(min_y, 0); i < min(max_y, (int)target_h); i++) {
        for (int j = max(min_x, 0); j < min(max_x, (int)target_w); j++) {
            for (int k = 0; k < sample_rate; k++) {
                float x = j + (float)(2 * k + 1) / (2 * sample_rate), y = i + (float)(2 * k + 1) / (2 * sample_rate);

                float u = (x - x0) / (x1 - x0), v = (y - y0) / (y1 - y0);
                // Color color = sampler->sample_bilinear(tex, u, v, 0);
                Color color = sampler->sample_trilinear(tex, u, v, 250 / (x1 - x0), 250 / (y1 - y0));

                supersample_target[4 * (j * sample_rate + i * target_w * sample_rate + k) + 0] = color.r;
                supersample_target[4 * (j * sample_rate + i * target_w * sample_rate + k) + 1] = color.g;
                supersample_target[4 * (j * sample_rate + i * target_w * sample_rate + k) + 2] = color.b;
                supersample_target[4 * (j * sample_rate + i * target_w * sample_rate + k) + 3] = color.a;
            }
        }
    }
}

// resolve samples to render target
void SoftwareRendererImp::resolve(void) {
    // Task 4:
    // Implement supersampling
    // You may also need to modify other functions marked with "Task 4".

    for (int i = 0; i < target_h; i++) {
        for (int j = 0; j < target_w; j++) {
            render_target[4 * (j + i * target_w) + 0] = 0;  // r
            render_target[4 * (j + i * target_w) + 1] = 0;  // g
            render_target[4 * (j + i * target_w) + 2] = 0;  // b
            render_target[4 * (j + i * target_w) + 3] = 0;  // a

            for (int k = 0; k < sample_rate; k++) {
                render_target[4 * (j + i * target_w) + 0] += supersample_target[4 * (j * sample_rate + i * target_w * sample_rate + k) + 0] / sample_rate;
                render_target[4 * (j + i * target_w) + 1] += supersample_target[4 * (j * sample_rate + i * target_w * sample_rate + k) + 1] / sample_rate;
                render_target[4 * (j + i * target_w) + 2] += supersample_target[4 * (j * sample_rate + i * target_w * sample_rate + k) + 2] / sample_rate;
                render_target[4 * (j + i * target_w) + 3] += supersample_target[4 * (j * sample_rate + i * target_w * sample_rate + k) + 3] / sample_rate;
            }
        }
    }

    memset(this->supersample_target, 255, 4 * this->target_w * this->target_h * sample_rate * sample_rate);
    return;
}

}  // namespace CMU462
