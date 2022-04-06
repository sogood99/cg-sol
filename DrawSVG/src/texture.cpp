#include "texture.h"

#include <assert.h>

#include <algorithm>
#include <iostream>

#include "color.h"

using namespace std;

namespace CMU462 {

inline void uint8_to_float(float dst[4], unsigned char* src) {
    uint8_t* src_uint8 = (uint8_t*)src;
    dst[0] = src_uint8[0] / 255.f;
    dst[1] = src_uint8[1] / 255.f;
    dst[2] = src_uint8[2] / 255.f;
    dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8(unsigned char* dst, float src[4]) {
    uint8_t* dst_uint8 = (uint8_t*)dst;
    dst_uint8[0] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[0])));
    dst_uint8[1] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[1])));
    dst_uint8[2] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[2])));
    dst_uint8[3] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[3])));
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {
    // NOTE:
    // This starter code allocates the mip levels and generates a level
    // map by filling each level with placeholder data in the form of a
    // color that differs from its neighbours'. You should instead fill
    // with the correct data!

    // Task 7: Implement this

    // check start level
    if (startLevel >= tex.mipmap.size()) {
        std::cerr << "Invalid start level";
    }

    // allocate sublevels
    int baseWidth = tex.mipmap[startLevel].width;
    int baseHeight = tex.mipmap[startLevel].height;
    int numSubLevels = (int)(log2f((float)max(baseWidth, baseHeight)));

    numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
    tex.mipmap.resize(startLevel + numSubLevels + 1);

    int width = baseWidth;
    int height = baseHeight;
    for (int i = 1; i <= numSubLevels; i++) {
        MipLevel& level = tex.mipmap[startLevel + i];

        // handle odd size texture by rounding down
        width = max(1, width / 2);
        assert(width > 0);
        height = max(1, height / 2);
        assert(height > 0);

        level.width = width;
        level.height = height;
        level.texels = vector<unsigned char>(4 * width * height);
    }

    // fill all 0 sub levels with interchanging colors (JUST AS A PLACEHOLDER)
    // Color colors[3] = {Color(1, 0, 0, 1), Color(0, 1, 0, 1), Color(0, 0, 1, 1)};
    // for (size_t i = 1; i < tex.mipmap.size(); ++i) {
    //     Color c = colors[i % 3];
    //     MipLevel& mip = tex.mipmap[i];

    //     for (size_t i = 0; i < 4 * mip.width * mip.height; i += 4) {
    //         float_to_uint8(&mip.texels[i], &c.r);
    //     }
    // }
    for (size_t i = 1; i < tex.mipmap.size(); i++) {
        MipLevel &last_mip = tex.mipmap[i - 1], &current_mip = tex.mipmap[i];
        int h = current_mip.height, w = current_mip.width;
        int l_w = last_mip.width;
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                for (int k = 0; k < 4; k++) {
                    current_mip.texels[4 * (j + w * i) + k] = 0.25 * (last_mip.texels[4 * (2 * j + 2 * l_w * i) + k] + last_mip.texels[4 * (2 * (j + 1) + 2 * l_w * i) + k] + last_mip.texels[4 * (2 * j + 2 * l_w * (i + 1)) + k] + last_mip.texels[4 * (2 * (j + 1) + 2 * l_w * (i + 1)) + k]);
                }
            }
        }
    }
}

Color Sampler2DImp::sample_nearest(Texture& tex,
                                   float u, float v,
                                   int level) {
    // Task 6: Implement nearest neighbour interpolation

    // return magenta for invalid level
    if (level >= tex.mipmap.size()) {
        return Color(1, 0, 1, 1);
    } else {
        Color color;
        int w = tex.mipmap[level].width, h = tex.mipmap[level].height;
        int i = floor(v * h), j = floor(u * w);
        std::vector<unsigned char>& tex_level = tex.mipmap[level].texels;

        color.r = tex_level[4 * (j + w * i) + 0];
        color.g = tex_level[4 * (j + w * i) + 1];
        color.b = tex_level[4 * (j + w * i) + 2];
        color.a = tex_level[4 * (j + w * i) + 3];
        return color;
    }
}

Color Sampler2DImp::sample_bilinear(Texture& tex,
                                    float u, float v,
                                    int level) {
    // Task 6: Implement bilinear filtering

    if (level >= tex.mipmap.size()) {
        // return magenta for invalid level
        return Color(1, 0, 1, 1);
    } else {
        Color color;
        int w = tex.mipmap[level].width, h = tex.mipmap[level].height;
        int i = floor(v * h), j = floor(u * w);
        if (i >= h - 1 || j >= w - 1) {
            // return magenta for invalid level
            return sample_nearest(tex, u, v, level);
        }
        float tx = u * w - floor(u * w), ty = v * h - floor(v * h);

        std::vector<unsigned char>& tex_level = tex.mipmap[level].texels;

        color.r = (1 - ty) * ((1 - tx) * tex_level[4 * (j + w * i) + 0] + tx * tex_level[4 * (j + 1 + w * i) + 0]) + ty * ((1 - tx) * tex_level[4 * (j + w * (i + 1)) + 0] + tx * tex_level[4 * (j + 1 + w * (i + 1)) + 0]);
        color.g = (1 - ty) * ((1 - tx) * tex_level[4 * (j + w * i) + 1] + tx * tex_level[4 * (j + 1 + w * i) + 1]) + ty * ((1 - tx) * tex_level[4 * (j + w * (i + 1)) + 1] + tx * tex_level[4 * (j + 1 + w * (i + 1)) + 1]);
        color.b = (1 - ty) * ((1 - tx) * tex_level[4 * (j + w * i) + 2] + tx * tex_level[4 * (j + 1 + w * i) + 2]) + ty * ((1 - tx) * tex_level[4 * (j + w * (i + 1)) + 2] + tx * tex_level[4 * (j + 1 + w * (i + 1)) + 2]);
        color.a = (1 - ty) * ((1 - tx) * tex_level[4 * (j + w * i) + 3] + tx * tex_level[4 * (j + 1 + w * i) + 3]) + ty * ((1 - tx) * tex_level[4 * (j + w * (i + 1)) + 3] + tx * tex_level[4 * (j + 1 + w * (i + 1)) + 3]);
        return color;
    }
}

Color Sampler2DImp::sample_trilinear(Texture& tex,
                                     float u, float v,
                                     float u_scale, float v_scale) {
    // Task 7: Implement trilinear filtering
    // u_scale corresponds to du/dx, v_scale dv/dy (not du/dy, dv/dx since only scaling and no rotation)

    float L = max(log2(max(u_scale, v_scale)), (float)0);
    float tL = L - floor(L);

    return (1 - tL) * sample_bilinear(tex, u, v, floor(L)) + tL * sample_bilinear(tex, u, v, ceil(L));
}

}  // namespace CMU462
