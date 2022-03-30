#include "rasterizer.h"

rasterizer::rasterizer(int w, int h) :width(w), height(h){
	frame_buf.resize(w * h);
	depth_buf.resize(w * h);
}

void rasterizer::clear(Buffers buff) {
    if ((buff & Buffers::Color) == Buffers::Color)
    {
        std::fill(frame_buf.begin(), frame_buf.end(), Eigen::Vector3f{ 0, 0, 255 });
    }
    if ((buff & Buffers::Depth) == Buffers::Depth)
    {
        std::fill(depth_buf.begin(), depth_buf.end(), std::numeric_limits<float>::infinity());
    }
}

void rasterizer::set_pixel(const Eigen::Vector3f& point, const Eigen::Vector3f& color)
{
    if (point.x() < 0 || point.x() >= width ||
        point.y() < 0 || point.y() >= height) return;
    auto ind = (height - 1 - point.y()) * width + point.x();
    frame_buf[ind] = color;
}