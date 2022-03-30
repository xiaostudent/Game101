#pragma once

#ifndef RASTERIZER
#define RASTERIZER
#include <vector>
#include<Eigen/Core>
#include<Eigen/Dense>

enum class Buffers
{
    Color = 1,
    Depth = 2
};

inline Buffers operator|(Buffers a, Buffers b)
{
    return Buffers((int)a | (int)b);
}

inline Buffers operator&(Buffers a, Buffers b)
{
    return Buffers((int)a & (int)b);
}

enum class Primitive
{
    Line,
    Triangle
};


class rasterizer {

public:
	rasterizer(int width, int height);
	void clear(Buffers buff);
    std::vector<Eigen::Vector3f>& frame_buffer() { return frame_buf; }
    void set_pixel(const Eigen::Vector3f& point, const Eigen::Vector3f& color);
private:
	std::vector<Eigen::Vector3f> frame_buf;
	std::vector<float> depth_buf;
    int width, height;
};


#endif
