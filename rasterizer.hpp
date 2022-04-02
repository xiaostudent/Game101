//
// Created by goksu on 4/6/19.
//

#pragma once

#include<Eigen/Core>
#include<Eigen/Dense>
#include <algorithm>
#include "global.hpp"
#include "Triangle.hpp"
#include <algorithm>
#include <map>
using namespace Eigen;

namespace rst
{
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

    /*
     * For the curious : The draw function takes two buffer id's as its arguments. These two structs
     * make sure that if you mix up with their orders, the compiler won't compile it.
     * Aka : Type safety
     * */
    struct pos_buf_id
    {
        int pos_id = 0;
    };

    struct ind_buf_id
    {
        int ind_id = 0;
    };

    struct col_buf_id
    {
        int col_id = 0;
    };

    class rasterizer
    {
    public:
        rasterizer(int w, int h);
        pos_buf_id load_positions(const std::vector<Eigen::Vector3f>& positions);
        ind_buf_id load_indices(const std::vector<Eigen::Vector3i>& indices);
        col_buf_id load_colors(const std::vector<Eigen::Vector3f>& colors);

        void set_model(const Eigen::Matrix4f& m);
        void set_view(const Eigen::Matrix4f& v);
        void set_projection(const Eigen::Matrix4f& p);

        void set_pixel(const Eigen::Vector3f& point, const Eigen::Vector3f& color, int samplingIndex);  //hjx
        Eigen::Vector3f& get_pixel(const Eigen::Vector3f& point, int samplingIndex);  //hjx

        void clear(Buffers buff);

        void draw(pos_buf_id pos_buffer, ind_buf_id ind_buffer, col_buf_id col_buffer, Primitive type);

        std::vector<Eigen::Vector3f>& frame_buffer() {
            Eigen::Vector3f tmp_color;
            for (size_t i = 0; i < width; i++)  //hjx
            {
                for (size_t j = 0; j < height; j++)
                {
                    int buf_index = get_index(i, j);
                    frame_buf[buf_index].x() = 0;
                    frame_buf[buf_index].y() = 0;
                    frame_buf[buf_index].z() = 0;
                    for (size_t k = 0; k < samplingCount; k++)
                    {
                        tmp_color = get_pixel(Vector3f(i, j,0), k);
                        frame_buf[buf_index].x() = frame_buf[buf_index].x() + tmp_color.x();
                        frame_buf[buf_index].y() = frame_buf[buf_index].y() + tmp_color.y();
                        frame_buf[buf_index].z() = frame_buf[buf_index].z() + tmp_color.z();
                    }
                    frame_buf[buf_index] /= samplingCount;
                }
            }
            return frame_buf; 
        }

    private:
        void draw_line(Eigen::Vector3f begin, Eigen::Vector3f end);
        void draw_point(const Eigen::Vector3f& point, const Eigen::Vector3f& color, int size);

        void rasterize_triangle(const Triangle& t);

        float getDepth(const Eigen::Vector3f& point, int samplingIndex);  //hjx
        void  setDepth(const Eigen::Vector3f& point, const float depth, int samplingIndex);  //hjx
        // VERTEX SHADER -> MVP -> Clipping -> /.W -> VIEWPORT -> DRAWLINE/DRAWTRI -> FRAGSHADER

    private:
        Eigen::Matrix4f model;
        Eigen::Matrix4f view;
        Eigen::Matrix4f projection;

        std::map<int, std::vector<Eigen::Vector3f>> pos_buf;
        std::map<int, std::vector<Eigen::Vector3i>> ind_buf;
        std::map<int, std::vector<Eigen::Vector3f>> col_buf;

        std::vector<Eigen::Vector3f> frame_buf;
        std::vector<Eigen::Vector3f> sampling_frame_buf;  //hjx
        std::vector<float> depth_buf;
        int get_index(int x, int y);

        int samplingCount = 4;  //hjx
        int width, height;
        std::vector<Eigen::Vector3f> sampling_list;

        int next_id = 0;
        int get_next_id() { return next_id++; }
    };
}
