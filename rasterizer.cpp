// clang-format off
//
// Created by goksu on 4/6/19.
//

#include <algorithm>
#include <vector>
#include "rasterizer.hpp"
#include <opencv2/opencv.hpp>
#include <math.h>
#include<Eigen/Core>
#include<Eigen/Dense>
#include <tuple>


rst::pos_buf_id rst::rasterizer::load_positions(const std::vector<Eigen::Vector3f> &positions)
{
    auto id = get_next_id();
    pos_buf.emplace(id, positions);

    return {id};
}

rst::ind_buf_id rst::rasterizer::load_indices(const std::vector<Eigen::Vector3i> &indices)
{
    auto id = get_next_id();
    ind_buf.emplace(id, indices);

    return {id};
}

rst::col_buf_id rst::rasterizer::load_colors(const std::vector<Eigen::Vector3f> &cols)
{
    auto id = get_next_id();
    col_buf.emplace(id, cols);

    return {id};
}

auto to_vec4(const Eigen::Vector3f& v3, float w = 1.0f)
{
    return Vector4f(v3.x(), v3.y(), v3.z(), w);
}

double product(Vector2f p1, Vector2f p2, Vector2f p3) {
    return (p2.x() - p1.x()) * (p3.y() - p1.y()) - (p2.y() - p1.y()) * (p3.x() - p1.x());
}

bool isInTriangle(Vector2f p1, Vector2f p2, Vector2f p3, Vector2f o) {
    if (product(p1, p2, p3) < 0) return isInTriangle(p1, p3, p2, o);
    if (product(p1, p2, o) > 0 && product(p2, p3, o) > 0 && product(p3, p1, o) > 0)
        return true;
    return false;
}


static bool insideTriangle(int x, int y, const Vector3f* _v)
{   
    // TODO : Implement this function to check if the point (x, y) is inside the triangle represented by _v[0], _v[1], _v[2]
    bool inside = false;
    inside = isInTriangle(Vector2f(_v[0].x(), _v[0].y()), Vector2f(_v[1].x(), _v[1].y()), Vector2f(_v[2].x(), _v[2].y()), Vector2f(x,y));
    return inside;
}

static std::tuple<float, float, float> computeBarycentric2D(float x, float y, const Vector3f* v)
{
    float c1 = (x*(v[1].y() - v[2].y()) + (v[2].x() - v[1].x())*y + v[1].x()*v[2].y() - v[2].x()*v[1].y()) / (v[0].x()*(v[1].y() - v[2].y()) + (v[2].x() - v[1].x())*v[0].y() + v[1].x()*v[2].y() - v[2].x()*v[1].y());
    float c2 = (x*(v[2].y() - v[0].y()) + (v[0].x() - v[2].x())*y + v[2].x()*v[0].y() - v[0].x()*v[2].y()) / (v[1].x()*(v[2].y() - v[0].y()) + (v[0].x() - v[2].x())*v[1].y() + v[2].x()*v[0].y() - v[0].x()*v[2].y());
    float c3 = (x*(v[0].y() - v[1].y()) + (v[1].x() - v[0].x())*y + v[0].x()*v[1].y() - v[1].x()*v[0].y()) / (v[2].x()*(v[0].y() - v[1].y()) + (v[1].x() - v[0].x())*v[2].y() + v[0].x()*v[1].y() - v[1].x()*v[0].y());
    return {c1,c2,c3};
}

void rst::rasterizer::draw(pos_buf_id pos_buffer, ind_buf_id ind_buffer, col_buf_id col_buffer, Primitive type)
{
    auto& buf = pos_buf[pos_buffer.pos_id];
    auto& ind = ind_buf[ind_buffer.ind_id];
    auto& col = col_buf[col_buffer.col_id];

    float f1 = (50 - 0.1) / 2.0;
    float f2 = (50 + 0.1) / 2.0;

    Eigen::Matrix4f mvp = projection * view * model;
    for (auto& i : ind)
    {
        Triangle t;
        Eigen::Vector4f v[] = {
                mvp * to_vec4(buf[i[0]], 1.0f),
                mvp * to_vec4(buf[i[1]], 1.0f),
                mvp * to_vec4(buf[i[2]], 1.0f)
        };
        //Homogeneous division
        for (auto& vec : v) {
            vec /= vec.w();
        }
        //Viewport transformation
        for (auto & vert : v)
        {
            vert.x() = 0.5*width*(vert.x()+1.0);
            vert.y() = 0.5*height*(vert.y()+1.0);
            vert.z() = vert.z() * f1 + f2;
        }

        for (int i = 0; i < 3; ++i)
        {
            t.setVertex(i, v[i].head<3>());
            t.setVertex(i, v[i].head<3>());
            t.setVertex(i, v[i].head<3>());
        }

        auto col_x = col[i[0]];
        auto col_y = col[i[1]];
        auto col_z = col[i[2]];

        t.setColor(0, col_x[0], col_x[1], col_x[2]);
        t.setColor(1, col_y[0], col_y[1], col_y[2]);
        t.setColor(2, col_z[0], col_z[1], col_z[2]);

        rasterize_triangle(t);
    }
}

void rst::rasterizer::draw_point(const Eigen::Vector3f& point, const Eigen::Vector3f& color, int size) {

    Eigen::Vector3f point00 = Eigen::Vector3f(point.x()-size, point.y()-size, 1.0f);

    for (size_t i = 0; i < 2*size; i++)
    {
        for (size_t j = 0; j < 2*size; j++)
        {
            set_pixel(Eigen::Vector3f(point00.x()+i, point00.y()+j, 1.0f), color);
        }
    }
}

//Screen space rasterization
void rst::rasterizer::rasterize_triangle(const Triangle& t) {
    draw_line(t.v[2], t.v[0]);
    draw_line(t.v[2], t.v[1]);
    draw_line(t.v[1], t.v[0]);

    auto v = t.toVector4();
    
    int minx = width-1, miny = height-1, maxx = 0, maxy = 0;
    for (auto& vert : v)
    {
        if (vert.x() < minx) {
            minx = vert.x();
        }

        if (vert.y() < miny) {
            miny = vert.y();
        }

        if (vert.x() > maxx) {
            maxx = vert.x();
        }

        if (vert.y() > maxy) {
            maxy = vert.y();
        }
    }

    minx = MIN(MAX(0, minx), width - 1);  // 范围0-width-1
    miny = MIN(MAX(0, miny), height - 1);  // 范围0-height-1
    maxx = MAX(MIN(maxx, width - 1),0);
    maxy = MAX(MIN(maxy, height - 1), 0);

    draw_point(Eigen::Vector3f(minx, miny, 1.0f), Eigen::Vector3f(t.color[0].x()*255, t.color[0].y() * 255, t.color[0].z() * 255),5);
    draw_point(Eigen::Vector3f(minx, maxy, 1.0f), Eigen::Vector3f(t.color[0].x() * 255, t.color[0].y() * 255, t.color[0].z() * 255), 5);
    draw_point(Eigen::Vector3f(maxx, miny, 1.0f), Eigen::Vector3f(t.color[0].x() * 255, t.color[0].y() * 255, t.color[0].z() * 255), 5);
    draw_point(Eigen::Vector3f(maxx, maxy, 1.0f), Eigen::Vector3f(t.color[0].x() * 255, t.color[0].y() * 255, t.color[0].z() * 255), 5);


    for (size_t i = minx; i <= maxx; i++)
    {
        for (size_t j = miny; j <= maxy; j++) {
            if (insideTriangle(i+0.5, j+0.5, t.v)) {
                auto my_turple = computeBarycentric2D(i + 0.5, j + 0.5, t.v);
                float alpha=std::get<0>(my_turple);
                float beta = std::get<0>(my_turple);
                float gamma = std::get<0>(my_turple);
                float w_reciprocal = 1.0/(alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
                float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
                z_interpolated *= w_reciprocal;
                float cur_z = depth_buf[get_index(i,j)];
                if (z_interpolated < cur_z)
                {
                    set_pixel(Eigen::Vector3f(i, j, 1.0f), Eigen::Vector3f(t.color[0].x() * 255, t.color[0].y() * 255, t.color[0].z() * 255));
                    depth_buf[get_index(i, j)] = z_interpolated;
                }
            }
        }
    }
    

    // TODO : Find out the bounding box of current triangle.
    // iterate through the pixel and find if the current pixel is inside the triangle

    // If so, use the following code to get the interpolated z value.
    //auto[alpha, beta, gamma] = computeBarycentric2D(x, y, t.v);
    //float w_reciprocal = 1.0/(alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
    //float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
    //z_interpolated *= w_reciprocal;

    // TODO : set the current pixel (use the set_pixel function) to the color of the triangle (use getColor function) if it should be painted.
}

void rst::rasterizer::set_model(const Eigen::Matrix4f& m)
{
    model = m;
}

void rst::rasterizer::set_view(const Eigen::Matrix4f& v)
{
    view = v;
}

void rst::rasterizer::set_projection(const Eigen::Matrix4f& p)
{
    projection = p;
}

void rst::rasterizer::clear(rst::Buffers buff)
{
    if ((buff & rst::Buffers::Color) == rst::Buffers::Color)
    {
        std::fill(frame_buf.begin(), frame_buf.end(), Eigen::Vector3f{0, 0, 0});
    }
    if ((buff & rst::Buffers::Depth) == rst::Buffers::Depth)
    {
        std::fill(depth_buf.begin(), depth_buf.end(), std::numeric_limits<float>::infinity());
    }
}

rst::rasterizer::rasterizer(int w, int h) : width(w), height(h)
{
    frame_buf.resize(w * h);
    depth_buf.resize(w * h);
}

int rst::rasterizer::get_index(int x, int y)
{
    return (height-1-y)*width + x;
}

void rst::rasterizer::set_pixel(const Eigen::Vector3f& point, const Eigen::Vector3f& color)
{
    //old index: auto ind = point.y() + point.x() * width;
    if (point.x() < 0 || point.x() >= width ||
        point.y() < 0 || point.y() >= height) {
        return;
    }
    auto ind = (height-1-point.y())*width + point.x();
    frame_buf[ind] = color;

}

// Bresenham's line drawing algorithm
// Code taken from a stack overflow answer: https://stackoverflow.com/a/16405254
void rst::rasterizer::draw_line(Eigen::Vector3f begin, Eigen::Vector3f end)
{
    auto x1 = begin.x();
    auto y1 = begin.y();
    auto x2 = end.x();
    auto y2 = end.y();

    Eigen::Vector3f line_color = { 255, 255, 255 };

    int x, y, dx, dy, dx1, dy1, px, py, xe, ye, i;

    dx = x2 - x1;
    dy = y2 - y1;
    dx1 = fabs(dx);//x的间隔长度
    dy1 = fabs(dy);
    px = 2 * dy1 - dx1;
    py = 2 * dx1 - dy1;

    if (dy1 <= dx1) //斜率小于或等于1
    {
        if (dx >= 0) //x1 ---------x2
        {
            x = x1;
            y = y1;
            xe = x2;
        }
        else
        {        //x2 ------x1
            x = x2;
            y = y2;
            xe = x1;
        }
        Eigen::Vector3f point = Eigen::Vector3f(x, y, 1.0f);
        set_pixel(point, line_color);
        for (i = 0; x < xe; i++)
        {
            x = x + 1;
            if (px < 0)  //2倍y的间隔小于x间隔
            {
                px = px + 2 * dy1;
            }
            else
            {
                if ((dx < 0 && dy < 0) || (dx > 0 && dy > 0))
                {
                    y = y + 1;
                }
                else
                {
                    y = y - 1;
                }
                px = px + 2 * (dy1 - dx1);
            }
            //            delay(0);
            Eigen::Vector3f point = Eigen::Vector3f(x, y, 1.0f);
            set_pixel(point, line_color);
        }
    }
    else//斜率大于1
    {
        if (dy >= 0)  // y2  --- y1
        {
            x = x1;
            y = y1;
            ye = y2;
        }
        else
        {
            x = x2;
            y = y2;
            ye = y1;
        }
        Eigen::Vector3f point = Eigen::Vector3f(x, y, 1.0f);
        set_pixel(point, line_color);
        for (i = 0; y < ye; i++)
        {
            y = y + 1;
            if (py <= 0)
            {
                py = py + 2 * dx1;
            }
            else
            {
                if ((dx < 0 && dy < 0) || (dx > 0 && dy > 0))
                {
                    x = x + 1;
                }
                else
                {
                    x = x - 1;
                }
                py = py + 2 * (dx1 - dy1);
            }
            //            delay(0);
            Eigen::Vector3f point = Eigen::Vector3f(x, y, 1.0f);
            set_pixel(point, line_color);
        }
    }
}

// clang-format on