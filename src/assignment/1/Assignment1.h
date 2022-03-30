#pragma once

#ifndef ASSIGNMENT1
#define ASSIGNMENT1

#include "base/Assignment.h"
#include "rasterizer.h"
#include <opencv2/opencv.hpp>

class Assignment1 : Assignment
{
public:
	Assignment1();
	~Assignment1();
	void  run() override;
};

Assignment1::Assignment1()
{
}

Assignment1::~Assignment1()
{

}

void Assignment1::run() {
    int key = 0;
    rasterizer r(700, 700);
    while (key != 27) {
        r.clear(Buffers::Color | Buffers::Depth);

        Eigen::Vector3f point = Eigen::Vector3f(1, 1, 1.0f);
        r.set_pixel(point, Eigen::Vector3f{ 255, 255, 255 });
        r.set_pixel(Eigen::Vector3f(1, 2, 1.0f), Eigen::Vector3f{ 255, 255, 255 });
        r.set_pixel(Eigen::Vector3f(1, 3, 1.0f), Eigen::Vector3f{ 255, 255, 255 });
        r.set_pixel(Eigen::Vector3f(1, 699, 1.0f), Eigen::Vector3f{ 255, 255, 255 });

        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::imshow("image", image);
        key = cv::waitKey(10);
    }
}

#endif // ! ASSIGNMENT1

