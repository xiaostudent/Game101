#include<cmath>
#include<Eigen/Core>
#include<Eigen/Dense>
#include<iostream>
#include "src/base/FileUtil.h"
#include <opencv2/opencv.hpp>
#include "src/assignment/1/Assignment1.h"
using namespace cv;
int main() {
    /*Mat image = imread(FileUtil::getRootPath().append("res/awesomeface.png"));
    if (image.empty()) {
        printf("could not load image...\n");
        return -1;

    }
    namedWindow("test_opencv_setup", 0);
    imshow("test_opencv_srtup", image);
    waitKey(0);*/

    auto assignment = new Assignment1();
    assignment->run();
    return 0;
}