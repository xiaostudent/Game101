#include<cmath>
#include<Eigen/Core>
#include<Eigen/Dense>
#include<iostream>
#include <opencv2/opencv.hpp>
using namespace cv;
int main() {

    Mat image = imread("D:/game/cocos2d-x/test/Game101/terrain5.png");
    if (image.empty()) {
        printf("could not load image...\n");
        return -1;

    }
    namedWindow("test_opencv_setup", 0);
    imshow("test_opencv_srtup", image);
    waitKey(0);
    return 0;
}