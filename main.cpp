#include <iostream>
#include <Eigen\Dense>

using namespace std;

typedef Eigen::Matrix<int, 3, 3> Matrix3i;

int main()
{
    /*
    Matrix的初始化方法
    Eigen::Matrix<int, 3, 3>
    int 代表Matrix的数据类型，3，3 分别代表 rows， cols
    Matrix3i m1;
    m1(0,0) = 1
    m1(0,1) = 2
    m1(0,2) = 3
    ...
    或者用 m1 << 1,2,3 ...
    */

    Matrix3i m1;
    m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    cout << "m1 = \n" << m1 << endl;

    Matrix3i m2;
    m2 << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    cout << "m2 = \n" << m2 << endl;

    cout << "m1 * m2 = \n" << (m1 * m2) << endl;

    while (true) {

    }

    return 0;
}