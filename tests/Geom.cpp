// Geom.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include "beizier_curve.h"

using Point = Point3d;

int test_BezierGPT() {
    // Define control points for a cubic Bezier curve
    std::vector<Point> controlPoints = {
        {0.0, 0.0, 0},
        {1.0, 2.0, 0},
        {3.0, 3.0, 0},
        {4.0, 0.0, 0}
    };

    BezierCurve_GPT<Point> bezier(controlPoints);

    double t = 0.5;
    Point point = bezier.calculate(t);
    Point firstDerivative = bezier.firstDerivative(t);
    Point secondDerivative = bezier.secondDerivative(t);

    std::cout << "Test derivative by GPT" << std::endl;
    std::cout << "Point on Bezier curve at t = " << t << ": (" << point.x() << ", " << point.y() << ")\n";
    std::cout << "First derivative of Bezier curve at t = " << t << ": (" << firstDerivative.x() << ", " << firstDerivative.y() << ")\n";
    std::cout << "Second derivative of Bezier curve at t = " << t << ": (" << secondDerivative.x() << ", " << secondDerivative.y() << ")\n";

    return 0;
}

int main()
{
    std::vector<Point> ctrlPts{
        Point(0,0,0),
        Point(5,10,0),
        Point(15,10,0),
        Point(20,0,0)
    };

    BeizierCurve<Point> curve(3, ctrlPts);

    Point eval = curve.EvalPoint(0.3);
    Point eval1 = curve.EvalPointDirect(0.3);
    Point  eval2 = curve.EvalPointByDeCasteljau(0.3);
    Point firstDeri = curve.EvalDerivation(1, 0.3);
    Point secondDeri = curve.EvalDerivation(2, 0.3);

    std::vector<Point> left;
    std::vector<Point> right;
    bool flag = curve.SubDivision(0.3, left, right);

    std::vector<Point> left1;
    std::vector<Point> right1;
    curve.SubDivision_1(0.3, left1, right1);

    std::vector<Point> elevate;
    curve.DegreeElevation(elevate);

    BeizierCurve<Point> curve_h(4, elevate);
    Point evalh = curve_h.EvalPoint(0.3);
    Point eval1h = curve_h.EvalPointDirect(0.3);
    Point  eval2h = curve_h.EvalPointByDeCasteljau(0.3);
    Point firstDerih = curve_h.EvalDerivation(1, 0.3);
    Point secondDerih = curve_h.EvalDerivation(2, 0.3);

    //test derivation
    // Define control points for a cubic Bezier curve
    std::vector<Point> controlPoints = {
        {0.0, 0.0, 0},
        {1.0, 2.0, 0},
        {3.0, 3.0, 0},
        {4.0, 0.0, 0}
    };
    BeizierCurve<Point> curve_deri(3, controlPoints);

    double t = 0.5;

    Point point = curve_deri.deCasteljau(controlPoints, t);
    Point firstDerivative = curve_deri.EvalDerivation(1, t);
    Point secondDerivative = curve_deri.EvalDerivation(2, t);

    std::cout << "Point on Bezier curve at t = " << t << ": (" << point.x() << ", " << point.y() << ")\n";
    std::cout << "First derivative of Bezier curve at t = " << t << ": (" << firstDerivative.x() << ", " << firstDerivative.y() << ")\n";
    std::cout << "Second derivative of Bezier curve at t = " << t << ": (" << secondDerivative.x() << ", " << secondDerivative.y() << ")\n";

    test_BezierGPT();
    std::cout << "Hello World!\n";
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
