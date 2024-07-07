// Geom.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include "bezier_curve.h"
#include "bspline_curve.h"
#include "bezier_curve_gpt.h"

using Point = Point3d;

int test_BSplineCurve()
{
	std::vector<Point> ctrlPts{
		Point(0,0,0),
		Point(5,10,0),
		Point(15,10,0),
		Point(20,0,0)
	};

    std::vector<double> knots{ 0,0,0,0,1,1,1,1 };

    BSplineCurve<Point> spline(3, knots, ctrlPts);
    double t = 0.65;
    Point eval = spline.EvalPointByBasis(t);
    eval = spline.EvalPointByDeBoor(t);

    double up = spline.ParamOfPoint(eval);

    std::vector<Point> deri = spline.CalDerivative(3, t);

    std::cout << deri[0];

    std::cout << "Point on BSpline curve at t = " << t << ": (" << eval.x() << ", " << eval.y() << ")\n";
    std::cout << "Point on BSpline curve by deBoor at t = " << t << ": (" << eval.x() << ", " << eval.y() << ")\n";

    return 0;
}

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

// 控制点类型
using ControlPoint = Point3d;
using ControlPoints = std::vector<ControlPoint>;


int test_insert_and_remove() {
    // 示例用法
    int degree = 2;
    std::vector<double> knotvector = { 0.0, 0.0, 0.0, 1.0,3.0, 3.0, 4.0, 4.0, 4.0 };
    ControlPoints ctrlpts = { {0.0, 0.0, 0}, {1.0, 2.0, 0}, {2.0, 3.0, 0}, {4.0, 4.0, 0}, {5.0, 2.0, 0} ,{6.0,1.0, 0} };
    double u = 2.5;

    double t = 0.3;
    BSplineCurve<Point> curv(degree, knotvector, ctrlpts);
    Point a = curv.EvalPointByBasis(t);

    ControlPoints result = curv.KnotInsertion(u);
    std::vector<double> knotvector_new = { 0.0, 0.0, 0.0, 1.0, 2.5, 3.0, 3.0, 4.0, 4.0, 4.0 };
    BSplineCurve<Point> curv_new(degree, knotvector_new, result);
    Point b = curv_new.EvalPointByBasis(t);
    // 输出结果
    for (const auto& pt : result) {
        std::cout << pt << " ";
        std::cout << std::endl;
    }

    BSplineCurve<Point> curv_remove = curv.KnotRemoval(3.0, 1);
    Point c = curv_remove.EvalPointByBasis(t);

    return 0;
}

int test_non_removable()
{
	int degree = 2;
	std::vector<double> knotvector = { 0.0, 0.0, 0.0, 1.0, 3.0, 4.0, 4.0, 4.0 };
	ControlPoints ctrlpts = { {0.0, 0.0, 0}, {1.0, 2.0, 0}, {2.0, 3.0, 0}, {4.0, 4.0, 0}, {5.0, 2.0, 0}};
    BSplineCurve<Point> curv(degree, knotvector, ctrlpts);
	BSplineCurve<Point> curv_remove = curv.KnotRemoval(3.0, 1);
	Point c = curv_remove.EvalPointByBasis(0.3);

    std::vector<BSplineCurve<Point>> bseg;
    curv.Decompose(bseg);

    BSplineCurve<Point> curv_h = curv.degreeElevate(1);
    Point d = curv_h.EvalPointByBasis(0.3);

    return 1;
}

int test_bezier_reduce()
{
	std::vector<Point> ctrlPts{
		Point(0,0,0),
		Point(5,10,0),
		Point(15,10,0),
		Point(20,0,0)
	};
   
	BeizierCurve<Point> curve(3, ctrlPts);
	std::vector<Point> elevate;
	curve.DegreeElevation(elevate);

    std::vector<double> knots{ 0,0,0,0,0,1,1,1,1,1};
	BSplineCurve<Point> curve_ele(4, knots, elevate);

    BSplineCurve<Point> curve_reduce;
    curve_ele.DegreeReduce();
    return 1;
}

int test_shit_knots()
{
	int degree = 2;
	std::vector<double> knotvector = { 0.0, 0.0, 0.0, 1.0, 3.0, 4.0, 4.0, 4.0 };
	ControlPoints ctrlpts = { {0.0, 0.0, 0}, {1.0, 2.0, 0}, {2.0, 3.0, 0}, {4.0, 4.0, 0}, {5.0, 2.0, 0} };
	BSplineCurve<Point> curv(degree, knotvector, ctrlpts);
    
    std::vector<Point> curve_tess;
    std::vector<double> us;
    curv.TessellateEqualKnot(curve_tess, us, 10);

    BSplineCurve<Point> curv_shift;
    curv.ReparamLinear(1.0, 6.0, curv_shift);
	std::vector<Point> curve_tess_shift;
    curv_shift.TessellateEqualKnot(curve_tess_shift, us, 10);
    return 1;

}

int main_()
{
    std::vector<Point> ctrlPts{
        Point(0,0,0),
        Point(5,10,0),
        Point(15,10,0),
        Point(20,0,0)
    };

    BeizierCurve<Point> curve(3, ctrlPts);

    Point eval = curve.EvalPoint(0.3);
    std::cout << "Point on Bezier curve at t = " << 0.3 << ": (" << eval.x() << ", " << eval.y() << ")\n";

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

    test_BSplineCurve();

    test_insert_and_remove();

    test_non_removable();

    test_bezier_reduce();

    test_shit_knots();

    std::cout << "Hello World!\n";
    return 1;
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
