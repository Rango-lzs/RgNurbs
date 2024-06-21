// Geom.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include "beizier_curve.h"
#include "bspline_curve.h"

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
    double t = 0.3;
    Point eval = spline.EvalPointByBasis(t);
    eval = spline.deBoor(t);

    std::vector<Point> deri = spline.CalDerivative(3, t);

    std::cout << deri;

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


#include <vector>
#include <iostream>
#include <algorithm>

// 用于深度拷贝的工具函数
template <typename T>
T deep_copy(const T& obj) {
    return obj;
}

// 控制点类型
using ControlPoint = Point3d;
using ControlPoints = std::vector<ControlPoint>;

// 查找重数函数 (伪代码实现，请根据实际情况调整)
int find_multiplicity(double u, const std::vector<double>& knotvector) {
    return std::count(knotvector.begin(), knotvector.end(), u);
}

// 查找跨度函数 (伪代码实现，请根据实际情况调整)
int find_span_linear(int degree, const std::vector<double>& knotvector, int num_ctrlpts, double u) {
    for (int i = degree; i < num_ctrlpts; ++i) {
        if (u >= knotvector[i] && u < knotvector[i + 1]) {
            return i;
        }
    }
    return num_ctrlpts - 1;
}

// 计算alpha值的函数 (伪代码实现，请根据实际情况调整)
double knot_insertion_alpha(double u, const std::vector<double>& knotvector, int degree, int k, int i, int L) {
    return (u - knotvector[k - degree + 1 + i]) / (knotvector[L + 1 + i] - knotvector[k - degree + 1 + i]);
}

// 节点插入算法
ControlPoints knot_insertion(int degree, const std::vector<double>& knotvector, const ControlPoints& ctrlpts, double u, int num = 1, int s = -1, int k = -1) {
    // 若未提供s和k参数，则计算它们
    if (s == -1) {
        s = find_multiplicity(u, knotvector);
    }
    if (k == -1) {
        k = find_span_linear(degree, knotvector, ctrlpts.size(), u);
    }

    // 初始化变量
    int np = ctrlpts.size();
    int nq = np + num;
    ControlPoints ctrlpts_new(nq);

    // 初始化一个长度为degree + 1的局部数组
    ControlPoints temp(degree + 1);

    // 保存未更改的控制点
    for (int i = 0; i < k - degree + 1; ++i) {
        ctrlpts_new[i] = ctrlpts[i];
    }
    for (int i = k - s; i < np; ++i) {
        ctrlpts_new[i + num] = ctrlpts[i];
    }

    // 开始填充将用于在节点插入期间更新控制点的临时局部数组
    for (int i = 0; i <= degree - s; ++i) {
        temp[i] = deep_copy(ctrlpts[k - degree + i]);
    }

    // 插入节点 "num" 次
    for (int j = 1; j <= num; ++j) {
        int L = k - degree + j;
        for (int i = 0; i <= degree - j - s; ++i) {
            double alpha = knot_insertion_alpha(u, knotvector,degree, k, i, L);
            temp[i] = alpha * temp[i + 1] + (1.0 - alpha) * temp[i];
            
        }
        ctrlpts_new[L] = deep_copy(temp[0]);
        ctrlpts_new[k + num - j - s] = deep_copy(temp[degree - j - s]);
    }

    // 加载剩余的控制点
    int L = k - degree + num;
    for (int i = L + 1; i < k - s; ++i) {
        ctrlpts_new[i] = deep_copy(temp[i - L]);
    }

    // 返回插入节点后的控制点
    return ctrlpts_new;
}

int insert() {
    // 示例用法
    int degree = 3;
    std::vector<double> knotvector = { 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0 };
    ControlPoints ctrlpts = { {0.0, 0.0, 0}, {1.0, 2.0, 0}, {2.0, 3.0, 0}, {4.0, 4.0, 0}, {5.0, 2.0, 0} };
    double u = 2.5;

    ControlPoints result = knot_insertion(degree, knotvector, ctrlpts, u);

    // 输出结果
    for (const auto& pt : result) {
        std::cout << pt.x() << " ";
        std::cout << std::endl;
    }

    return 0;
}