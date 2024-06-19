#ifndef BSPLINE_CURVE_HH
#define BSPLINE_CURVE_HH

#include <vector>

//BSpline curve was define by the degree, knots and control points
/*
* 1. 参数点求值(直接法 和 DeBoor法)
* 2. 导数计算，支持一阶和二阶
* 3. 节点插入
* 4. 升阶，降阶
* 5. Subdivide
*/

template<class Point>
class BSplineCurve
{
public:
	BSplineCurve(int degree, const std::vector<double>& knots, const std::vector<Point>& ctrlPts);

	//通过基函数计算点
	EvalPointByBasis(double t);
	{

	}


	//通过deBoor方法计算点
	EvalPointByDeBoor(double t)
	{

	}
	
	//求一阶，二阶导数
	Point CalDerivative(int order)
	{

	}

private:
	int m_degree;
	std::vector<double> m_knots;	//size m
	std::vector<Point> m_ctrl_pts;  //size n    m = n + m_degree +1;
};

#endif
