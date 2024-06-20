#ifndef BSPLINE_CURVE_HH
#define BSPLINE_CURVE_HH

#include "bspline_basis.h"
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
	BSplineCurve(int degree, const std::vector<double>& knots, const std::vector<Point>& ctrlPts)
		:m_degree(degree)
		,m_knots(knots)
		,m_ctrl_pts(ctrlPts)
	{

	}

	//通过基函数计算点, 对于t:[ui,ui+1),仅N(k,p) k = i, i-1, i-2, ... i-p 
	//共p+1个基函数为非零的
	Point EvalPointByBasis(double t)
	{
		Point pt;
		int span = BSplineBasis::FindSpan(m_knots, m_degree, t);
		std::vector<double> N = BSplineBasis::BasisFunctions(span,m_degree, m_knots, t);

		for (int i = 0; i <= m_degree; i++)
		{
			pt = pt+ N[i] * m_ctrl_pts[span - m_degree+ i];
		}
		return pt;
	}


	//通过deBoor方法计算点
	Point EvalPointByDeBoor(double t)
	{
		return Point();
	}
	
	//求一阶，二阶导数
	Point CalDerivative(int order)
	{
		return Point();
	}

private:
	int m_degree;
	std::vector<double> m_knots;	//size m
	std::vector<Point> m_ctrl_pts;  //size n    m = n + m_degree +1;
};

#endif
