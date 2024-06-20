#ifndef BSPLINE_CURVE_HH
#define BSPLINE_CURVE_HH

#include "bspline_basis.h"
#include <vector>

//BSpline curve was define by the degree, knots and control points
/*
* 1. ��������ֵ(ֱ�ӷ� �� DeBoor��)
* 2. �������㣬֧��һ�׺Ͷ���
* 3. �ڵ����
* 4. ���ף�����
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

	//ͨ�������������, ����t:[ui,ui+1),��N(k,p) k = i, i-1, i-2, ... i-p 
	//��p+1��������Ϊ�����
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


	//ͨ��deBoor���������
	Point EvalPointByDeBoor(double t)
	{
		return Point();
	}
	
	//��һ�ף����׵���
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
