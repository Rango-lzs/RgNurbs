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
		, m_knots(knots)
		, m_ctrl_pts(ctrlPts)
	{

	}

	//ͨ�������������, ����t:[ui,ui+1),��N(k,p) k = i, i-1, i-2, ... i-p 
	//��p+1��������Ϊ�����
	Point EvalPointByBasis(double t)
	{
		Point pt;
		int span = BSplineBasis::FindSpan(m_knots, m_degree, t);
		std::vector<double> N = BSplineBasis::BasisFunctions(span, m_degree, m_knots, t);

		for (int i = 0; i <= m_degree; i++)
		{
			pt = pt + N[i] * m_ctrl_pts[span - m_degree + i];
		}
		return pt;
	}


	//ref: https://en.wikipedia.org/wiki/De_Boor%27s_algorithm
	/*Evaluates S(x).
	Arguments
		-------- -
		span : Index of knot interval that contains x.
		u : Position.
		knots : Array of knot positions, needs to be padded as described above.
		c : Array of control points.
		p : Degree of B - spline.
	*/
	//����Beizier��deCas�㷨����ͣ�Ľ���cut cornor������deCas�������Կ�����DeBoor����������
	//Deboor�������и������ÿ�������Ǳ仯�ģ�deCas���䣬��������Beizie�Ľڵ��������ظ��ģ����������
	//�и�ϵ����һ����
	Point deBoor(double u)
	{
		int span = BSplineBasis::FindSpan(m_knots, m_degree, u);
		int p = m_degree;
		std::vector<Point> d(p + 1);
		for (int i = 0; i < p +1 ; i++)
		{
			d[i] = m_ctrl_pts[i + span - p];
		}

		for (int r = 1; r < p + 1; ++r)
		{
			for (int j = p; j > r - 1; --j)
			{
				double alpha = (u - m_knots[j + span - p]) / (m_knots[j + 1 + span - r] - m_knots[j + span - p]);
				d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j];
			}
		}

		return d[p];
	}

	//ͨ��deBoor���������
	Point EvalPointByDeBoor(double t)
	{
		return Point(0,0,0);
	}

	//��һ�ף����׵���
	std::vector<Point> CalDerivative(int order, double u)
	{
		std::vector<Point> ptRet(order +1);
		int span = BSplineBasis::FindSpan(m_knots, m_degree, u);
		std::vector<std::vector<double>> dbdu = BSplineBasis::BasisFunctionsDerivatives(span, m_degree, order, m_knots, u);

		order = std::min(order, m_degree);
		for (int i = 0; i <= order; ++i)
		{
			Point du;
			int p = m_degree;
			for (int j = 0; j < p+1; ++j)
			{
				du = du + m_ctrl_pts[span - p + j] * dbdu[i][j];
			}
			ptRet[i] = du;
		}

		return ptRet;
	}

private:
	int m_degree;
	std::vector<double> m_knots;	//size m
	std::vector<Point> m_ctrl_pts;  //size n    m = n + m_degree +1;
};

#endif
