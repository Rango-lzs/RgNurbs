#ifndef BSPLINE_CURVE_HH
#define BSPLINE_CURVE_HH

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
	BSplineCurve(int degree, const std::vector<double>& knots, const std::vector<Point>& ctrlPts);

	//ͨ�������������
	EvalPointByBasis(double t);
	{

	}


	//ͨ��deBoor���������
	EvalPointByDeBoor(double t)
	{

	}
	
	//��һ�ף����׵���
	Point CalDerivative(int order)
	{

	}

private:
	int m_degree;
	std::vector<double> m_knots;	//size m
	std::vector<Point> m_ctrl_pts;  //size n    m = n + m_degree +1;
};

#endif
