#include "beizier_curve.h"
#include "nurbs_export.h"
#include <assert.h>


namespace
{
	int comb(int n, int r)
	{
		if (r == 0 || r == n)
			return 1;
		else
			return comb(n - 1, r - 1) + comb(n - 1, r);
	}

	std::vector<double> calBernstein(int degree, double param)
	{
		std::vector<double> B(degree + 1);
		B[0] = 1.0;

		double t1 = 1.0 - param;
		for (int j = 1; j <= degree; j++)
		{
			double saved = 0.0;
			for (int k = 0; k < j; k++)
			{
				double temp = B[k];
				B[k] = saved + t1 * temp;
				saved = param * temp;
			}
			B[j] = saved;
		}
		return B;
	}

	// Linearly interpolate between two points
	template<class Point>
	Point lerp(const Point& p0, const Point& p1, double t) {
		return (1 - t) * p0 + t * p1;
	}
}

template<class Point>
struct BeizierCurve<Point>::DataRep
{
public:
	DataRep(int degree, const std::vector<Point>& ctrlPts)
		:m_degree(degree)
		, m_ctrlPts(ctrlPts)
	{

	}
	int m_degree;
	std::vector<Point> m_ctrlPts;
};

template<class Point>
BeizierCurve<Point>::BeizierCurve(int degree, std::vector<Point> ctrlPts)
	:m_degree(degree)
	,m_ctrlPts(ctrlPts)
{
	assert(degree == ctrlPts.size() - 1);
}

//
template<class Point>
Point BeizierCurve<Point>::EvalPoint(double param)
{
	std::vector<double> B = calBernstein(m_degree, param);
	Point temp;
	for (int k = 0; k <= m_degree; k++)
	{
		temp = temp + B[k] * m_ctrlPts[k];
	}
	return temp;
}

template<class Point>
Point BeizierCurve<Point>::EvalPointDirect(double t)
{
	Point temp;
	int n = m_degree;
	for (int k = 0; k <= m_degree; k++)
	{
		temp = temp + comb(n, k) * pow(t,k) * (pow((1 - t) ,(n - k))) * m_ctrlPts[k];
	}
	return temp;
}

//����������
template<class Point>
Point BeizierCurve<Point>::EvalPointByDeCasteljau(double param)
{
	std::vector<Point> temp = m_ctrlPts;
	for (int k = 1; k <= m_degree; k++)
	{
		for (int i = 0; i <= m_degree - k; i++)
		{
			temp[i] = (1.0 - param) * temp[i] + param * temp[i + 1];
		}
	}
	return temp[0];
}

// De Casteljau's algorithm to compute a point on a Bezier curve
template<class Point>
Point BeizierCurve<Point>::deCasteljau(const std::vector<Point>& controlPoints, double t) {
	std::vector<Point> points = controlPoints;
	int n = points.size();

	for (int r = 1; r < n; ++r) {
		for (int i = 0; i < n - r; ++i) {
			points[i] = lerp(points[i], points[i + 1], t);
		}
	}

	return points[0];
}

// Compute the derivative up to second order using De Casteljau's algorithm
template<class Point>
Point BeizierCurve<Point>::EvalDerivation(int order, double t)
{
	int n = m_degree;
	std::vector<Point> firstDeriCtrlPts(n);
	for (int i = 0; i < n; ++i) {
		Point dp = n * (m_ctrlPts[i + 1] + (-1.0) * m_ctrlPts[i]);

		firstDeriCtrlPts[i] = dp;
	}
	if (order == 1)
	{
		return deCasteljau(firstDeriCtrlPts, t);
	}

	int m = firstDeriCtrlPts.size() - 1;
	std::vector<Point> secondDeriCtrlPts(m);

	for (int i = 0; i < m; ++i) {
		Point ddp = m * (firstDeriCtrlPts[i + 1] + (-1.0) * firstDeriCtrlPts[i]);
		secondDeriCtrlPts[i] = ddp;
	}

	if (order == 2)
	{
		return deCasteljau(secondDeriCtrlPts, t);
	}

	return Point();
}

template<class Point>
bool BeizierCurve<Point>::SubDivision(double param, std::vector<Point>& left, std::vector<Point>& right)
{
	left.resize(m_degree + 1, Point());
	right.resize(m_degree + 1, Point());

	left[0] = m_ctrlPts[0];
	right[m_degree] = m_ctrlPts[m_degree];

	std::vector<Point> temp = m_ctrlPts;
	//���������ε��ƣ�ÿ��ǰ��һ����ȡ��β��
	for (int k = 1; k <= m_degree; k++)
	{
		for (int i = 0; i <= m_degree - k; i++)
		{
			temp[i] = (1.0 - param) * temp[i] + param * temp[i + 1];
		}

		left[k] = temp[0];
		right[m_degree - k] = temp[m_degree - k];
	}
	return true;
}

template<class Point>
bool BeizierCurve<Point>::SubDivision_1(double param, std::vector<Point>& cv_left, std::vector<Point>& cv_right)
{
	cv_left.resize(m_ctrlPts.size());
	cv_right.resize(m_ctrlPts.size());

	int p = m_ctrlPts.size() - 1;
	//�� cv ������ cv_right,��Ϊÿ�ε�������ĵ���������1��
	//��˿���ֱ����cv_right��¼��ȫ�������������е�β�㼯��
	std::copy(m_ctrlPts.begin(), m_ctrlPts.end(), cv_right.begin());
	cv_left[0] = cv_right[0];

	//p�ε���
	for (int i = 0; i < p; ++i)
	{
		for (int j = 0; j < p - i; ++j)
			cv_right[j] = (1.0 - param) * cv_right[j] + param * cv_right[j + 1];
		//�����Ƶ�Ϊÿ�ε������׵�
		cv_left[i + 1] = cv_right[0];
	}

	return true;
}

template<class Point>
bool BeizierCurve<Point>::DegreeElevation(std::vector<Point>& nweCtrlPts)
{
	int n = m_degree + 1;
	nweCtrlPts.resize(m_degree + 2);

	nweCtrlPts.front() = m_ctrlPts.front();
	nweCtrlPts.back() = m_ctrlPts.back();

	for (int i = 1; i < n; ++i)
	{
		double ratio = ((double)i) / n;
		nweCtrlPts[i] = ratio * m_ctrlPts[i - 1] + (1.0 - ratio) * m_ctrlPts[i];
	}

	return true;
}

template class RG_API BeizierCurve<Point3d>;
template class RG_API BeizierCurve<Point2d>;
