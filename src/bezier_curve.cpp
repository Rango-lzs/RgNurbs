#include "bezier_curve.h"
#include "nurbs_export.h"
#include "rg_assert.h"

#define Degree m_dataRep->m_degree
#define CtrlPts m_dataRep->m_ctrlPts
#define Knots m_dataRep->m_knots

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
BeizierCurve<Point>::BeizierCurve()
	:m_dataRep(std::make_unique<DataRep>(0, std::vector<Point>()))
{

}

template<class Point>
BeizierCurve<Point>::BeizierCurve(int degree, std::vector<Point> ctrlPts)
	:m_dataRep(std::make_unique<DataRep>(degree, ctrlPts))
{
	RgAssert(degree == ctrlPts.size() - 1,"Invalid Parameter");
}


template<class Point>
BeizierCurve<Point>::BeizierCurve(const BeizierCurve& other)
	:m_dataRep(std::make_unique<DataRep>(other.Degree, other.CtrlPts))
{
	RgAssert(Degree == CtrlPts.size() - 1, "Invalid Parameter");
}


template<class Point>
BeizierCurve<Point>& BeizierCurve<Point>::operator=(const BeizierCurve& other)
{
	m_dataRep = std::make_unique<DataRep>(other.Degree, other.CtrlPts);
	return *this;
}


template<class Point>
BeizierCurve<Point>::~BeizierCurve()
{

}

//
template<class Point>
Point BeizierCurve<Point>::EvalPoint(double param)
{
	std::vector<double> B = calBernstein(Degree, param);
	Point temp;
	for (int k = 0; k <= Degree; k++)
	{
		temp = temp + B[k] * CtrlPts[k];
	}
	return temp;
}

template<class Point>
Point BeizierCurve<Point>::EvalPointDirect(double t)
{
	Point temp;
	int n = Degree;
	for (int k = 0; k <= Degree; k++)
	{
		temp = temp + comb(n, k) * pow(t,k) * (pow((1 - t) ,(n - k))) * CtrlPts[k];
	}
	return temp;
}

//递推三角形
template<class Point>
Point BeizierCurve<Point>::EvalPointByDeCasteljau(double param)
{
	std::vector<Point> temp = CtrlPts;
	for (int k = 1; k <= Degree; k++)
	{
		for (int i = 0; i <= Degree - k; i++)
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
	int n = Degree;
	std::vector<Point> firstDeriCtrlPts(n);
	for (int i = 0; i < n; ++i) {
		Point dp = n * (CtrlPts[i + 1] + (-1.0) * CtrlPts[i]);

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
	left.resize(Degree + 1, Point());
	right.resize(Degree + 1, Point());

	left[0] = CtrlPts[0];
	right[Degree] = CtrlPts[Degree];

	std::vector<Point> temp = CtrlPts;
	//沿着三角形递推，每往前推一步，取首尾点
	for (int k = 1; k <= Degree; k++)
	{
		for (int i = 0; i <= Degree - k; i++)
		{
			temp[i] = (1.0 - param) * temp[i] + param * temp[i + 1];
		}

		left[k] = temp[0];
		right[Degree - k] = temp[Degree - k];
	}
	return true;
}

template<class Point>
bool BeizierCurve<Point>::SubDivision_1(double param, std::vector<Point>& cv_left, std::vector<Point>& cv_right)
{
	cv_left.resize(CtrlPts.size());
	cv_right.resize(CtrlPts.size());

	int p = CtrlPts.size() - 1;
	//将 cv 拷贝到 cv_right,因为每次迭代计算的点数依次少1，
	//因此可以直接用cv_right记录，全部算完后就是所有的尾点集合
	std::copy(CtrlPts.begin(), CtrlPts.end(), cv_right.begin());
	cv_left[0] = cv_right[0];

	//p次迭代
	for (int i = 0; i < p; ++i)
	{
		for (int j = 0; j < p - i; ++j)
			cv_right[j] = (1.0 - param) * cv_right[j] + param * cv_right[j + 1];
		//左侧控制点为每次迭代的首点
		cv_left[i + 1] = cv_right[0];
	}

	return true;
}

template<class Point>
bool BeizierCurve<Point>::DegreeElevation(std::vector<Point>& nweCtrlPts)
{
	int n = Degree + 1;
	nweCtrlPts.resize(Degree + 2);

	nweCtrlPts.front() = CtrlPts.front();
	nweCtrlPts.back() = CtrlPts.back();

	for (int i = 1; i < n; ++i)
	{
		double ratio = ((double)i) / n;
		nweCtrlPts[i] = ratio * CtrlPts[i - 1] + (1.0 - ratio) * CtrlPts[i];
	}

	return true;
}

template class RG_API BeizierCurve<Point3d>;
template class RG_API BeizierCurve<Point2d>;
