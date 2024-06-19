#ifndef BEIZIER_CURVE_HH
#define BEIZIER_CURVE_HH

#include "point.h"
#include <vector>

template<class T>
class BeizierCurve
{
	static_assert(std::is_same_v<T, Point3d>||std::is_same_v<T, Point2d>, "only support point3d and point2d");
	using Point = T;
public:
	//beizier curve was total defined by control points
	BeizierCurve(int degree, std::vector<Point> ctrlPts);

    //
	Point EvalPoint(double param);

	Point EvalPointDirect(double param);

	//递推三角形
	Point EvalPointByDeCasteljau(double param)
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

	// Linearly interpolate between two points
	Point lerp(const Point& p0, const Point& p1, double t) {
		return (1 - t) * p0 + t * p1;
	}

	// De Casteljau's algorithm to compute a point on a Bezier curve
	Point deCasteljau(const std::vector<Point>& controlPoints, double t) {
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
	Point EvalDerivation(int order, double t)
	{		
			int n = m_degree;
			std::vector<Point> firstDeriCtrlPts(n);
			for (int i = 0; i < n; ++i) {
				Point dp = n *( m_ctrlPts[i + 1] + (-1.0)*m_ctrlPts[i]);

				firstDeriCtrlPts[i] =  dp;
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

	bool SubDivision(double param, std::vector<Point>& left, std::vector<Point>& right)
	{
		left.resize(m_degree + 1, Point());
		right.resize(m_degree + 1, Point());

		left[0] = m_ctrlPts[0];
		right[m_degree] = m_ctrlPts[m_degree];

		std::vector<Point> temp = m_ctrlPts;
		//沿着三角形递推，每往前推一步，取首尾点
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

	bool SubDivision_1(double param, std::vector<Point>& cv_left, std::vector<Point>& cv_right)
	{
		cv_left.resize(m_ctrlPts.size());
		cv_right.resize(m_ctrlPts.size());

		int p = m_ctrlPts.size() - 1;
		//将 cv 拷贝到 cv_right,因为每次迭代计算的点数依次少1，
		//因此可以直接用cv_right记录，全部算完后就是所有的尾点集合
		std::copy(m_ctrlPts.begin(), m_ctrlPts.end(), cv_right.begin());
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

	bool DegreeElevation(std::vector<Point>& nweCtrlPts)
	{
		int n = m_degree + 1;
		nweCtrlPts.resize(m_degree+2);

		nweCtrlPts.front() = m_ctrlPts.front();
		nweCtrlPts.back() = m_ctrlPts.back();

		for (int i = 1; i < n; ++i)
		{
			double ratio = ((double)i) / n;
			nweCtrlPts[i] = ratio * m_ctrlPts[i - 1] + (1.0 - ratio) * m_ctrlPts[i];
		}

		return true;
	}

private:
    int m_degree;
    std::vector<Point> m_ctrlPts;
};

template <class T>
class BezierCurve_GPT {
public:
	using Point = T;
	BezierCurve_GPT(const std::vector<Point>& controlPoints) : controlPoints(controlPoints) {}

	// Calculate the point on the Bezier curve at parameter t
	Point calculate(double t) const {
		int n = controlPoints.size() - 1;
		Point result = { 0.0, 0.0, 0.0 };

		for (int i = 0; i <= n; ++i) {
			double coefficient = binomialCoefficient(n, i) * pow(1 - t, n - i) * pow(t, i);
			result = result + coefficient * controlPoints[i];
		}

		return result;
	}

	// Calculate the first derivative of the Bezier curve at parameter t
	Point firstDerivative(double t) const {
		int n = controlPoints.size() - 1;
		Point result = { 0.0, 0.0 , 0.0};

		for (int i = 0; i < n; ++i) {
			double coefficient = n * binomialCoefficient(n - 1, i) * pow(1 - t, (n - 1) - i) * pow(t, i);
			result=  result +  coefficient * (controlPoints[i + 1] + (-1.0)* controlPoints[i]);
		}

		return result;
	}

	// Calculate the second derivative of the Bezier curve at parameter t
	Point secondDerivative(double t) const {
		int n = controlPoints.size() - 1;
		Point result = { 0.0, 0.0 , 0.0};

		for (int i = 0; i < n - 1; ++i) {
			double coefficient = n * (n - 1) * binomialCoefficient(n - 2, i) * pow(1 - t, (n - 2) - i) * pow(t, i);
			result = result + coefficient * (controlPoints[i + 2] + (-2.0) * controlPoints[i + 1] + controlPoints[i]);
		}

		return result;
	}

private:
	std::vector<Point> controlPoints;

	// Calculate the binomial coefficient "n choose k"
	int binomialCoefficient(int n, int k) const {
		if (k > n) return 0;
		if (k == 0 || k == n) return 1;
		k = std::min(k, n - k); // take advantage of symmetry
		int c = 1;
		for (int i = 0; i < k; ++i) {
			c = c * (n - i) / (i + 1);
		}
		return c;
	}
};

#endif
