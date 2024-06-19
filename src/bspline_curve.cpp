//#include "bspline_curve.h"
//
//
//std::vector<double> basis_funcs(int span, int degree, const std::vector<double>& knots, double t)
//{
//	std::vector<double> N(degree + 1);
//	N[0] = 1.0;
//
//	std::vector<double> left(degree + 1);
//	std::vector<double> right(degree + 1);
//
//	for (int j = 1; j <= degree; j++)
//	{
//		left[j] = t - m_knots[span + 1 - j];
//		right[j] = m_knots[span + j] - t;
//
//		double saved = 0.0;
//
//		for (int r = 0; r < j; r++)
//		{
//			double temp = N[r] / (right[r + 1] + left[j - r]);
//			N[r] = saved + right[r + 1] * temp;
//			saved = left[j - r] * temp;
//		}
//		N[j] = saved;
//	}
//	return N;
//}
//
////ref: https://en.wikipedia.org/wiki/De_Boor%27s_algorithm
//def deBoor(k: int, x : int, t, c, p : int) :
//	"""Evaluates S(x).
//
//	Arguments
//	-------- -
//	k : Index of knot interval that contains x.
//	x : Position.
//	t : Array of knot positions, needs to be padded as described above.
//	c : Array of control points.
//	p : Degree of B - spline.
//	"""
//	d = [c[j + k - p] for j in range(0, p + 1)]
//
//	for r in range(1, p + 1) :
//		for j in range(p, r - 1, -1) :
//			alpha = (x - t[j + k - p]) / (t[j + 1 + k - r] - t[j + k - p])
//			d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j]
//
//			return d[p]