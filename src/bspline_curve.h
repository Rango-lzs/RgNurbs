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
	using PointArray = std::vector<Point>;

public:
	BSplineCurve(int degree, const std::vector<double>& knots, const std::vector<Point>& ctrlPts)
		:m_degree(degree)
		, m_knots(knots)
		, m_ctrl_pts(ctrlPts)
	{

	}

	//通过基函数计算点, 对于t:[ui,ui+1),仅N(k,p) k = i, i-1, i-2, ... i-p 
	//共p+1个基函数为非零的
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
	//类似Beizier的deCas算法，不停的进行cut cornor操作，deCas方法可以看做是DeBoor方法的特例
	//Deboor方法中切割比例在每个区间是变化的，deCas不变，这是由于Beizie的节点向量是重复的，计算出来的
	//切割系数是一样的
	Point deBoor(double u)
	{
		int span = BSplineBasis::FindSpan(m_knots, m_degree, u);
		int p = m_degree;
		std::vector<Point> d(p + 1);
		for (int i = 0; i < p + 1; i++)
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

	//通过deBoor方法计算点
	Point EvalPointByDeBoor(double t)
	{
		return Point(0, 0, 0);
	}

	//求一阶，二阶导数
	PointArray CalDerivative(int order, double u)
	{
		std::vector<Point> ptRet(order + 1);
		int span = BSplineBasis::FindSpan(m_knots, m_degree, u);
		std::vector<std::vector<double>> dbdu = BSplineBasis::BasisFunctionsDerivatives(span, m_degree, order, m_knots, u);

		order = std::min(order, m_degree);
		for (int i = 0; i <= order; ++i)
		{
			Point du;
			int p = m_degree;
			for (int j = 0; j < p + 1; ++j)
			{
				du = du + m_ctrl_pts[span - p + j] * dbdu[i][j];
			}
			ptRet[i] = du;
		}

		return ptRet;
	}

	// 节点插入算法
	PointArray KnotInsertion(double u, int num = 1) {

		int degree = m_degree;
		int s = BSplineBasis::FindMultiplicity(u, m_knots);
		int k = BSplineBasis::FindSpan(m_knots, degree, u);

		// 初始化变量
		int np = m_ctrl_pts.size();
		int nq = np + num;
		PointArray ctrlpts_new(nq);

		// 初始化一个长度为degree + 1的局部数组
		PointArray temp(degree + 1);

		// 保存未更改的控制点
		for (int i = 0; i < k - degree + 1; ++i) {
			ctrlpts_new[i] = m_ctrl_pts[i];
		}
		for (int i = k - s; i < np; ++i) {
			ctrlpts_new[i + num] = m_ctrl_pts[i];
		}

		// 开始填充将用于在节点插入期间更新控制点的临时局部数组
		for (int i = 0; i <= degree - s; ++i) {
			temp[i] = m_ctrl_pts[k - degree + i];
		}

		// 插入节点 "num" 次
		for (int j = 1; j <= num; ++j) {
			int L = k - degree + j;
			for (int i = 0; i <= degree - j - s; ++i) {
				double alpha = BSplineBasis::knot_insertion_alpha(u, m_knots, k, i, L);
				//double alpha = (u - m_knots[L + i]) / (m_knots[i + k + 1] - m_knots[L + i]);
				temp[i] = alpha * temp[i + 1] + (1.0 - alpha) * temp[i];

			}
			ctrlpts_new[L] = temp[0];
			ctrlpts_new[k + num - j - s] = temp[degree - j - s];
		}

		// 加载剩余的控制点
		int L = k - degree + num;
		for (int i = L + 1; i < k - s; ++i) {
			ctrlpts_new[i] = temp[i - L];
		}

		// 返回插入节点后的控制点
		return ctrlpts_new;
	}

	BSplineCurve<Point> KnotRemoval(double u, int times)
	{
		double tol = 1.0e-6;  //Refer to Eq 5.30 for the meaning
		int s = BSplineBasis::FindMultiplicity(u, m_knots);
		int r = BSplineBasis::FindSpan(m_knots, m_degree, u);
		int ord = m_degree + 1;
		int degree = m_degree;

		//Edge case
		if (times < 1)
		{
			return *this;
		}

		//Initialize variables
		int first = r - degree;
		int last = r - s;

		//Don't change input variables, prepare new ones for updating
		PointArray ctrlpts_new = m_ctrl_pts;
		PointArray temp(2 * degree + 1);

		int t = 0;
		for (t = 0; t < times; t++)
		{
			int off = first - 1;
			temp[0] = m_ctrl_pts[off];
			temp[last + 1 - off] = m_ctrl_pts[last + 1];
			int i = first;
			int j = last;
			int ii = 1;
			int jj = last - off;
			while (j - i >= t)  // j -i >= ?? need to confirm
			{
				double alfi = (u - m_knots[i]) / (m_knots[i + ord + t] - m_knots[i]);
				double alfj = (u - m_knots[j - t]) / (m_knots[j + ord] - m_knots[j - t]);
				temp[ii] = (m_ctrl_pts[i] - (1.0 - alfi) * temp[ii - 1]) / alfi;
				temp[jj] = (m_ctrl_pts[j] - alfj * temp[jj + 1]) / (1.0 - alfj);
				++i; 
				++ii;
				--j; 
				--jj;
			}

			/*if (j - i < t)
			{
				if (MathUtils::IsLessThanOrEqual(temp[ii - 1].Distance(temp[jj + 1]), tol))
				{
					remflag = true;
				}
			}
			else
			{
				double alphai = (removeKnot - knotVector[i]) / (knotVector[i + order + t] - knotVector[i]);
				if (MathUtils::IsLessThanOrEqual(controlPoints[i].Distance(alphai * temp[ii + t + 1] + (1.0 - alphai) * temp[ii - 1]), tol))
				{
					remflag = true;
				}
			}

			if (!remflag)
			{
				break;
			}*/


			i = first; 
			j = last;
			while (j - i > t)
			{
				ctrlpts_new[i] = temp[i - off];
				ctrlpts_new[j] = temp[j - off];
				++i;
				--j;
			}

			--first;
			++last;
		}

		if (t == 0)
		{
			throw;
		}


		int j = (2 * r - s - degree) / 2;
		int i = j;

		for (int k = 1; k < t; k++)
		{
			if (k % 2 == 1)
			{
				i = i + 1;
			}
			else
			{
				j = j - 1;
			}
		}

		int n = m_ctrl_pts.size() - 1;
		for (int k = i + 1; k <= n; k++)
		{
			ctrlpts_new[j] = m_ctrl_pts[k];
			j = j + 1;
		}
		for (int i = 0; i < t; i++)
		{
			ctrlpts_new.pop_back();
		}

		//shift the knots
		std::vector<double> knots_new = m_knots;
		int m = n + degree + 1;
		for (int k = r + 1; k <= m; k++)
		{
			knots_new[k - times] = knots_new[k];
		}
		knots_new.resize(m + 1 - times);

		return BSplineCurve<Point>(m_degree, knots_new, ctrlpts_new);
	}

private:
	int m_degree;
	std::vector<double> m_knots;	//size m
	PointArray m_ctrl_pts;  //size n    m = n + m_degree +1;
};

#endif
