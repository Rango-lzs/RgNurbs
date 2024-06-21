#ifndef BSPLINE_BASIS_HH
#define BSPLINE_BASIS_HH

#include "vector"

//the class that provide some algorithm of b-spline basis
//所有的计算都是基于节点向量计算
class BSplineBasis
{
public:
	//二分查找计算
	//注意在区间[u(p),u(m-p)]都是最多要算p+1个基函数， 如果在i<p 或者i>m-p 则需要计算的基函数少于p+1个 
	//所以这个时候直接计算开头的p+1个基函数 和 结尾的 p+1个基函数 一定是对的
	//所以这个函数叫findspan不是很准确，而是找到需要计算的最后一个基函数的索引
	static int FindSpan(const std::vector<double>& knots, int p , double u)
	{
		int m = knots.size() - 1;
		int n = m - p - 1;
		if (u == knots[n + 1]) //should deal with floating equal？
		{
			return n;
		}

		int low = p;
		int high = n + 1;
		int mid = (low + high) / 2;
		while (u < knots[mid] || u >= knots[mid + 1])
		{
			if (u < knots[mid])
			{
				high = mid;
			}
			else
			{
				low = mid;
			}
			mid = (low + high) / 2;
		}

		return mid;
	}


	static std::vector<double> BasisFunctions(int span, int degree, const std::vector<double>& knots, double u)
	{
		/*VALIDATE_ARGUMENT(span >= 0, "spanIndex", "SpanIndex must greater than or equals zero.");
		VALIDATE_ARGUMENT(degree >= 0, "degree", "Degree must greater than or equals zero.");
		VALIDATE_ARGUMENT(knots.size() > 0, "knotVector", "KnotVector size must greater than zero.");
		VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knots), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
		VALIDATE_ARGUMENT_RANGE(u, knots[0], knots[knots.size() - 1]);*/

		std::vector<double> basisFunctions(degree + 1);
		basisFunctions[0] = 1.0;

		std::vector<double> left(degree + 1);
		std::vector<double> right(degree + 1);

		for (int j = 1; j <= degree; j++)
		{
			left[j] = u - knots[span + 1 - j];
			right[j] = knots[span + j] - u;

			double saved = 0.0;

			for (int r = 0; r < j; r++)
			{
				double temp = basisFunctions[r] / (right[r + 1] + left[j - r]);
				basisFunctions[r] = saved + right[r + 1] * temp;
				saved = left[j - r] * temp;
			}
			basisFunctions[j] = saved;
		}
		return basisFunctions;
	}

	//计算0到order阶的, 非零基函数导数
	static std::vector<std::vector<double>> BasisFunctionsDerivatives(int span, int degree, int order, const std::vector<double>& knots, double u)
	{
		/*VALIDATE_ARGUMENT(span >= 0, "spanIndex", "SpanIndex must greater than or equals zero.");
		VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
		VALIDATE_ARGUMENT(derivative <= degree, "derivative", "Derivative must not greater than degree.");
		VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
		VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
		VALIDATE_ARGUMENT_RANGE(paramT, knotVector[0], knotVector[knotVector.size() - 1]);*/

		std::vector<std::vector<double>> derivatives(order + 1, std::vector<double>(degree + 1));
		std::vector<std::vector<double>> ndu(degree + 1, std::vector<double>(degree + 1));

		ndu[0][0] = 1.0;

		std::vector<double> left(degree + 1);
		std::vector<double> right(degree + 1);

		double saved = 0.0;
		double temp = 0.0;

		for (int j = 1; j <= degree; j++)
		{
			left[j] = u - knots[span + 1 - j];
			right[j] = knots[span + j] - u;

			saved = 0.0;
			for (int r = 0; r < j; r++)
			{
				ndu[j][r] = right[r + 1] + left[j - r];
				temp = ndu[r][j - 1] / ndu[j][r];

				ndu[r][j] = saved + right[r + 1] * temp;
				saved = left[j - r] * temp;
			}
			ndu[j][j] = saved;
		}

		for (int j = 0; j <= degree; j++)
		{
			derivatives[0][j] = ndu[j][degree];
		}

		std::vector<std::vector<double>> a(2, std::vector<double>(degree + 1));
		for (int r = 0; r <= degree; r++)
		{
			int s1 = 0;
			int s2 = 1;
			a[0][0] = 1.0;

			for (int k = 1; k <= order; k++)
			{
				double d = 0.0;
				int rk = r - k;
				int pk = degree - k;

				if (r >= k)
				{
					a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
					d = a[s2][0] * ndu[rk][pk];
				}

				int j1 = 0;
				int j2 = 0;

				if (rk >= -1)
					j1 = 1;
				else
					j1 = -rk;

				if (r - 1 <= pk)
					j2 = k - 1;
				else
					j2 = degree - r;

				for (int j = j1; j <= j2; j++)
				{
					a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
					d += a[s2][j] * ndu[rk + j][pk];
				}
				if (r <= pk)
				{
					a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
					d += a[s2][k] * ndu[r][pk];
				}
				derivatives[k][r] = d;

				int temp = s1;
				s1 = s2;
				s2 = temp;
			}
		}

		int r = degree;
		for (int k = 1; k <= order; k++)
		{
			for (int j = 0; j <= degree; j++)
			{
				derivatives[k][j] *= r;
			}
			r *= degree - k;
		}
		return derivatives;
	}
};

#endif