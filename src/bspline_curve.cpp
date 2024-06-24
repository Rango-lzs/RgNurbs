#include "bspline_curve.h"
#include "point.h"

//template<class Point>
//using PointArray = BSplineCurve<Point>::PointArray;
#define Degree m_dataRep->m_degree
#define CtrlPts m_dataRep->m_ctrl_pts
#define Knots m_dataRep->m_knots

namespace
{
	int Bin(int n, int r)
	{
		if (r == 0 || r == n)
			return 1;
		else
			return Bin(n - 1, r - 1) + Bin(n - 1, r);
	}
}

template<class Point>
struct BSplineCurve<Point>::DataRep
{
public:
	DataRep(int degree, const std::vector<double>& knots, const std::vector<Point>& ctrlPts)
		:m_degree(degree)
		, m_knots(knots)
		, m_ctrl_pts(ctrlPts)
	{

	}
	int m_degree;
	std::vector<double> m_knots;	//size m
	PointArray m_ctrl_pts;  //size n    m = n + Degree +1;
};

template<class Point>
BSplineCurve<Point>::BSplineCurve()
	:m_dataRep(std::make_unique<DataRep>(0, std::vector<double>(), PointArray()))
{

}

template<class Point>
BSplineCurve<Point>::BSplineCurve(int degree, const std::vector<double>& knots, const std::vector<Point>& ctrlPts)
	:m_dataRep(std::make_unique<DataRep>(degree, knots, ctrlPts))
{

}

template<class Point>
BSplineCurve<Point>::~BSplineCurve()
{

}

template<class Point>
BSplineCurve<Point>::BSplineCurve(const BSplineCurve& other)
{
	m_dataRep.release();
	m_dataRep = std::make_unique<DataRep>(other.Degree, other.Knots, other.CtrlPts);
}

template<class Point>
Point BSplineCurve<Point>::EvalPointByBasis(double t)
{
	Point pt;
	int span = BSplineBasis::FindSpan(Knots, Degree, t);
	std::vector<double> N = BSplineBasis::BasisFunctions(span, Degree, Knots, t);

	for (int i = 0; i <= Degree; i++)
	{
		pt = pt + N[i] * CtrlPts[span - Degree + i];
	}
	return pt;
}

template<class Point>
Point BSplineCurve<Point>::EvalPointByDeBoor(double u)
{
	int span = BSplineBasis::FindSpan(Knots, Degree, u);
	int p = Degree;
	std::vector<Point> d(p + 1);
	for (int i = 0; i < p + 1; i++)
	{
		d[i] = CtrlPts[i + span - p];
	}

	for (int r = 1; r < p + 1; ++r)
	{
		for (int j = p; j > r - 1; --j)
		{
			double alpha = (u - Knots[j + span - p]) / (Knots[j + 1 + span - r] - Knots[j + span - p]);
			d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j];
		}
	}

	return d[p];
}

template<class Point>
typename BSplineCurve<Point>::PointArray BSplineCurve<Point>::CalDerivative(int order, double u)
{
	std::vector<Point> ptRet(order + 1);
	int span = BSplineBasis::FindSpan(Knots, Degree, u);
	std::vector<std::vector<double>> dbdu = BSplineBasis::BasisFunctionsDerivatives(span, Degree, order, Knots, u);

	order = std::min(order, Degree);
	for (int i = 0; i <= order; ++i)
	{
		Point du;
		int p = Degree;
		for (int j = 0; j < p + 1; ++j)
		{
			du = du + CtrlPts[span - p + j] * dbdu[i][j];
		}
		ptRet[i] = du;
	}

	return ptRet;
}

template<class Point>
typename BSplineCurve<Point>::PointArray BSplineCurve<Point>::KnotInsertion(double u, int num = 1) {

	int degree = Degree;
	int s = BSplineBasis::FindMultiplicity(u, Knots);
	int k = BSplineBasis::FindSpan(Knots, degree, u);

	// 初始化变量
	int np = CtrlPts.size();
	int nq = np + num;
	PointArray ctrlpts_new(nq);

	// 初始化一个长度为degree + 1的局部数组
	PointArray temp(degree + 1);

	// 保存未更改的控制点
	for (int i = 0; i < k - degree + 1; ++i) {
		ctrlpts_new[i] = CtrlPts[i];
	}
	for (int i = k - s; i < np; ++i) {
		ctrlpts_new[i + num] = CtrlPts[i];
	}

	// 开始填充将用于在节点插入期间更新控制点的临时局部数组
	for (int i = 0; i <= degree - s; ++i) {
		temp[i] = CtrlPts[k - degree + i];
	}

	// 插入节点 "num" 次
	for (int j = 1; j <= num; ++j) {
		int L = k - degree + j;
		for (int i = 0; i <= degree - j - s; ++i) {
			double alpha = BSplineBasis::knot_insertion_alpha(u, Knots, k, i, L);
			//double alpha = (u - Knots[L + i]) / (Knots[i + k + 1] - Knots[L + i]);
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

template<class Point>
BSplineCurve<Point> BSplineCurve<Point>::KnotRemoval(double u, int times)
{
	double tol = 1.0e-6;  //Refer to Eq 5.30 for the meaning
	int s = BSplineBasis::FindMultiplicity(u, Knots);
	int r = BSplineBasis::FindSpan(Knots, Degree, u);
	int order = Degree + 1;
	int degree = Degree;

	//Edge case
	if (times < 1)
	{
		return *this;
	}

	//Initialize variables
	int first = r - degree;
	int last = r - s;

	//Don't change input variables, prepare new ones for updating
	PointArray ctrlpts_new = CtrlPts;
	PointArray temp(2 * degree + 1);

	int t = 0;
	for (t = 0; t < times; t++)
	{
		int off = first - 1;
		temp[0] = CtrlPts[off];
		temp[last + 1 - off] = CtrlPts[last + 1];
		int i = first;
		int j = last;
		int ii = 1;
		int jj = last - off;
		bool bRemable = false;

		while (j - i >= t)  // j -i >= ?? need to confirm
		{
			double alfi = (u - Knots[i]) / (Knots[i + order + t] - Knots[i]);
			double alfj = (u - Knots[j - t]) / (Knots[j + order] - Knots[j - t]);
			temp[ii] = (CtrlPts[i] - (1.0 - alfi) * temp[ii - 1]) / alfi;
			temp[jj] = (CtrlPts[j] - alfj * temp[jj + 1]) / (1.0 - alfj);
			++i;
			++ii;
			--j;
			--jj;
		}

		if (j - i < t)
		{
			if (temp[ii - 1].DistanceTo(temp[jj + 1]) <= tol)
			{
				bRemable = true;
			}
		}
		else
		{
			double alphai = (u - Knots[i]) / (Knots[i + order + t] - Knots[i]);
			if (CtrlPts[i].DistanceTo(alphai * temp[ii + t + 1] + (1.0 - alphai) * temp[ii - 1]) <= tol)
			{
				bRemable = true;
			}
		}

		if (!bRemable)
		{
			//throw;
		}


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

	int n = CtrlPts.size() - 1;
	for (int k = i + 1; k <= n; k++)
	{
		ctrlpts_new[j] = CtrlPts[k];
		j = j + 1;
	}
	for (int i = 0; i < t; i++)
	{
		ctrlpts_new.pop_back();
	}

	//shift the knots
	std::vector<double> knots_new = Knots;
	int m = n + degree + 1;
	for (int k = r + 1; k <= m; k++)
	{
		knots_new[k - times] = knots_new[k];
	}
	knots_new.resize(m + 1 - times);

	return BSplineCurve<Point>(Degree, knots_new, ctrlpts_new);
}


template<class Point>
BSplineCurve<Point> BSplineCurve<Point>::KnotRefine(std::vector<double>& us)
{
	int degree = Degree;
	const std::vector<double>& knots = Knots;
	const PointArray& ctrlPts = CtrlPts;
	double tol = 1.0e-4;

	int n = ctrlPts.size() - 1;
	int m = n + degree + 1;
	int r = us.size() - 1;

	int a = BSplineBasis::FindSpan(knots, degree, us[0]);
	int b = BSplineBasis::FindSpan(knots, degree, us[r]) + 1;

	std::vector<double> knots_new(m + r + 2);
	for (int j = 0; j <= a; j++)
	{
		knots_new[j] = knots[j];
	}
	for (int j = b + degree; j <= m; j++)
	{
		knots_new[j + r + 1] = knots[j];
	}

	PointArray ctrlPts_new(n + r + 2);
	for (int j = 0; j <= a - degree; j++)
	{
		ctrlPts_new[j] = ctrlPts[j];
	}
	for (int j = b - 1; j <= n; j++)
	{
		ctrlPts_new[j + r + 1] = ctrlPts[j];
	}

	int i = b + degree - 1;
	int k = b + degree + r;
	for (int j = r; j >= 0; j--)
	{
		while (knots_new[j] <= knots[i] && i > a)
		{
			ctrlPts_new[k - degree - 1] = ctrlPts[i - degree - 1];
			knots_new[k] = knots[i];
			k = k - 1;
			i = i - 1;
		}

		ctrlPts_new[k - degree - 1] = ctrlPts_new[k - degree];
		for (int l = 1; l <= degree; l++)
		{
			int ind = k - degree + l;
			double alpha = knots_new[k + l] - us[j];
			if (abs(alpha) <= tol)
			{
				ctrlPts_new[ind - 1] = ctrlPts_new[ind];
			}
			else
			{
				alpha = alpha / (knots_new[k + l] - knots[i - degree + l]);
				ctrlPts_new[ind - 1] = alpha * ctrlPts_new[ind - 1] + (1.0 - alpha) * ctrlPts_new[ind];
			}
		}
		knots_new[k] = us[j];
		k = k - 1;
	}
	return BSplineCurve<Point>(Degree, knots_new, ctrlPts_new);
}

template <class Point>
void BSplineCurve<Point>::Decompose(std::vector<BSplineCurve<Point>>& bezierSegs) const
{
	int degree = Degree;
	int i, m, a, b, nb, mult, j, r, save, s, k;
	double numer, alpha;
	std::vector<double> alphas(degree + 1);

	// all the curves will have the same vector
	std::vector<double> bezier_knot;
	bezier_knot.resize(2 * (degree + 1));
	for (i = 0; i < bezier_knot.size() / 2; i++)
		bezier_knot[i] = 0;
	for (i = bezier_knot.size() / 2; i < bezier_knot.size(); i++)
		bezier_knot[i] = 1;

	bezierSegs.resize(CtrlPts.size() - degree);
	for (i = 0; i < bezierSegs.size(); i++) {
		bezierSegs[i].Knots = bezier_knot;
		bezierSegs[i].Degree = degree;
		bezierSegs[i].CtrlPts.resize(degree +1);
	}

	auto ctrlPts = CtrlPts;
	auto knots = Knots;
	m = ctrlPts.size() + degree;
	a = degree;
	b = degree + 1;
	nb = 0;

	for (i = 0; i <= degree; i++)
		bezierSegs[nb].CtrlPts[i] = ctrlPts[i];
	while (b < m) {
		i = b;
		while (b < m && knots[b + 1] <= knots[b]) b++;
		mult = b - i + 1;
		if (mult < degree) {
			numer = knots[b] - knots[a]; // th enumerator of the alphas
			for (j = degree; j > mult; j--) // compute and store the alphas
				alphas[j - mult - 1] = numer / (knots[a + j] - knots[a]);
			r = degree - mult; // insert knot r times
			for (j = 1; j <= r; j++) {
				save = r - j;
				s = mult + j; // this many new points
				for (k = degree; k >= s; k--) {
					alpha = alphas[k - s];
					bezierSegs[nb].CtrlPts[k] = alpha * bezierSegs[nb].CtrlPts[k] + (1.0 - alpha) * bezierSegs[nb].CtrlPts[k - 1];
				}
				if (b < m) // control point of
					bezierSegs[nb + 1].CtrlPts[save] = bezierSegs[nb].CtrlPts[degree]; // next segment
			}
		}
		nb++;
		if (b < m) { // initialize for next segment
			for (i = degree - mult; i <= degree; i++)
				bezierSegs[nb].CtrlPts[i] = ctrlPts[b - degree + i];
			a = b;
			b++;
		}
	}
	bezierSegs.resize(nb);
}

template <class Point>
BSplineCurve<Point> BSplineCurve<Point>::degreeElevate(int t) {
	if (t <= 0) {
		throw;
	}
	//NurbsCurve<T, N> c(*this);

	int i, j, k;
	int n = CtrlPts.size() - 1;
	int p = Degree;
	int m = n + p + 1;
	int ph = p + t;
	int ph2 = ph / 2;

	std::vector<std::vector<double>> bezalfs(p + t + 1, std::vector<double>(p + 1));
	//Matrix<T> bezalfs(p + t + 1, p + 1); // coefficients for degree elevating the Bezier segment
	PointArray bpts(p + 1); // pth-degree Bezier control points of the current segment
	PointArray ebpts(p + t + 1); // (p+t)th-degree Bezier control points of the  current segment
	PointArray Nextbpts(p - 1); // leftmost control points of the next Bezier segment
	std::vector<double> alphas(p - 1); // knot instertion alphas.

	// Compute the binomial coefficients
	//std::vector<std::vector<double>> Bin(ph + 1, std::vector<double>(ph2 + 1));
	//binomialCoef(Bin);

	// Compute Bezier degree elevation coefficients
	double inv, mpi;
	bezalfs[0][0] = bezalfs[ph][p] = 1.0;
	for (i = 1; i <= ph2; i++) {
		inv = 1.0 / Bin(ph, i);
		mpi = std::min(p, i);
		for (j = std::max(0, i - t); j <= mpi; j++) {
			bezalfs[i][j] = inv * Bin(p, j) * Bin(t, i - j);
		}
	}

	for (i = ph2 + 1; i < ph; i++) {
		mpi = std::min(p, i);
		for (j = std::max(0, i - t); j <= mpi; j++)
			bezalfs[i][j] = bezalfs[ph - i][p - j];
	}

	int mh = ph;
	int kind = ph + 1;
	double ua = Knots[0];
	double ub = 0.0;
	int r = -1;
	int oldr;
	int a = p;
	int b = p + 1;
	int cind = 1;
	int rbz, lbz = 1;
	int mul, save, s;
	double alf;
	int first, last, kj;
	double den, bet, gam, numer;

	//resize(c.P.n() + c.P.n() * t, ph); // Allocate more control points than necessary
	int size_pre = CtrlPts.size() + CtrlPts.size() * t;
	PointArray ctrlPts_new(size_pre);
	ctrlPts_new[0] = CtrlPts[0];

	std::vector<double> knots_new(size_pre + ph + 1);
	for (int i = 0; i <= ph; i++)
	{
		knots_new[i] = ua;
	}

	for (int i = 0; i <= p; i++)
	{
		bpts[i] = CtrlPts[i];
	}

	ctrlPts_new[0] = CtrlPts[0];
	for (i = 0; i <= ph; i++) {
		knots_new[i] = ua;
	}

	// Initialize the first Bezier segment

	for (i = 0; i <= p; i++)
		bpts[i] = CtrlPts[i];

	while (b < m) { // Big loop thru knot vector
		i = b;
		while (b < m && Knots[b] >= Knots[b + 1]) // for some odd reasons... == doesn't work
			b++;
		mul = b - i + 1;
		mh += mul + t;
		ub = Knots[b];
		oldr = r;
		r = p - mul;
		if (oldr > 0)
			lbz = (oldr + 2) / 2;
		else
			lbz = 1;
		if (r > 0)
			rbz = ph - (r + 1) / 2;
		else
			rbz = ph;
		if (r > 0) { // Insert knot to get Bezier segment
			numer = ub - ua;
			for (k = p; k > mul; k--) {
				alphas[k - mul - 1] = numer / (Knots[a + k] - ua);
			}
			for (j = 1; j <= r; j++) {
				save = r - j; s = mul + j;
				for (k = p; k >= s; k--) {
					bpts[k] = alphas[k - s] * bpts[k] + (1.0 - alphas[k - s]) * bpts[k - 1];
				}
				Nextbpts[save] = bpts[p];
			}
		}

		for (i = lbz; i <= ph; i++) { // Degree elevate Bezier,  only the points lbz,...,ph are used
			ebpts[i] = Point();
			mpi = std::min(p, i);
			for (j = std::max(0, i - t); j <= mpi; j++)
				ebpts[i] += bezalfs[i][j] * bpts[j];
		}

		if (oldr > 1) { // Must remove knot u=c.U[a] oldr times
		  // if(oldr>2) // Alphas on the right do not change
		  //	alfj = (ua-U[kind-1])/(ub-U[kind-1]) ;
			first = kind - 2; last = kind;
			den = ub - ua;
			bet = (ub - knots_new[kind - 1]) / den;
			for (int tr = 1; tr < oldr; tr++) { // Knot removal loop
				i = first; j = last;
				kj = j - kind + 1;
				while (j - i > tr) { // Loop and compute the new control points for one removal step
					if (i < cind) {
						alf = (ub - knots_new[i]) / (ua - knots_new[i]);
						ctrlPts_new[i] = alf * ctrlPts_new[i] + (1.0 - alf) * ctrlPts_new[i - 1];
					}
					if (j >= lbz) {
						if (j - tr <= kind - ph + oldr) {
							gam = (ub - knots_new[j - tr]) / den;
							ebpts[kj] = gam * ebpts[kj] + (1.0 - gam) * ebpts[kj + 1];
						}
						else {
							ebpts[kj] = bet * ebpts[kj] + (1.0 - bet) * ebpts[kj + 1];
						}
					}
					++i; --j; --kj;
				}
				--first; ++last;
			}
		}

		if (a != p) // load the knot u=c.U[a]
			for (i = 0; i < ph - oldr; i++) {
				knots_new[kind++] = ua;
			}
		for (j = lbz; j <= rbz; j++) { // load control points onto the curve
			ctrlPts_new[cind++] = ebpts[j];
		}

		if (b < m) { // Set up for next pass thru loop
			for (j = 0; j < r; j++)
				bpts[j] = Nextbpts[j];
			for (j = r; j <= p; j++)
				bpts[j] = CtrlPts[b - p + j];
			a = b;
			b++;
			ua = ub;
		}
		else {
			for (i = 0; i <= ph; i++)
				knots_new[kind + i] = ub;
		}
	}

	ctrlPts_new.resize(mh - ph);
	int order_h = ph;
	knots_new.resize(mh + 1);
	return BSplineCurve<Point>(ph, knots_new, ctrlPts_new);
	//resize(mh - ph, ph); // Resize to the proper number of control points
}

template<class Point>
bool BSplineCurve<Point>::BezDegreeReduce(BSplineCurve<Point>& result)
{
	int degree = Degree;
	const std::vector<double>& knots = Knots;
	PointArray& ctrlPts = CtrlPts;

	double tol = 0.01;

	/*int size = controlPoints.size();
	bool isBezier = ValidationUtils::IsValidBezier(degree, size);
	if (!isBezier) return false;*/

	int r = floor((degree - 1) / 2);
	PointArray ctrlPts_new(degree);
	ctrlPts_new[0] = ctrlPts[0];
	ctrlPts_new[degree - 1] = ctrlPts[degree];

	std::vector<double> alpha(degree);
	for (int i = 0; i < degree; i++)
	{
		alpha[i] = double(i) / double(degree);
	}
	double error = 0.0;
	if (degree % 2 == 0)
	{
		for (int i = 1; i < r + 1; i++)
		{
			ctrlPts_new[i] = (ctrlPts[i] - alpha[i] * ctrlPts_new[i - 1]) / (1 - alpha[i]);
		}
		for (int i = degree - 2; i > r; i--)
		{
			ctrlPts_new[i] = (ctrlPts[i + 1] - (1 - alpha[i + 1]) * ctrlPts_new[i + 1]) / alpha[i + 1];
		}
		error = (ctrlPts[r + 1].DistanceTo(0.5 * (ctrlPts_new[r] + ctrlPts_new[r + 1])));
		double c = Bin(degree, r + 1);
		error = error * (c * pow(0.5, r + 1) * pow(1 - 0.5, degree - r - 1));
	}
	else
	{
		for (int i = 1; i < r; i++)
		{
			ctrlPts_new[i] = (ctrlPts[i] - alpha[i] * ctrlPts_new[i - 1]) / (1 - alpha[i]);
		}
		for (int i = degree - 2; i > r; i--)
		{
			ctrlPts_new[i] = (ctrlPts[i + 1] - (1 - alpha[i + 1]) * ctrlPts_new[i + 1]) / alpha[i + 1];
		}
		Point PLr = (ctrlPts[r] - alpha[r] * ctrlPts_new[r - 1]) / (1 - alpha[r]);
		Point PRr = (ctrlPts[r + 1] - (1 - alpha[r + 1]) * ctrlPts_new[r + 1]) / alpha[r + 1];
		ctrlPts_new[r] = 0.5 * (PLr + PRr);
		error = PLr.DistanceTo(PRr);
		double maxU = (degree - std::sqrt(degree)) / (2 * degree);
		error = error * 0.5 * (1 - alpha[r]) * (Bin(degree, r) * pow(maxU, r) * pow(1 - maxU, r + 1) * (1 - 2 * maxU));
	}
	if (error > tol) return false;

	//auto map = KnotVectorUtils::GetKnotMultiplicityMap(knotVector);
	std::vector<double> knots_new(2*degree);
	for (int i = 0; i < degree;i++)
	{
		knots_new[i] = knots.front();
		knots_new[2 * degree - i - 1] = knots.back();
	}
	
	result.Degree = degree - 1;
	result.Knots = knots_new;
	result.CtrlPts = ctrlPts_new;
	return true;
}


template<class Point>
BSplineCurve<Point> BSplineCurve<Point>::DegreeReduce()
{
	int p = Degree;
	int m = Knots.size() - 1;
	PointArray bpts(p + 1);
	PointArray Nextbpts(p - 1);
	PointArray rbpts(p);
	std::vector<double> alphas(p - 1);
	std::vector<double> e(m);

	int ph = p - 1;
	int mh = ph;
	int kind = ph + 1;
	int r = -1;
	int a = p;
	int b = p + 1;
	int cind = 1;
	int mult = p;

	PointArray ctrlPts_new(CtrlPts.size());
	std::vector<double> knots_new(mh + 1);
	ctrlPts_new[0] = CtrlPts[0];
	for (int i = 0; i <= ph; i++)
	{
		knots_new[i] = Knots[0]; // left p+1 knots
	}
	
	for (int i = 0; i < m; i++)
	{
		e[i] = 0.0;
	}

	//loop over the knots
	while (b < m)
	{
		int i = b;
		while (b < m && Knots[b] == Knots[b + 1]) {
			b = b + 1;
		}
		mult = b - i + 1;
		mh = mh + mult - 1;
		int oldr = r;
		r = p - mult;

		int lbz = 0;
		if (oldr > 0) {
			lbz = (oldr + 2) / 2;
		}
		else
		{
			lbz = 1;
		}

		//插入节点，提取beizier
		if (r > 0)
		{
			double numer = Knots[b] - Knots[a];
			for (int k = p; k >= mult; k--)
			{
				alphas[k - mult - 1] = numer / (Knots[a + k] - Knots[a]);
			}

			for (int j = 1; j <= r; j++)
			{
				int save = r - j;
				s = mult + j;
				for (int k = p; k >= s; k--)
				{
					bpts[k] = alphas[k - s] * bpts[k] + (1 - alphas[k - s]) * bpts[k - 1];
				}
				Nextbpts[save] = bpts[0];
			}
		}

		//Beizier 降解
		BezDegreeReduce(bpts, rbpts, MaxErr);
		e[a] = e[a] + MaxErr;
		if (e[a] > Tol)
			return false;

		//消去节点U[a] oldr次

	}

}


template class RG_API BSplineCurve<Point3d>;
template class RG_API BSplineCurve<Point2d>;