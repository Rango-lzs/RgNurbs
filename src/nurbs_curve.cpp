#include "nurbs_curve.h"
#include "bspline_curve.h"


#define Degree m_dataRep->m_degree
#define CtrlPts m_dataRep->m_ctrl_pts
#define Knots m_dataRep->m_knots
#define Weights m_dataRep->m_weights

#define FetchData(degree,knots,ctrlPts,weights) \
	degree = Degree;\
	knots = Knots;\
	ctrlPts = CtrlPts;\
	weights = Weights;

template<class HPoint>
struct NurbsCurve<HPoint>::DataRep
{
public:
	DataRep(int degree, const std::vector<double>& knots, const std::vector<Point>& ctrlPts, std::vector<double>& weights)
		:m_degree(degree)
		, m_knots(knots)
		, m_ctrl_pts(ctrlPts)
		,m_weights(weights)
	{

	}
	int m_degree;
	std::vector<double> m_knots;	//size m
	HPointArray m_ctrl_pts;  //size n    m = n + Degree +1;
	std::vector<double> m_weights;
};

template<class HPoint>
NurbsCurve<HPoint>::NurbsCurve()
	:m_dataRep(std::make_unique<DataRep>(0, std::vector<double>(), HPointArray(), std::vector<double>()))
{

}

template<class HPoint>
NurbsCurve<HPoint>::NurbsCurve(int degree, const std::vector<double>& knots, const std::vector<HPoint>& ctrlPts, std::vector<double>& weights)
	:m_dataRep(std::make_unique<DataRep>(degree, knots, ctrlPts, weights))
{

}

template<class Point>
NurbsCurve<Point>::~NurbsCurve()
{

}

template<class Point>
NurbsCurve<Point>::NurbsCurve(const NurbsCurve& other)
{
	m_dataRep.release();
	m_dataRep = std::make_unique<DataRep>(other.Degree, other.Knots, other.CtrlPts,other.Weights);
}

template<class Point>
NurbsCurve<Point>& NurbsCurve<Point>::operator=(const NurbsCurve& other)
{
	m_dataRep.release();
	m_dataRep = std::make_unique<DataRep>(other.Degree, other.Knots, other.CtrlPts,other.Weights);
	return *this;
}

template<class HPoint>
typename HPoint::Point NurbsCurve<HPoint>::EvalPoint(double t)
{
	BSplineCurve<HPoint> bspline(degree, knots, ctrlPts);
	HPoint pt = bspline.EvalPointByBasis(t);
	return pt.ToPoint();
}