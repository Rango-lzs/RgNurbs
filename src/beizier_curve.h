#ifndef BEIZIER_CURVE_HH
#define BEIZIER_CURVE_HH

#include "point.h"
#include <vector>

class BeizierCurve
{
public:
	//beizier curve was total defined by control points
	BeizierCurve(int degree, std::vector<Point3d> ctrlPts);

    //
    Point3d EvalPoint(double param);

    Point3d EvalPointDirect(double param);

    Point3d EvalPointByDeCasteljau(double param)
    {
        std::vector<Point3d> temp =m_ctrlPts;
        for (int k = 1; k <= m_degree; k++)
        {
            for (int i = 0; i <= m_degree - k; i++)
            {
                temp[i] = (1.0 - param) * temp[i] + param * temp[i + 1];
            }
        }
        return temp[0];
    }

private:
    int m_degree;
    std::vector<Point3d> m_ctrlPts;
};

#endif
