

#ifndef RG_MATH_H
#define RG_MATH_H

#include <math.h>

#define RG_ZERO          1.0e-12                          /* System zero */
#define RG_ISZERO(x)      (-RG_ZERO < (x) && (x) < RG_ZERO)
#define RG_DEG_TO_RAD    0.01745329251994329577           /* (1.0/180.0) * pi */
#define RG_INFINITY      1.0e+20                          /* Computer infinity */
#define RG_PI            3.141592653589793238463          /* pi (180 deg in radians) */
#define RG_2PI           6.283185307179586476925          /* 2.0 * pi */
#define RG_PI_INV        0.318309886183790671538          /* 1.0 / pi */
#define RG_PI_DIV_2      1.570796326794896619231          /* pi / 2.0 */
#define RG_PI_DIV_4      0.785398163339744830962          /* pi / 4.0 */
#define RG_PI_DIV_6      0.523598775598298873077          /* pi / 6.0 */
#define RG_SQRT2         1.414213562373095048801          /* sqrt(2) */

namespace RgMath
{
    template<typename T>
    struct Lim1
    {
        T min;
        T max;
    };
    template<typename T>
    inline bool operator==(const Lim1<T>& lim1, const Lim1<T>& lim2) { return lim1.min == lim2.min && lim1.max == lim2.max; }

    template<typename T>
    inline bool operator!=(const Lim1<T>& lim1, const Lim1<T>& lim2) { return !(lim1 == lim2); }

    template<typename T>
    struct Lim2
    {
        Lim1<T> x;
        Lim1<T> y;      
    };
    template<typename T>
    inline bool operator==(const Lim2<T> lim1, const Lim2<T>& lim2) { return lim1.x == lim2.x && lim1.y == lim2.y; }

    template<typename T>
    inline bool operator!=(const Lim2<T>& lim1, const Lim2<T>& lim2) { return !(lim1 == lim2); }

    template<typename T>
    struct Lim3
    {
        Lim1<T> x;
        Lim1<T> y;
        Lim1<T> z;     
    };
    template<typename T>
    inline bool operator==(const Lim3<T>& lim1, const Lim3<T>& lim2) { return lim1.x == lim2.x && lim1.y == lim2.y&& lim1.z == lim2.z; }

    template<typename T>
    inline bool operator!=(const Lim3<T>& lim1, const Lim3<T>& lim2) { return !(lim1 == lim2); }

    template<typename T>
    void Lim2Init(Lim2<T> &lim)
    {
        lim.x.min = (T)RG_INFINITY;
        lim.x.max = -(T)RG_INFINITY;
        lim.y.min = (T)RG_INFINITY;
        lim.y.max = -(T)RG_INFINITY;
    }

    template<typename T>
    void Lim3Init(Lim3<T> &lim)
    {
        lim.x.min = (T)RG_INFINITY;
        lim.x.max = -(T)RG_INFINITY;
        lim.y.min = (T)RG_INFINITY;
        lim.y.max = -(T)RG_INFINITY;
        lim.z.min = (T)RG_INFINITY;
        lim.z.max = -(T)RG_INFINITY;
    }

    template<typename T>
    struct Point2
    {
        constexpr Point2(T _x = 0, T _y = 0) :x(_x), y(_y) {}

        Point2(T p[2])
        {
            x = p[0];
            y = p[1];
        }

        Point2(const T &p)
        {
            x = p.x;
            y = p.y;
        }

        Point2(const T *p)
        {
            x = p[0];
            y = p[1];
        }

        T x;
        T y;

        T& operator[](int idx)
        {
            return *(data() + idx);
        }

        T* data()
        {
            return reinterpret_cast<T*> (this);
        }

        T Length() const
        {
            return sqrt(x * x + y * y);
        }

        T Distance(const Point2<T> &p) const
        {
            return sqrt((x - p.x) * (x - p.x) + (y - p.y) * (y - p.y));
        }

        bool IsSameTol(const Point2<T> &p, T tol = (T)RG_ZERO) const
        {
            if (fabs(x - p.x) > tol)
                return false;
            if (fabs(y - p.y) > tol)
                return false;
            return true;
        }
    };

    template<typename T>
    inline Point2<T> operator +(const Point2<T>& p1, const Point2<T>& p2) noexcept { return { p1.x + p2.x, p1.y + p2.y }; }

    template<typename T>
    inline Point2<T> operator -(const Point2<T>& p1, const Point2<T>& p2) noexcept { return { p1.x - p2.x, p1.y - p2.y }; }

    template<typename T>
    inline Point2<T> operator *(T d, const Point2<T>& p1) noexcept { return { p1.x * d, p1.y * d }; }

    template<typename T>
    inline Point2<T> operator /(const Point2<T>& p1, T d) noexcept { return (1.0 / d) * p1; }

    template<typename T>
    inline T operator *(const Point2<T>& p1, const Point2<T>& p2) noexcept { return { p1.x * p2.y - p1.y * p2.x }; }

    template<typename T>
    inline T operator %(const Point2<T>& p1, const Point2<T>& p2) noexcept { return { p1.x * p2.x + p1.y * p2.y }; }

    template<typename T>
    inline Point2<T> operator -(const Point2<T>& p1) noexcept { return { -p1.x, -p1.y }; }

    template<typename T>
    inline Point2<T> operator +=(Point2<T>& p, const Point2<T>& p1) noexcept { p.x += p1.x; p.y += p1.y; return p; }

    template<typename T>
    inline Point2<T> operator -=(Point2<T>& p, const Point2<T>& p1) noexcept { p.x -= p1.x; p.y -= p1.y; return p; }

    template<typename T>
    inline Point2<T> operator *=(Point2<T>& p, T d) noexcept { p.x *= d; p.y *= d; return p; }

    template<typename T>
    inline Point2<T> operator /=(Point2<T>& p, T d) noexcept { p.x /= d; p.y /= d; return p; }

    template<typename T>
    inline bool operator ==(const Point2<T>& p, const Point2<T>& other) noexcept { return p.x == other.x && p.y == other.y; }

    template<typename T>
    inline bool operator<(const Point2<T>& p1, const Point2<T>& p2) noexcept
    {
        if (p1.x < p2.x || (p1.x == p2.x && p1.y < p2.y))
            return true;
        else
            return false;
    }

    template<typename T>
    struct Point3
    {
        constexpr Point3(T _x = 0, T _y = 0, T _z = 0) : x(_x), y(_y), z(_z) {}

        Point3(T p[3])
        {
            x = p[0];
            y = p[1];
            z = p[2];
        }

        Point3(const T &p)
        {
            x = p.x;
            y = p.y;
            z = p.z;
        }

        Point3(const T *p)
        {
            x = p[0];
            y = p[1];
            z = p[2];
        }

        T x;
        T y;
        T z;

        T& operator[](int idx)
        {
            return *(const_cast<T *>(data()) + idx);
        }

        const T* data() const
        {
            return reinterpret_cast<const T*> (this);
        }

        void Revese()
        {
            x *= -1;
            y *= -1;
            z *= -1;
        }

        T Length() const
        {
            return sqrt(x * x + y * y + z * z);
        }

        bool Normalize()
        {
            // refer to "diV3Normalize"
            constexpr T V_GEOM_ZERO = (T)RG_ZERO;
            T length = Length();
            if (length < V_GEOM_ZERO)
            {
                return false;
            }
            this->x /= length;
            this->y /= length;
            this->z /= length;
            return true;
        }

        T Distance(const Point3 &p) const
        {
            return sqrt((x - p.x) * (x - p.x) + (y - p.y) * (y - p.y) + (z - p.z) * (z - p.z));
        }

        T Distance2(const Point3 &p) const
        {
            return (x - p.x) * (x - p.x) + (y - p.y) * (y - p.y) + (z - p.z) * (z - p.z);
        }

        bool IsSameTol(const Point3 &p, T tol) const
        {
            if (fabs(x - p.x) > tol)
                return false;
            if (fabs(y - p.y) > tol)
                return false;
            if (fabs(z - p.z) > tol)
                return false;
            return true;
        }

        template<typename ToType = T>
        ZwArray<ToType, 3> ToArray() const
        {
            return { (ToType)x, (ToType)y, (ToType)z };
        }
    };

    template<typename T>
    T Dot(const Point3<T> &p1, const Point3<T> &p2)
    {
        return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
    }

    template<typename T>
    Point3<T> Cross(const Point3<T> &p1, const Point3<T> &p2)
    {
        return {
            p1.y * p2.z - p1.z * p2.y,
            p1.z * p2.x - p1.x * p2.z,
            p1.x * p2.y - p1.y * p2.x
        };
    }

    template<typename T>
    T Norm(const Point3<T> &p)
    {
        return sqrt(Dot(p, p));
    }

    template<typename T>
    T NormSq(const Point3<T> &v)
    {
        return Dot(v, v);
    }

    template<typename T>
    T Angle(const Point3<T> &a, const Point3<T> &b)
    {
        T cosTheta = Dot(a, b);
        T sinTheta = Norm(Cross(a, b));
        return (T)atan2(sinTheta, cosTheta);
    }

    template<typename T>
    inline Point3<T> operator +(const Point3<T>& p1, const Point3<T>& p2) noexcept { return { p1.x + p2.x, p1.y + p2.y, p1.z + p2.z }; }

    template<typename T>
    inline Point3<T> operator -(const Point3<T>& p1, const Point3<T>& p2) noexcept { return { p1.x - p2.x, p1.y - p2.y, p1.z - p2.z }; }

    template<typename T>
    inline Point3<T> operator *(T d, const Point3<T>& p1) noexcept { return { d * p1.x, d * p1.y, d * p1.z, }; }

    template<typename T>
    inline Point3<T> operator *(const Point3<T>& p1, T d) noexcept { return { d * p1.x, d * p1.y, d * p1.z, }; }

    template<typename T>
    inline Point3<T> operator *(const Point3<T>& p1, int d) noexcept { return { (T)d * p1.x, (T)d * p1.y, (T)d * p1.z, }; }

    template<typename T>
    inline Point3<T> operator *(const Point3<T>& p1, const Point3<T>& p2) noexcept { return { p1.y * p2.z - p1.z * p2.y, p1.z * p2.x - p1.x * p2.z, p1.x * p2.y - p1.y * p2.x, }; }

    template<typename T>
    inline Point3<T> operator /(const Point3<T>& p1, T d) noexcept { return ((T)1.0 / d) * p1; }

    template<typename T>
    inline Point3<T> operator /(const Point3<T>& p1, int d) noexcept { return ((T)1.0 / (T)d) * p1; }

    template<typename T>
    inline T operator %(const Point3<T>& p1, const Point3<T>& p2) noexcept { return { p1.x * p2.x + p1.y * p2.y + p1.z * p2.z }; }

    template<typename T>
    inline Point3<T> operator -(const Point3<T>& p1) noexcept { return { -p1.x, -p1.y, -p1.z }; }

    template<typename T>
    inline Point3<T> operator +=(Point3<T>& p, const Point3<T>& p1) noexcept { p.x += p1.x; p.y += p1.y; p.z += p1.z; return p; }

    template<typename T>
    inline Point3<T> operator -=(Point3<T>& p, const Point3<T>& p1) noexcept { p.x -= p1.x; p.y -= p1.y; p.z -= p1.z; return p; }

    template<typename T>
    inline Point3<T> operator *=(Point3<T>& p, T d) noexcept { p.x *= d; p.y *= d; p.z *= d; return p; }

    template<typename T>
    inline Point3<T> operator *=(Point3<T>& p, int d) noexcept { p.x *= (T)d; p.y *= (T)d; p.z *= (T)d; return p; }

    template<typename T>
    inline Point3<T> operator /=(Point3<T>& p, T d) noexcept { p.x /= d; p.y /= d; p.z /= d; return p; }

    template<typename T>
    inline Point3<T> operator /=(Point3<T>& p, int d) noexcept { p.x /= (T)d; p.y /= (T)d; p.z /= (T)d; return p; }

    template<typename T>
    inline bool operator ==(const Point3<T>& p1, const Point3<T>& p2)
    {
        return p1.x == p2.x && p1.y == p2.y && p1.z == p2.z;
    }

    template<typename T>
    inline bool operator !=(const Point3<T>& p1, const Point3<T>& p2)
    {
        return p1.x != p2.x || p1.y != p2.y || p1.z != p2.z;
    }

    template<typename T>
    inline bool operator <(const Point3<T>& p1, const Point3<T>& p2) noexcept
    {
        if (p1.x < p2.x || (p1.x == p2.x && p1.y < p2.y) || (p1.x == p2.x && p1.y == p2.y && p1.z < p2.z))
            return true;
        else
            return false;
    }

    template<typename T>
    struct Point5
    {
        T u;
        T v;
        T x;
        T y;
        T z;
    };

    template<typename T>
    struct Point8
    {
        T u;
        T v;
        T x;
        T y;
        T z;
        T i;
        T j;
        T k;
    };

    template<typename T>
    struct Mat3
    {
        char identity;              // 1 if identity matrix, else 0.
        T xx, yx, zx, xt;
        T xy, yy, zy, yt;
        T xz, yz, zz, zt;

        void IdentityMatrix()
        {
            // refer to "VgIdentMat"
            memset(this, 0, sizeof(*this));
            this->identity = 1;
            this->xx = 1.0;
            this->yy = 1.0;
            this->zz = 1.0;
        }

        void UpdateIdentityFlag()
        {
            // refer to "VmUpdateIdent"
            this->identity = 0;
            constexpr T V_GEOM_ZERO = (T)1.0e-12;
            if ((fabs(this->xx - 1.0) < V_GEOM_ZERO) && (fabs(this->yy - 1.0) < V_GEOM_ZERO)
                && (fabs(this->zz - 1.0) < V_GEOM_ZERO) && (fabs(this->xt) < V_GEOM_ZERO)
                && (fabs(this->yt) < V_GEOM_ZERO) && (fabs(this->zt) < V_GEOM_ZERO))
            {
                this->identity = 1;
            }
        }

        template<typename GLType>
        void TranToGLMatrix(GLType *dst)
        {
            // copy from DD_MATRIX_COPY
            // The OpenGL matrix is stored in column major order
            dst[0] = (GLType)this->xx; dst[4] = (GLType)this->yx; dst[8] =  (GLType)this->zx;  dst[12] = (GLType)this->xt;
            dst[1] = (GLType)this->xy; dst[5] = (GLType)this->yy; dst[9] =  (GLType)this->zy;  dst[13] = (GLType)this->yt;
            dst[2] = (GLType)this->xz; dst[6] = (GLType)this->yz; dst[10] = (GLType)this->zz;  dst[14] = (GLType)this->zt;
            dst[3] = (GLType)0.0;      dst[7] = (GLType)0.0;      dst[11] = (GLType)0.0;       dst[15] = (GLType)1.0;
        }

        void XformDir(Point3<T> *dir) const
        {
            // refer to "VxMat3Dir"
            if (this->identity == 0)
            {
                Point3<T> vec;
                vec = *dir;
                dir->x = vec.x * this->xx + vec.y * this->yx + vec.z * this->zx;
                dir->y = vec.x * this->xy + vec.y * this->yy + vec.z * this->zy;
                dir->z = vec.x * this->xz + vec.y * this->yz + vec.z * this->zz;
            }
        }

        void XformInvDir(Point3<T> *invDir) const
        {
            // refer to "VxMat3InvDir"
            if (this->identity == 0)
            {
                Point3<T> vec = *invDir;
                invDir->x = vec.x * this->xx + vec.y * this->xy + vec.z * this->xz;
                invDir->y = vec.x * this->yx + vec.y * this->yy + vec.z * this->yz;
                invDir->z = vec.x * this->zx + vec.y * this->zy + vec.z * this->zz;
            }
        }

        void XformPnt(Point3<T> *pnt) const
        {
            // refer to "VxMat3Pnt"
            if (this->identity == 0)
            {
                Point3<T> pt = *pnt;
                pnt->x = this->xt + pt.x*this->xx + pt.y*this->yx + pt.z*this->zx;
                pnt->y = this->yt + pt.x*this->xy + pt.y*this->yy + pt.z*this->zy;
                pnt->z = this->zt + pt.x*this->xz + pt.y*this->yz + pt.z*this->zz;
            }
        }

        void XformInvPnt(Point3<T> *pnt) const
        {
            // refer to "VxMat3InvPnt"
            if (this->identity == 0)
            {
                T dx, dy, dz;

                dx = pnt->x - this->xt;
                dy = pnt->y - this->yt;
                dz = pnt->z - this->zt;

                pnt->x = this->xx*dx + this->xy*dy + this->xz*dz;
                pnt->y = this->yx*dx + this->yy*dy + this->yz*dz;
                pnt->z = this->zx*dx + this->zy*dy + this->zz*dz;
            }
        }

        // res = this * mul
        int MultMat(const Mat3 *mul, Mat3 *res) const
        {
            // refer to "diMultMat"
            if (this->identity) { *res = *mul; return(0); }
            if (mul->identity) { *res = *this; return(0); }
            Mat3 temp;
            temp.xx = this->xx*mul->xx + this->yx*mul->xy + this->zx*mul->xz;
            temp.xy = this->xy*mul->xx + this->yy*mul->xy + this->zy*mul->xz;
            temp.xz = this->xz*mul->xx + this->yz*mul->xy + this->zz*mul->xz;
            temp.yx = this->xx*mul->yx + this->yx*mul->yy + this->zx*mul->yz;
            temp.yy = this->xy*mul->yx + this->yy*mul->yy + this->zy*mul->yz;
            temp.yz = this->xz*mul->yx + this->yz*mul->yy + this->zz*mul->yz;
            temp.zx = this->xx*mul->zx + this->yx*mul->zy + this->zx*mul->zz;
            temp.zy = this->xy*mul->zx + this->yy*mul->zy + this->zy*mul->zz;
            temp.zz = this->xz*mul->zx + this->yz*mul->zy + this->zz*mul->zz;
            temp.xt = this->xx*mul->xt + this->yx*mul->yt + this->zx*mul->zt + this->xt;
            temp.yt = this->xy*mul->xt + this->yy*mul->yt + this->zy*mul->zt + this->yt;
            temp.zt = this->xz*mul->xt + this->yz*mul->yt + this->zz*mul->zt + this->zt;
            temp.identity = 0;
            *res = temp;
            return(0);
        }

        void Origin(Point3<T>&pt) const
        {
            pt.x = this->xt;
            pt.y = this->yt;
            pt.z = this->zt;
        }

        void XDir(Point3<T>&pt) const
        {
            pt.x = this->xx;
            pt.y = this->xy;
            pt.z = this->xz;
        }

        void YDir(Point3<T>&pt) const
        {
            pt.x = this->yx;
            pt.y = this->yy;
            pt.z = this->yz;
        }

        void ZDir(Point3<T>&pt) const
        {
            pt.x = this->zx;
            pt.y = this->zy;
            pt.z = this->zz;
        }

        void SetXdir(Point3<T>&pt)
        {
            this->xx = pt.x;
            this->xy = pt.y;
            this->xz = pt.z;
        }

        void SetYdir(Point3<T>&pt)
        {
            this->yx = pt.x;
            this->yy = pt.y;
            this->yz = pt.z;
        }

        void SetZdir(Point3<T>&pt)
        {
            this->zx = pt.x;
            this->zy = pt.y;
            this->zz = pt.z;
        }

        void SetOrigin(Point3<T>&pt)
        {
            this->xt = pt.x;
            this->yt = pt.y;
            this->zt = pt.z;
        }
    };

    template<typename T>
    struct Line2
    {
        T a;
        T b;
        T c;

        constexpr Line2(T _a = 0, T _b = 0, T _c = 0) : a(_a), b(_b), c(_c) {}

        Line2(T params[3])
        {
            SetByEquationParams(params);
        }

        Line2(const Point2<T> &pt1, const Point2<T> &pt2)
        {
            SetByTwoPoints(pt1, pt2);
        }

        void SetByEquationParams(T params[3])
        {
            a = params[0];
            b = params[1];
            c = params[2];
        }

        void SetByTwoPoints(const Point2<T> &pt1, const Point2<T> &pt2)
        {
            a = pt1.y - pt2.y;
            b = pt2.x - pt1.x;
            c = -a * pt1.x - b * pt1.y;
        }

        void SetByPointAndDirection(const Point2<T> &pt, const Point2<T> &dir)
        {
            Point2<T> pt2 = pt + dir;
            SetByTwoPoints(pt, pt2);
        }

        bool Normalize()
        {
            T dirLen = sqrt(a * a + b * b);
            if (dirLen < (T)RG_ZERO)
            {
                return false;
            }
            a /= dirLen;
            b /= dirLen;
            c /= dirLen;
            return true;
        }

        Point2<T> Direction() const
        {
            return Point2<T>(b, -a);
        }

        T DistanceToPoint(const Point2<T>& pt, bool isNormalizedLine = false) const
        {
            if (isNormalizedLine)
                return fabs(a * pt.x + b * pt.y + c);
            else
                return fabs(a * pt.x + b * pt.y + c) / sqrt(a * a + b * b);
        }

        bool IsPointOnLine(const Point2<T>& pt, bool isNormalizedLine = false, T tol = (T)RG_ZERO) const
        {
            return DistanceToPoint(pt, isNormalizedLine) <= tol;
        }

        bool AreTwoPntsOnSameSide(const Point2<T>& pt1, const Point2<T>& pt2) const
        {
            T dis_1 = a * pt1.x + b * pt1.y + c;
            T dis_2 = a * pt2.x + b * pt2.y + c;
            return std::signbit(dis_1) == std::signbit(dis_2);
        }

        Point2<T> ProjectionOfPnt(const Point2<T>& pt, bool isNormalizedLine = false, T tol = (T)RG_ZERO) const
        {
            if (IsPointOnLine(pt, isNormalizedLine, tol))
                return pt;

            Point2<T> res;
            res.x = (b*b*pt.x - a * b*pt.y - a * c);
            res.y = (a*a*pt.y - a * b*pt.x - b * c);
            if (!isNormalizedLine)
            {
                T len = sqrt(a * a + b * b);
                res.x /= len;
                res.y /= len;
            }
            return res;
        }

        Line2<T> GetVerticalLine(const Point2<T>& cross_pnt) const
        {
            Line2<T> res;
            res.SetByPointAndDirection(cross_pnt, Point2<T>(a, b));
            return std::move(res);
        }
    };

    template<typename T>
    struct LineSegment2
    {
        Point2<T> start_pnt;
        Point2<T> end_pnt;
        Line2<T> equ;

        LineSegment2(const Point2<T> &startPnt, const Point2<T> &endPnt)
        {
            SetByTwoPoints(startPnt, endPnt);
        }

        void SetByTwoPoints(const Point2<T> &startPnt, const Point2<T> &endPnt)
        {
            start_pnt = startPnt;
            end_pnt = endPnt;
            equ.SetByTwoPoints(start_pnt, end_pnt);
            equ.Normalize();
        }

        T MinDistanceToPoint(const Point2<T> &pt)
        {
            Line2<T> verLine = equ.GetVerticalLine(pt);
            if (verLine.AreTwoPntsOnSameSide(start_pnt, end_pnt))
            {
                T res_1 = Point2<T>(pt - start_pnt).Length();
                T res_2 = Point2<T>(pt - end_pnt).Length();
                return res_1 < res_2 ? res_1 : res_2;
            }
            else
            {
                return equ.DistanceToPoint(pt, true);
            }
        }

        bool IsPointOnSegment(const Point2<T>& pt, T tol = (T)RG_ZERO) const
        {
            if (equ.DistanceToPoint(pt, true) > tol)
                return false;

            Point2<T> dir1 = pt - start_pnt;
            Point2<T> dir2 = pt - end_pnt;
            if (dir1.IsSameTol({ 0,0 }, tol)
                || dir2.IsSameTol({ 0,0 }, tol))
                return true;

            if (fabs(start_pnt.x - end_pnt.x) > tol)
            {
                return (std::signbit(dir1.x) != std::signbit(dir2.x));
            }
            else
            {
                return (std::signbit(dir1.y) != std::signbit(dir2.y));
            }
        }
    };

    template<typename T>
    struct Plane
    {
        T a;
        T b;
        T c;
        T d;

        constexpr Plane(T _a = 0, T _b = 0, T _c = 0, T _d = 0) : a(_a), b(_b), c(_c), d(_d) {}

        Plane(T params[4])
        {
            SetByEquationParams(params);
        }

        Plane(const Point3<T> &pt1, const Point3<T> &pt2, const Point3<T> &pt3)
        {
            SetByThreePoints(pt1, pt2, pt3);
        }

        Plane(const Point3<T> &pt, const Point3<T> &dir)
        {
            SetByPointAndDirection(pt, dir);
        }

        void SetByEquationParams(T params[4])
        {
            a = params[0];
            b = params[1];
            c = params[2];
            d = params[3];
        }

        void SetByThreePoints(const Point3<T> &pt1, const Point3<T> &pt2, const Point3<T> &pt3)
        {
            Point3<T> dir = Cross(pt2 - pt1, pt3 - pt2);
            a = dir.x;
            b = dir.y;
            c = dir.z;
            d = -Dot(pt1, dir);
        }

        void SetByPointAndDirection(const Point3<T> &pt, const Point3<T> &dir)
        {
            a = dir.x;
            b = dir.y;
            c = dir.z;
            d = -Dot(pt, dir);
        }

        Point3<T> Direction() const
        {
            return Point3<T>(a, b, c);
        }

        Point3<T> PointOnPlane() const
        {
            Point3<T> pt(0, 0, 0);
            if (!RG_ISZERO(a))
            {
                pt.x = -d / a;
            }
            else if (!RG_ISZERO(b))
            {
                pt.y = -d / b;
            }
            else if (!RG_ISZERO(c))
            {
                pt.z = -d / c;
            }
            return pt;
        }

        bool Normalize()
        {
            T dirLen = Direction().Length();
            if (dirLen < (T)RG_ZERO)
            {
                return false;
            }
            a /= dirLen;
            b /= dirLen;
            c /= dirLen;
            d /= dirLen;
            return true;
        }

        void Revese()
        {
            a *= -1;
            b *= -1;
            c *= -1;
        }

        bool IsSameTo(Plane<T> p, T tol = (T)RG_ZERO) const
        {
            Plane<T> p0(a, b, c, d);
            p0.Normalize();
            p.Normalize();

            if (fabs(p0.a - p.a) <= tol
                && fabs(p0.b - p.b) <= tol
                && fabs(p0.c - p.c) <= tol
                && fabs(p0.d - p.d) <= tol)
                return true;

            return false;
        }

        bool IsIntersectWith(Plane<T> p, T tol = (T)RG_ZERO) const
        {
            return !IsParallelTo(p, tol);
        }

        bool IsParallelTo(Plane<T> p, T tol = (T)RG_ZERO) const
        {
            Plane<T> p0(a, b, c, d);
            p0.Normalize();
            p.Normalize();

            if (fabs(p0.a - p.a) <= tol
                && fabs(p0.b - p.b) <= tol
                && fabs(p0.c - p.c) <= tol
                && fabs(p0.d - p.d) > tol)
                return true;

            if (fabs(p0.a + p.a) <= tol
                && fabs(p0.b + p.b) <= tol
                && fabs(p0.c + p.c) <= tol
                && fabs(p0.d - p.d) > tol)
                return true;

            return false;
        }

        T DistanceToPoint(const Point3<T>& pt, bool isNormalizedPlane = false) const
        {
            if (isNormalizedPlane)
                return fabs(a * pt.x + b * pt.y + c * pt.z + d);
            else
                return fabs(a * pt.x + b * pt.y + c * pt.z + d) / sqrt(a * a + b * b + c * c);
        }

        T SignedDistanceToPoint(const Point3<T>& pt, bool isNormalizedPlane = false) const
        {
            if (isNormalizedPlane)
                return (a * pt.x + b * pt.y + c * pt.z + d);
            else
                return (a * pt.x + b * pt.y + c * pt.z + d) / sqrt(a * a + b * b + c * c);
        }

        bool IsPointOnPlane(const Point3<T>& pt, T tol = (T)RG_ZERO, bool isNormalizedPlane = false) const
        {
            return DistanceToPoint(pt, isNormalizedPlane) <= tol;
        }

        bool IsIntersectWith(const Lim3<T>& lim3, T tol = (T)RG_ZERO) const
        {
            Plane<T> p(a, b, c, d);
            p.Normalize();

            Point3<T> vtxs[8];
            vtxs[0] = { lim3.x.min, lim3.y.min, lim3.z.min };
            vtxs[1] = { lim3.x.max, lim3.y.min, lim3.z.min };
            vtxs[2] = { lim3.x.max, lim3.y.max, lim3.z.min };
            vtxs[3] = { lim3.x.min, lim3.y.max, lim3.z.min };
            vtxs[4] = { lim3.x.min, lim3.y.min, lim3.z.max };
            vtxs[5] = { lim3.x.max, lim3.y.min, lim3.z.max };
            vtxs[6] = { lim3.x.max, lim3.y.max, lim3.z.max };
            vtxs[7] = { lim3.x.min, lim3.y.max, lim3.z.max };

            T signedDistance = p.SignedDistanceToPoint(vtxs[0], true);
            if (fabs(signedDistance) <= tol)
                return true;

            bool inFrontOfPlane = signedDistance > 0;
            for (int i = 1; i <= 7; i++)
            {
                signedDistance = p.SignedDistanceToPoint(vtxs[i], true);
                if (fabs(signedDistance) <= tol)
                    return true;

                if ((signedDistance > 0) != inFrontOfPlane)
                    return true;
            }
            return false;
        }
    };

    template<typename T>
    struct Triangle
    {
        Point3<T> pts[3]; /* triangle                                               */
    };

    template<typename T>
    struct Segment
    {
        Point3<T> pts[2]; /* segment                                               */
    };

    template<typename T>
    void Lim2AddPoint(Lim2<T> &lim, const Point2<T> &point)
    {
        if (lim.x.min > point.x) lim.x.min = point.x;
        if (lim.x.max < point.x) lim.x.max = point.x;
        if (lim.y.min > point.y) lim.y.min = point.y;
        if (lim.y.max < point.y) lim.y.max = point.y;
    }

    template<typename T>
    void Lim3AddPoint(Lim3<T> &lim, const Point3<T> &point)
    {
        if (lim.x.min > point.x) lim.x.min = point.x;
        if (lim.x.max < point.x) lim.x.max = point.x;
        if (lim.y.min > point.y) lim.y.min = point.y;
        if (lim.y.max < point.y) lim.y.max = point.y;
        if (lim.z.min > point.z) lim.z.min = point.z;
        if (lim.z.max < point.z) lim.z.max = point.z;
    }

    template<typename T>
    T Lim3Length(const Lim3<T> &lim)
    {
        return ZwPoint3(lim.x.max - lim.x.min, lim.y.max - lim.y.min, lim.z.max - lim.z.min).Length();
    }
}

#endif

