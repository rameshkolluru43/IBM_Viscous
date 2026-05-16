#ifndef LOCATORFUNCTIONS_H
#define LOCATORFUNCTIONS_H

/*****************************************************************************************/
/** Locator Functions */
/*****************************************************************************************/
/** If point between the region */
template <typename T>
T Between(T x1, T y1, T x2, T y2, T px, T py)
{
    if ((((x1 <= px) && (px <= x2)) || ((x2 <= px) && (px <= x1))) && (((y1 <= py) && (py <= y2)) || ((y2 <= py) && (py <= y1))))
        return 1;
    else
        return 0;
}
/** If point on the equation of the line */
template <typename T>
T PointOnLineEq(T x1, T y1, T x2, T y2, T px, T py)
{
    if ((y1 * (x2 - x1) + (y2 - y1) * (px - x1) - py * (x2 - x1)) == 0)
        return 1;
    else
        return 0;
}
/** If point on the line segment*/
template <typename T>
T PointOnLine(T x1, T y1, T x2, T y2, T px, T py)
{
    double A = PointOnLineEq(x1, y1, x2, y2, px, py);
    double B = Between(x1, y1, x2, y2, px, py);
    if ((A == 1) && (B == 1))
        return 1;
    else
        return 0;
}
/** Projection distance between the point and the Line segment(internal projection) */
template <typename T>
T Projector(T x1, T y1, T x2, T y2, T px, T py)
{
    double disty = 100, distx = 100, qx, qy, Marker = 0;
    if (((y2 <= py) && (py <= y1)) || ((y1 <= py) && (py <= y2)))
    {
        if (y1 != y2)
        {
            qx = (x2 + (py - y2) / (y1 - y2) * (x1 - x2));
            qy = py;
            distx = sqrt((qx - px) * (qx - px) + (qy - py) * (qy - py));
        }
        else
        {
            distx = Minimum2(fabs(x1 - px), fabs(x2 - px));
        }
        Marker = Marker + 1;
    }
    if (((x2 <= px) && (px <= x1)) || ((x1 <= px) && (px <= x2)))
    {
        if (x1 != x2)
        {
            qy = (y1 + (y2 - y1) * (px - x1) / (x2 - x1));
            qx = px;
            disty = sqrt((qx - px) * (qx - px) + (qy - py) * (qy - py));
        }
        else
        {
            disty = Minimum2(fabs(y1 - py), fabs(y2 - py));
        }
        Marker = Marker + 1;
    }
    double distc = sqrt((x1 - px) * (x1 - px) + (y1 - py) * (y1 - py));
    double distd = sqrt((x2 - px) * (x2 - px) + (y2 - py) * (y2 - py));
    if (Marker != 0)
    {
        return Minimum2(distx, disty);
    }
    else
    {
        return Minimum2(distc, distd);
    }
}

/** Perpendicular distance between the point and the Line segment */
template <typename T>
T Radius(T x1, T y1, T x2, T y2, T px, T py)
{
    double dist, A, B, C, qx, qy, RMS;
    A = -(y2 - y1);
    B = (x2 - x1);
    RMS = (A * A + B * B);
    C = x1 * (y2 - y1) - y1 * (x2 - x1);
    qx = (B * B * px - A * B * py - A * C) / RMS;
    if (B != 0)
    {
        qy = -(A * qx / B) - (C / B);
    }
    else
    {
        qy = py;
    }
    if ((((x1 <= qx) && (qx <= x2)) || ((x2 <= qx) && (qx <= x1))) && (((y1 <= qy) && (qy <= y2)) || ((y2 <= qy) && (qy <= y1))))
    {
        dist = sqrt((qx - px) * (qx - px) + (qy - py) * (qy - py));
        return dist;
    }
    else
    {
        return 1000;
    }
}
/** To get Image point */
template <typename T>
T FootX(T x1, T y1, T x2, T y2, T px, T py)
{
    double A, B, C, qx, RMS;
    A = -(y2 - y1);
    B = (x2 - x1);
    RMS = (A * A + B * B);
    C = x1 * (y2 - y1) - y1 * (x2 - x1);
    qx = (B * B * px - A * B * py - A * C) / RMS;
    return qx;
}
template <typename T>
T FootY(T x1, T y1, T x2, T y2, T px, T py)
{
    double A, B, C, qx, qy, RMS;
    A = -(y2 - y1);
    B = (x2 - x1);
    RMS = (A * A + B * B);
    C = x1 * (y2 - y1) - y1 * (x2 - x1);
    qx = (B * B * px - A * B * py - A * C) / RMS;
    if (B != 0)
    {
        qy = -(A * qx / B) - (C / B);
    }
    else
    {
        qy = py;
    }
    return qy;
}
/** Robust ray-casting test: does a horizontal ray to +x from (px,py) cross
 * the segment (x1,y1)-(x2,y2)?  Standard PNPOLY half-open interval check.
 * Returns 1 on a crossing, 0 otherwise.  Combined with a polygon-closing
 * loop, an even/odd count of crossings determines inside/outside.
 */
template <typename T>
T HitOrMiss(T x1, T y1, T x2, T y2, T px, T py)
{
    const bool straddles_y = ((y1 > py) != (y2 > py));
    if (!straddles_y)
        return 0;
    const T x_intersect = (x2 - x1) * (py - y1) / (y2 - y1) + x1;
    return (px < x_intersect) ? 1 : 0;
}
/**Bilinear Interpolation*/
template <typename T>
T Bilinear(T x11, T y11, T x21, T y21, T x12, T y12, T x22, T y22, T px, T py, T q11, T q21, T q12, T q22)
{
    double q, r1, r2;
    r1 = (x22 - px) * q11 / (x22 - x11) + (px - x11) * q21 / (x22 - x11);
    r2 = (x22 - px) * q12 / (x22 - x11) + (px - x11) * q22 / (x22 - x11);
    q = (y22 - py) * r1 / (y22 - y11) + (py - y11) * r2 / (y22 - y11);
    return q;
}

#endif // LOCATORFUNCTIONS_H
