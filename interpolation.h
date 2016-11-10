#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <QString>
#include <QVector>
#include <QPointF>
#include <QDebug>

typedef QVector<QPointF> V_Nodes;
//-------------------------------------------------------------------------------------------------
class Interpolation
{
public:
    struct Spline_Factors {
        Spline_Factors() {}
        Spline_Factors(double _a, double _b, double _c, double _d)
            :   a(_a), b(_b), c(_c), d(_d) {}
        double a, b, c, d;
    };
    struct KSI_ETA {
        KSI_ETA() {}
        KSI_ETA(double _ksi, double _eta) : ksi(_ksi), eta(_eta) {}
        double ksi, eta;
    };

public:
    Interpolation();
    double forward(const double value, const int exp, const V_Nodes &) throw (QString);
    double inverse(const double value, const int exp, const V_Nodes &) throw (QString);

    double spline(const V_Nodes &, const double value) const throw (QString);

    void build_spline_tbl(const V_Nodes &);

private:
    double get_Mult(const double value, const int nummer, const V_Nodes & environ) const;
    double get_Diff(const int from, const int to, const V_Nodes & environ);

    QVector<Spline_Factors> all_factors;
};

//-------------------------------------------------------------------------------------------------
#endif // INTERPOLATION_H
