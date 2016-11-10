#include "interpolation.h"

const double EPS = 0.0000000000000000001;
//-------------------------------------------------------------------------------------------------
QString toString(const double value) {
    return QString::number(value);
}

QString toString(const Interpolation::Spline_Factors & value) {
    return " {" + QString::number(value.a) + ", " + QString::number(value.b) + ", "
                + QString::number(value.c) + ", " + QString::number(value.d) + "}\n";
}

QString toString(const Interpolation::KSI_ETA & value) {
    return  " {" + QString::number(value.ksi) + ", " + QString::number(value.eta) + "}\n";
}

QString toString(const QPointF & value) {
    return  " {" + QString::number(value.x()) + ", " + QString::number(value.y()) + "}\n";
}

template <typename Value>
void debug_print(const QVector<Value> & out_vector) {
    QString str_vector = "";
    typename QVector<Value>::const_iterator It = out_vector.begin();
    for (; It != out_vector.end(); ++It) {
        str_vector += toString(*It);
        str_vector += ' ';
    }
    qDebug() << str_vector;
}

//-------------------------------------------------------------------------------------------------
Interpolation::Interpolation()
{}

//-------------------------------------------------------------------------------------------------
double Interpolation::forward(const double value, const int exp, const V_Nodes & environ)
    throw (QString)
{
    if (exp <= 1 || exp != environ.count())
        throw QString("Неверное кол-во элементов окружения. Этап calculate.");

    double summ = 0;
    double summand = 0;

    for (int cnt = 0; cnt < exp; ++cnt) {
        summand = get_Mult(value, cnt, environ) * get_Diff(0, cnt, environ);
        summ += summand;
    }

    return summ;
}

//-------------------------------------------------------------------------------------------------
double Interpolation::get_Mult(const double value, const int nummer, const V_Nodes & environ) const {
    double result = 1;
    double factor = 1;

    for (int cnt = 0; cnt < nummer; ++cnt) {
        factor = value - environ[cnt].x();
        result *= factor;
    }

    return result;
}

//-------------------------------------------------------------------------------------------------
double Interpolation::get_Diff(const int from, const int to, const V_Nodes & environ) {
    if (from == to) {
        return environ[from].y();
    }

    double result = 0;

    if (to - from > 1)
        result = get_Diff(from + 1, to, environ) - get_Diff(from, to - 1, environ);
    else if (to - from == 1) {
        result = environ[to].y() - environ[from].y();
    }

    result /= (environ[to].x() - environ[from].x());

    return result;
}

//-------------------------------------------------------------------------------------------------
double Interpolation::inverse(const double value, const int exp, const V_Nodes & environ) throw (QString) {
    V_Nodes inverse_list;

    V_Nodes::const_iterator It = environ.begin();
    for (; It < environ.end(); ++It)
        inverse_list.append(QPointF(It->y(), It->x()));

    try {
        return this->forward(value, exp, inverse_list);
    } catch (QString) {
        throw;
    }
}

double Interpolation::spline(const V_Nodes & all_nodes, const double value) const throw (QString) {
    if (all_factors.count() > 0) {
        for (int I = 0; I < all_nodes.size() - 1; ++I) {
            if (qAbs(all_nodes[I].x() - value) < EPS) {
                return all_nodes[I].y();
            }

            if (all_nodes[I].x() < value && all_nodes[I + 1].x() > value) {
                double interval_len = value - all_nodes[I].x();
                const Spline_Factors & cur_Factor = all_factors[I];

                double result = cur_Factor.a + cur_Factor.b * interval_len + cur_Factor.c * interval_len * interval_len
                              + cur_Factor.d * interval_len * interval_len * interval_len;

                //qDebug() << toString(cur_Factor) << "    " << result;
                return result;
            }
        }
        int last = all_nodes.count() - 1;
        if (qAbs(all_nodes[last].x() - value) < EPS) {
            return all_nodes[last].y();
        }
        throw QString("Значение находится вне таблицы");
    } else
        throw QString("Таблица коэффициентов сплайна не вычислена");
    return 0;
}

void Interpolation::build_spline_tbl(const V_Nodes & all_nodes) {
    all_factors.clear();
    int N = all_nodes.count() - 1;
    double H = all_nodes[1].x() - all_nodes[0].x();
    double KSI = 0, ETA = 0, prev_KSI = 0;
    double A = H, B = 0, D = 0, F = 0;

    QVector<double> Hs;
    Hs.append(H);
    QVector<KSI_ETA> interms;
    interms.append(KSI_ETA(KSI, ETA));

    for (int I = 2; I <= N; ++I) {
        A = H;
        H = all_nodes[I].x() - all_nodes[I - 1].x();
        B = -2. * (H + A);
        D = H;
        F = -3. * ((all_nodes[I].y() - all_nodes[I - 1].y()) / H - (all_nodes[I - 1].y() - all_nodes[I - 2].y()) / A);

        prev_KSI = KSI;
        KSI = D / (B - A * prev_KSI);
        ETA = (A * ETA + F) / (B - A * prev_KSI);

        interms.append(KSI_ETA(KSI, ETA));
        Hs.append(H);
    }

    A = 0; B = 0; D = 0;
    double C = 0, prev_C = C;
    for (int I = 1; I <= N; ++I)
        all_factors.append(Spline_Factors(0, 0, 0, 0));

    for (unsigned int I = N; I > 0; --I) {
        A = all_nodes[I - 1].y();
        prev_C = C;
        C = interms[I - 1].ksi * C + interms[I - 1].eta;
        B = (all_nodes[I].y() - A) / Hs[I - 1] - Hs[I - 1] / 3 * (2 * C + prev_C);
        D = (prev_C - C) / (3 * Hs[I - 1]);

        all_factors[I - 1] = Spline_Factors(A, B, C, D);
    }

    /*
    debug_print<QPointF>(all_nodes);
    debug_print<double>(Hs);
    debug_print<KSI_ETA>(interms);
    debug_print<Spline_Factors>(all_factors);
    */
}
//-------------------------------------------------------------------------------------------------
