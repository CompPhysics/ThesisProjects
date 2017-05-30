#ifndef PM_H
#define PM_H

class PM
{
public:
    PM();
    double  makeStateSpace  ();
    double  assym           (int p, int q, int r, int s);
    double  assym_single    (int p, int q);
    double  h0              (int p);
    double  f               (int p);
    int     kUnique1        (int p, int s1);
    int     kUnique2        (int p, int q, int s1, int s2);
    int     kUnique3        (int p, int q, int r, int s1, int s2, int s3);
    int     kUnique4        (int p, int q, int r, int s, int s1, int s2, int s3, int s4);
    int     kUnique5        (int p, int q, int r, int s, int t, int s1, int s2, int s3, int s4, int s5);
};

#endif // PM_H
