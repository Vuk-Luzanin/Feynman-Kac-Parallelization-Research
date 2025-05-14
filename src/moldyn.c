#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "util.h"

#define MM 15
#define NPART 4 * MM *MM *MM

double epot;
double vir;
double count;

void dfill(int n, double val, double a[], int ia)
{
    for (int i = 0; i < (n - 1) * ia + 1; i += ia)
    {
        a[i] = val;
    }
}

void domove(int n3, double x[], double vh[], double f[], double side)
{
    for (int i = 0; i < n3; i++)
    {
        x[i] += vh[i] + f[i];
        // Periodic boundary conditions
        if (x[i] < 0.0)
            x[i] += side;
        if (x[i] > side)
            x[i] -= side;
        // Partial velocity updates
        vh[i] += f[i];
        // Initialise forces for the next iteration
        f[i] = 0.0;
    }
}

void dscal(int n, double sa, double sx[], int incx)
{
    if (incx == 1)
    {
        for (int i = 0; i < n; i++)
            sx[i] *= sa;
    }
    else
    {
        int j = 0;
        for (int i = 0; i < n; i++)
        {
            sx[j] *= sa;
            j += incx;
        }
    }
}

void fcc(double x[], int mm, double a)
{
    int ijk = 0;

    for (int lg = 0; lg < 2; lg++)
        for (int i = 0; i < mm; i++)
            for (int j = 0; j < mm; j++)
                for (int k = 0; k < mm; k++)
                {
                    x[ijk] = i * a + lg * a * 0.5;
                    x[ijk + 1] = j * a + lg * a * 0.5;
                    x[ijk + 2] = k * a;
                    ijk += 3;
                }

    for (int lg = 1; lg < 3; lg++)
        for (int i = 0; i < mm; i++)
            for (int j = 0; j < mm; j++)
                for (int k = 0; k < mm; k++)
                {
                    x[ijk] = i * a + (2 - lg) * a * 0.5;
                    x[ijk + 1] = j * a + (lg - 1) * a * 0.5;
                    x[ijk + 2] = k * a + a * 0.5;
                    ijk += 3;
                }
}

void forces(int npart, double x[], double f[], double side, double rcoff)
{
    vir = 0.0;
    epot = 0.0;
    double sideh = 0.5 * side;
    double rcoffs = rcoff * rcoff;

#pragma omp parallel for collapse(2) default(none) shared(npart, x, f, side, sideh, rcoffs) reduction(+                  \
                                                                                                      : vir) reduction(+ \
                                                                                                                       : epot)
    for (int i = 0; i < npart * 3; i += 3)
    {
        for (int j = i + 3; j < npart * 3; j += 3)
        {
            double xx = x[i] - x[j];
            double yy = x[i + 1] - x[j + 1];
            double zz = x[i + 2] - x[j + 2];
            if (xx < -sideh)
                xx += side;
            if (xx > sideh)
                xx -= side;
            if (yy < -sideh)
                yy += side;
            if (yy > sideh)
                yy -= side;
            if (zz < -sideh)
                zz += side;
            if (zz > sideh)
                zz -= side;
            double rd = xx * xx + yy * yy + zz * zz;

            if (rd <= rcoffs)
            {
                double rrd = 1.0 / rd;
                double rrd2 = rrd * rrd;
                double rrd3 = rrd2 * rrd;
                double rrd4 = rrd2 * rrd2;
                double rrd6 = rrd2 * rrd4;
                double rrd7 = rrd6 * rrd;
                epot += rrd6 - rrd3;
                double r148 = rrd7 - 0.5 * rrd4;
                vir -= rd * r148;
                double forcex = xx * r148;
                double forcey = yy * r148;
                double forcez = zz * r148;
#pragma omp atomic
                f[i] += forcex;
#pragma omp atomic
                f[j] -= forcex;
#pragma omp atomic
                f[i + 1] += forcey;
#pragma omp atomic
                f[j + 1] -= forcey;
#pragma omp atomic
                f[i + 2] += forcez;
#pragma omp atomic
                f[j + 2] -= forcez;
            }
        }
    }
}

double mkekin(int npart, double f[], double vh[], double hsq2, double hsq)
{
    double sum = 0.0;

    for (int i = 0; i < 3 * npart; i++)
    {
        f[i] *= hsq2;
        vh[i] += f[i];
        sum += vh[i] * vh[i];
    }

    return sum / hsq;
}

// Sample Maxwell distribution at temperature tref
void mxwell(double vh[], int n3, double h, double tref)
{
    int npart = n3 / 3;
    double ekin = 0.0;
    double sp = 0.0;

    srand48(4711);
    double tscale = 16.0 / ((double)npart - 1.0);

    for (int i = 0; i < n3; i += 2)
    {
        double v1;
        double v2;
        double s = 2.0;
        while (s >= 1.0)
        {
            v1 = 2.0 * drand48() - 1.0;
            v2 = 2.0 * drand48() - 1.0;
            s = v1 * v1 + v2 * v2;
        }
        double r = sqrt(-2.0 * log(s) / s);
        vh[i] = v1 * r;
        vh[i + 1] = v2 * r;
    }

    for (int i = 0; i < n3; i += 3)
        sp += vh[i];
    sp /= (double)npart;
    for (int i = 0; i < n3; i += 3)
    {
        vh[i] -= sp;
        ekin += vh[i] * vh[i];
    }

    sp = 0.0;
    for (int i = 1; i < n3; i += 3)
        sp += vh[i];
    sp /= (double)npart;
    for (int i = 1; i < n3; i += 3)
    {
        vh[i] -= sp;
        ekin += vh[i] * vh[i];
    }

    sp = 0.0;
    for (int i = 2; i < n3; i += 3)
        sp += vh[i];
    sp /= (double)npart;
    for (int i = 2; i < n3; i += 3)
    {
        vh[i] -= sp;
        ekin += vh[i] * vh[i];
    }

    double sc = h * sqrt(tref / (tscale * ekin));
    for (int i = 0; i < n3; i++)
        vh[i] *= sc;
}

void prnout(int move, double ekin, double epot, double tscale, double vir, double vel, double count, int npart, double den)
{
    double ek = 24.0 * ekin;
    epot *= 4.0;
    double etot = ek + epot;
    double temp = tscale * ekin;
    double pres = den * 16.0 * (ekin - vir) / (double)npart;
    vel /= (double)npart;
    double rp = (count / (double)npart) * 100.0;
    printf(" %6d%12.4f%12.4f%12.4f%10.4f%10.4f%10.4f%6.1f\n", move, ek, epot, etot, temp, pres, vel, rp);
}

double velavg(int npart, double vh[], double vaver, double h)
{
    double vaverh = vaver * h;
    double vel = 0.0;

    count = 0.0;
    for (int i = 0; i < npart * 3; i += 3)
    {
        double sq = sqrt(vh[i] * vh[i] + vh[i + 1] * vh[i + 1] + vh[i + 2] * vh[i + 2]);
        if (sq > vaverh)
            count++;
        vel += sq;
    }
    vel /= h;

    return vel;
}

int main(void)
{
    double x[NPART * 3];
    double vh[NPART * 3];
    double f[NPART * 3];

    double den = 0.83134;
    double side = pow((double)NPART / den, 0.3333333);
    double tref = 0.722;
    double rcoff = (double)MM / 4.0;
    double h = 0.064;
    int irep = 10;
    int istop = 20;
    int movemx = 20;

    double a = side / (double)MM;
    double hsq = h * h;
    double hsq2 = hsq * 0.5;
    double tscale = 16.0 / ((double)NPART - 1.0);
    double vaver = 1.13 * sqrt(tref / 24.0);

    printf("TEST: num_threads=%ld\n", get_num_threads());

    // Generate fcc lattice for atoms inside box
    fcc(x, MM, a);
    // Initialize velocities and forces (which are zero in fcc positions)
    mxwell(vh, 3 * NPART, h, tref);
    dfill(3 * NPART, 0.0, f, 1);

    double time = 0.0;

    for (int move = 1; move <= movemx; ++move)
    {
        double start = omp_get_wtime();
        // Move the particles and partially update velocities
        domove(3 * NPART, x, vh, f, side);
        // Compute forces in the new positions and accumulate the virial
        // and potential energy.
        forces(NPART, x, f, side, rcoff);
        // Scale forces, complete update of velocities and compute k.e.
        double ekin = mkekin(NPART, f, vh, hsq2, hsq);
        // Average the velocity and temperature scale if desired
        double vel = velavg(NPART, vh, vaver, h);
        if (move < istop && fmod(move, irep) == 0)
        {
            double sc = sqrt(tref / (tscale * ekin));
            dscal(3 * NPART, sc, vh, 1);
            ekin = tref / tscale;
        }
        time += omp_get_wtime() - start;
        prnout(move, ekin, epot, tscale, vir, vel, count, NPART, den);
    }

    printf("%d    %lf\n", movemx, time);

    return 0;
}
