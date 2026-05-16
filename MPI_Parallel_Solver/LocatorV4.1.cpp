/*****************************************************************************************/
/** RAY TRACING CODE - Code outputs original Geom files for final interpolation in solver*/
/** Shrinath K S */
/** Aero- IISC */
/*****************************************************************************************/
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cstdlib>

// Constants - Better naming to avoid conflicts
const int MAX_ARRAY_SIZE = 3000;
const int MAX_ARRAY_SIZE_Y = 3000;
#include "MathFunctions.h"
#include "LocatorFunctions.h"
using namespace std;

static int i, j, k, l __attribute__((unused)), Hits = 0;
static double Xmin = -5.0, Xmax = 15.0, Ymin = -5.0, Ymax = 5.0, Nx = 790, Ny = 1185, Nc;
static double Kx[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], Ky[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y];
static double X[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], Y[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], Cx[MAX_ARRAY_SIZE], Cy[MAX_ARRAY_SIZE], Fl[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], Bl[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y];
static double Cxmax = 0.0, Cymax = 0.0, Cxmin = Xmax, Cymin = Ymax, dist, Geomx[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], Geomy[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y];
static double dx, dy, PtonBndry, RadiA, RadiB, Radi, minradius = 0, Rad[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], Gst[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], Mark[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y];
static double Qx, Qy, ImX[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], ImY[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], ImNx[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], ImNy[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], ImBx[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], ImBy[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y];
static double ImGx[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], ImGy[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y];
int main(int argc, char *argv[])
{
    if (argc >= 3)
    {
        Nx = std::atof(argv[1]);
        Ny = std::atof(argv[2]);
    }
    if (argc >= 7)
    {
        Xmin = std::atof(argv[3]);
        Xmax = std::atof(argv[4]);
        Ymin = std::atof(argv[5]);
        Ymax = std::atof(argv[6]);
    }
    if (argc >= 3)
    {
        std::cout << "Grid: Nx=" << Nx << " Ny=" << Ny
                  << " domain [" << Xmin << "," << Xmax << "] x ["
                  << Ymin << "," << Ymax << "]\n";
    }
    /*****************************************************************************************/
    /** GEOMETRY */
    /*****************************************************************************************/
    ifstream Geometry;
    Geometry.open("Geometry.dat");
    cout << "Reading Points" << endl;
    while (!Geometry.eof())
    {
        Geometry >> Nc;
        for (int i = 1; i <= Nc; i++)
        {
            Geometry >> Cx[i] >> Cy[i];
            cout << Cx[i] << "\t" << Cy[i] << "\t" << i << "/" << Nc << endl;
        }
    }
    ofstream Geom("Geom.dat");
    Geom << "TITLE = Boundary" << endl
         << "VARIABLES = X, Y,Flag" << endl;
    Geom << "Zone T = Omega I = " << Nc << endl;
    for (k = 1; k <= (Nc); k++)
    {
        Geom << Cx[k] << "\t" << Cy[k] << "\t" << "0.0" << endl;
    }
    Geom.close();

    for (k = 1; k <= Nc; k++)
    {
        if (Cxmax < Cx[k])
        {
            Cxmax = Cx[k];
        }
        if (Cxmin > Cx[k])
        {
            Cxmin = Cx[k];
        }
        if (Cymax < Cy[k])
        {
            Cymax = Cy[k];
        }
        if (Cymin > Cy[k])
        {
            Cymin = Cy[k];
        }
    }
    cout << "CxMAX " << Cxmax << "  CxMin " << Cxmin << "\t" << endl;
    /*****************************************************************************************/
    /** CARTESIAN GRID */
    /*****************************************************************************************/
    ofstream Grido("Grid.dat");
    Grido << "TITLE = Flow" << endl
          << "VARIABLES = X, Y,Boundary, Flag" << endl;
    Grido << "Zone T = Omega I = " << Nx + 1 << " J = " << Ny + 1 << endl;
    dx = (Xmax - Xmin) / Nx;
    dy = (Ymax - Ymin) / Ny;
    for (j = 1; j <= Ny + 1; j++)
    {
        for (i = 1; i <= Nx + 1; i++)
        {
            X[i][j] = Xmin + (i - 1) * dx;
            Y[i][j] = Ymin + (j - 1) * dy;
            Grido << X[i][j] << "\t" << Y[i][j] << "\t" << "0.0" << "\t" << "0.0" << endl;
            // cout<<X[i][j]<<"\t"<<Y[i][j]<<"\t"<<dx<<"\t"<<"\n";
        }
    }
    Grido.close();
    /***********************************************************************************************/
    /** CELL CENTRE CALCULATION */
    /***********************************************************************************************/
    for (i = 1; i <= Nx; i++)
    {
        for (j = 1; j <= Ny; j++)
        {
            Kx[i][j] = (X[i][j] + X[i + 1][j] + X[i + 1][j + 1] + X[i][j + 1]) / 4.0;
            Ky[i][j] = (Y[i][j] + Y[i + 1][j] + Y[i + 1][j + 1] + Y[i][j + 1]) / 4.0;
        }
    }
    /*****************************************************************************************/
    /** CELL CENTRE LOCATOR */
    /*****************************************************************************************/
    /* Polygon segments are walked with closure: edge k connects Cx[k] to Cx[k % Nc + 1]. */
    for (i = 1; i <= Nx; i++)
    {
        for (j = 1; j <= Ny; j++)
        {
            for (k = 1; k <= Nc; k++)
            {
                int kp1 = (k % (int)Nc) + 1;
                PtonBndry = PointOnLine(Cx[k], Cy[k], Cx[kp1], Cy[kp1], Kx[i][j], Ky[i][j]);
                if (PtonBndry == 1)
                {
                    Fl[i][j] = 1;
                    Bl[i][j] = 1.0;
                    break;
                }
            }
        }
    }
    /* Inside/outside via PNPOLY ray casting; loop closes the polygon. */
    for (i = 1; i <= Nx; i++)
    {
        for (j = 1; j <= Ny; j++)
        {
            if (Bl[i][j] == 0)
            {
                Hits = 0;
                for (k = 1; k <= Nc; k++)
                {
                    int kp1 = (k % (int)Nc) + 1;
                    Hits = Hits + HitOrMiss(Cx[k], Cy[k], Cx[kp1], Cy[kp1], Kx[i][j], Ky[i][j]);
                }
                if ((Hits % 2) == 0)
                {
                    Fl[i][j] = 1.0;   /* outside body -> fluid */
                }
                else
                {
                    Fl[i][j] = -1.0;  /* inside body  -> solid */
                }
            }
        }
    }
    /*****************************************************************************************/
    /** Radius  */
    /*****************************************************************************************/
    for (i = 1; i <= Nx; i++)
    {
        for (j = 1; j <= Ny; j++)
        {
            minradius = 10;
            for (k = 1; k <= Nc; k++)
            {
                int kp1 = (k % (int)Nc) + 1;
                RadiA = Radius(Cx[k], Cy[k], Cx[kp1], Cy[kp1], Kx[i][j], Ky[i][j]);
                RadiB = Projector(Cx[k], Cy[k], Cx[kp1], Cy[kp1], Kx[i][j], Ky[i][j]);
                Radi = Minimum2(RadiA, RadiB);
                if (minradius > Radi)
                {
                    minradius = Radi;
                    Mark[i][j] = k;
                }
            }
            Rad[i][j] = minradius;
        }
    }
    /*****************************************************************************************/
    /** Ghost Cell marker  */
    /*****************************************************************************************/
    for (i = 1; i <= Nx; i++)
    {
        for (j = 1; j <= Ny; j++)
        {
            if ((Rad[i][j] >= (1.0 * dx)) && (Fl[i][j] == 1))
            {
                Gst[i][j] = 0;
            }
            if ((Rad[i][j] < (1.0 * dx)) && (Fl[i][j] == 1))
            {
                Gst[i][j] = 0;
            }
            if ((Rad[i][j] < (1.0 * dx)) && (Fl[i][j] == -1))
            {
                Gst[i][j] = 1;
                Rad[i][j] = -Rad[i][j];
            }
            if ((Rad[i][j] >= (1.0 * dx)) && (Fl[i][j] == -1))
            {
                Gst[i][j] = 0;
                Rad[i][j] = -Rad[i][j];
            }
        }
    }
    /*****************************************************************************************/
    /** Ghost Cell Normals  */
    /*****************************************************************************************/
    for (i = 1; i <= Nx; i++)
    {
        for (j = 1; j <= Ny; j++)
        {
            if (Gst[i][j] == 1)
            {
                k = Mark[i][j];
                int kp1 = (k % (int)Nc) + 1;
                dist = sqrt((Cx[kp1] - Cx[k]) * (Cx[kp1] - Cx[k]) + (Cy[kp1] - Cy[k]) * (Cy[kp1] - Cy[k]));
                ImNx[i][j] = -(Cy[kp1] - Cy[k]) / dist;
                ImNy[i][j] = (Cx[kp1] - Cx[k]) / dist;
            }
        }
    }
    /*****************************************************************************************/
    /** IMAGE POINTS  */
    /*****************************************************************************************/
    int counter = 0;
    for (i = 1; i <= Nx; i++)
    {
        for (j = 1; j <= Ny; j++)
        {
            if (Gst[i][j] == 1)
            {
                k = Mark[i][j];
                int kp1 = (k % (int)Nc) + 1;
                Qx = FootX(Cx[k], Cy[k], Cx[kp1], Cy[kp1], Kx[i][j], Ky[i][j]);
                Geomx[i][j] = Qx;
                Qy = FootY(Cx[k], Cy[k], Cx[kp1], Cy[kp1], Kx[i][j], Ky[i][j]);
                Geomy[i][j] = Qy;
                ImX[i][j] = 2 * Qx - Kx[i][j];
                ImY[i][j] = 2 * Qy - Ky[i][j];
                counter = counter + 1;
            }
        }
    }
    /*****************************************************************************************/
    /** IMAGE NEIGHBOUR POINTS - O(1) per ghost cell                                          */
    /*                                                                                       */
    /* Kx, Ky are cell centres of a uniform Cartesian grid:                                  */
    /*   Kx[i][j] = Xmin + (i - 0.5) * dx,  Ky[i][j] = Ymin + (j - 0.5) * dy                 */
    /* The cell whose lower-left corner brackets a query point (px, py) is therefore         */
    /*   im = floor((px - (Xmin + 0.5*dx)) / dx) + 1                                         */
    /*   jm = floor((py - (Ymin + 0.5*dy)) / dy) + 1                                         */
    /* Result clamped to [1, Nx-1] x [1, Ny-1] so the bilinear stencil (im..im+1, jm..jm+1)  */
    /* always stays in range. This replaces an O(N^4) scan with O(1) per ghost.              */
    /*****************************************************************************************/
    {
        const double inv_dx = 1.0 / dx;
        const double inv_dy = 1.0 / dy;
        const double x0 = Xmin + 0.5 * dx;
        const double y0 = Ymin + 0.5 * dy;
        const int max_i = (int)Nx - 1;
        const int max_j = (int)Ny - 1;
        for (i = 1; i <= Nx; i++)
        {
            for (j = 1; j <= Ny; j++)
            {
                if (Gst[i][j] == 1)
                {
                    int im = (int)std::floor((ImX[i][j] - x0) * inv_dx) + 1;
                    int jm = (int)std::floor((ImY[i][j] - y0) * inv_dy) + 1;
                    if (im < 1)      im = 1;
                    if (jm < 1)      jm = 1;
                    if (im > max_i)  im = max_i;
                    if (jm > max_j)  jm = max_j;
                    ImBx[i][j] = im;
                    ImBy[i][j] = jm;

                    int gim = (int)std::floor((Geomx[i][j] - x0) * inv_dx) + 1;
                    int gjm = (int)std::floor((Geomy[i][j] - y0) * inv_dy) + 1;
                    if (gim < 1)     gim = 1;
                    if (gjm < 1)     gjm = 1;
                    if (gim > max_i) gim = max_i;
                    if (gjm > max_j) gjm = max_j;
                    ImGx[i][j] = gim;
                    ImGy[i][j] = gjm;
                }
            }
        }
    }
    cout << "Ghost cells have been Located" << endl;
    /*for( k=1;k<Nc;k++){
        for( i=1;i<=Nx;i++){
            for( j=1;j<=Ny;j++){
             if ((Kx[i][j]<=Cx[k])&&(Kx[i+1][j]>Cx[k])&&((Ky[i][j]<=Cy[k])&&(Ky[i][j+1]>Cy[k])))
                {
                 ImGx[k]=i;ImGy[k]=j; cout<< Kx[i][j]<<" "<<Cx[k]<<" "<<i<<" "<<j<<" "<<k<<endl;
                }
          }
       }
    }*/
    // cout<<"Geom Points Sorted with Ghost cell"<<endl;
    /*****************************************************************************************/
    /** Output  */
    /*****************************************************************************************/
    ofstream Grids("GridPoints.dat");
    Grids << Nx + 1 << "\t" << Ny + 1 << "\t" << endl;
    for (j = 1; j <= Ny + 1; j++)
    {
        for (i = 1; i <= Nx + 1; i++)
        {
            Grids << X[i][j] << "\t" << Y[i][j] << endl;
        }
    }
    Grids.close();
    cout << "Written Grid points.dat" << endl;

    ofstream Boundary("Boundary.dat");
    Boundary << "TITLE = Boundary" << endl
             << "VARIABLES = X, Y, Flag, Radii,Gst" << endl;
    Boundary << "Zone T = Omega I = " << Nx << " J = " << Ny << endl;
    for (j = 1; j <= Ny; j++)
    {
        for (i = 1; i <= Nx; i++)
        {
            Boundary << Kx[i][j] << "\t" << Ky[i][j] << "\t" << Fl[i][j] << "\t" << Rad[i][j] << "\t" << Gst[i][j] << endl;
        }
    }
    Boundary.close();
    cout << "Written Boundary points.dat" << endl;
    /*****************************************************************************************/
    ofstream Ibound("ImmersedBoundary.dat");
    Ibound << Nx << "\t" << Ny << "\t" << endl;
    for (j = 1; j <= Ny; j++)
    {
        for (i = 1; i <= Nx; i++)
        {
            Ibound << Kx[i][j] << "\t" << Ky[i][j] << "\t" << Fl[i][j] << "\t" << Gst[i][j] << "\t" << Rad[i][j] << "\t" << ImNx[i][j] << "\t" << ImNy[i][j] << "\t" << ImX[i][j] << "\t" << ImY[i][j] << "\t" << ImBx[i][j] << "\t" << ImBy[i][j] << "\t" << Geomx[i][j] << "\t" << Geomy[i][j] << "\t" << ImGx[i][j] << "\t" << ImGy[i][j] << endl;
        }
    }
    cout << "Written Immersed Boundary points.dat" << endl;
    Ibound.close();
    ofstream IGeom("ImmersedGeometry.dat");
    IGeom << counter << "\t" << endl;
    for (j = 1; j <= Ny; j++)
    {
        for (i = 1; i <= Nx; i++)
        {
            if (Gst[i][j] == 1)
            {
                IGeom << Geomx[i][j] << "\t" << Geomy[i][j] << "\t" << ImGx[i][j] << "\t" << ImGy[i][j] << endl;
            }
        }
    }
    IGeom.close();
    cout << "Writing IGeom.dat" << endl;
    /*****************************************************************************************/
    cout << "***************END OF RUN*****************";
    return 0;
}
