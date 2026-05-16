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

static int i, j, k, l, Hits = 0;
static double Xmin = -5.0, Xmax = 15.0, Ymin = -5.0, Ymax = 5.0, Nx = 790, Ny = 1185, Nc, junk, Px, Py;
static double Kx[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], Ky[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y];
static double X[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], Y[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], Cx[MAX_ARRAY_SIZE], Cy[MAX_ARRAY_SIZE], Area[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], Fl[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], Bl[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y];
static double Cxmax = 0.0, Cymax = 0.0, Cxmin = Xmax, Cymin = Ymax, dist, Geomx[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], Geomy[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y];
static double dx, dy, OutPt, PtonBndry, RadiA, RadiB, Radi, minradius = 0, Rad[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], Gst[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], Mark[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y];
static double X1, X2, Y1, Y2, Qx, Qy, ImX[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], ImY[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], ImNx[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], ImNy[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], ImBx[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], ImBy[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y];
static double ImGx[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y], ImGy[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE_Y];
int main()
{
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
    /*** CELL CENTRE CALCULATION
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
    for (i = 1; i <= Nx; i++)
    {
        for (j = 1; j <= Ny; j++)
        {
            for (k = 1; k < Nc; k++)
            {
                PtonBndry = PointOnLine(Cx[k], Cy[k], Cx[k + 1], Cy[k + 1], Kx[i][j], Ky[i][j]);
                if (PtonBndry == 1)
                {
                    Fl[i][j] = 1;
                    Bl[i][j] = 1.0;
                    break;
                }
            }
        }
    }
    for (i = 1; i <= Nx; i++)
    {
        for (j = 1; j <= Ny; j++)
        {
            if (Bl[i][j] == 0)
            {
                Hits = 0;
                l = Nc - 1;
                for (k = 1; k < (Nc); k++)
                {
                    Hits = Hits + HitOrMiss(Cx[l], Cy[l], Cx[k], Cy[k], Kx[i][j], Ky[i][j]);
                    l = k;
                }
                if ((Hits == 0) || ((Hits != 0) && ((Hits % 2) == 0.0)))
                {
                    Fl[i][j] = 1.0;
                }
                if (Hits != 0 && (Hits % 2) != 0)
                {
                    Fl[i][j] = -1.0;
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
            for (k = 1; k < Nc; k++)
            {
                RadiA = Radius(Cx[k], Cy[k], Cx[k + 1], Cy[k + 1], Kx[i][j], Ky[i][j]);
                RadiB = Projector(Cx[k], Cy[k], Cx[k + 1], Cy[k + 1], Kx[i][j], Ky[i][j]);
                Radi = Minimum2(RadiA, RadiB);
                // cout<<" RadA : "<<RadiA<<" RadB : "<<RadiB<<"Radi : "<<Radi<<"\t"<<endl;
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
                dist = sqrt((Cx[k + 1] - Cx[k]) * (Cx[k + 1] - Cx[k]) + (Cy[k + 1] - Cy[k]) * (Cy[k + 1] - Cy[k]));
                ImNx[i][j] = -(Cy[k + 1] - Cy[k]) / dist;
                ImNy[i][j] = (Cx[k + 1] - Cx[k]) / dist;
                // cout<< Mark[i][j]<<" "<<ImNx[i][j]<<"  "<<ImNy[i][j]<<"  "<<endl;
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
                Qx = FootX(Cx[k], Cy[k], Cx[k + 1], Cy[k + 1], Kx[i][j], Ky[i][j]);
                Geomx[i][j] = Qx;
                Qy = FootY(Cx[k], Cy[k], Cx[k + 1], Cy[k + 1], Kx[i][j], Ky[i][j]);
                Geomy[i][j] = Qy;
                ImX[i][j] = 2 * Qx - Kx[i][j];
                ImY[i][j] = 2 * Qy - Ky[i][j];
                counter = counter + 1;
                // cout<< Rad[i][j]<<"  "<<Mark[i][j]<<" "<<ImX[i][j]<<"  "<<ImY[i][j]<<"  "<<endl;
            }
        }
    }
    /*****************************************************************************************/
    /** IMAGE NEIGHBOUR POINTS  */
    /*****************************************************************************************/
    for (i = 1; i <= Nx; i++)
    {
        for (j = 1; j <= Ny; j++)
        {
            // i = 3; j =2;
            if (Gst[i][j] == 1)
            {
                for (int im = 1; im <= Nx; im++)
                {
                    for (int jm = 1; jm <= Ny; jm++)
                    {
                        if ((Kx[im][jm] <= ImX[i][j]) && (Kx[im + 1][jm] > ImX[i][j]) && ((Ky[im][jm] <= ImY[i][j]) && (Ky[im][jm + 1] > ImY[i][j])))
                        {
                            ImBx[i][j] = im;
                            ImBy[i][j] = jm;
                        }
                        if (((Kx[im][jm] <= Geomx[i][j]) && (Kx[im + 1][jm] > Geomx[i][j])) && ((Ky[im][jm] <= Geomy[i][j]) && (Ky[im][jm + 1] > Geomy[i][j])))
                        {
                            ImGx[i][j] = im;
                            ImGy[i][j] = jm;
                        }
                    }
                }
            } // cout<<"Wow "<<im<<" "<<jm<<" "<<i<<" "<<j<<endl;
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
