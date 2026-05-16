/*****************************************************************************************/
/** VTK Converter for CFD Solution Files */
/** Converts .dat format solution files to VTK format for visualization */
/** Author: AI Assistant for Shrinath K S */
/*****************************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <iomanip>

using namespace std;

struct GridPoint
{
    double x, y;
    double rho, mach, pressure, u, v;
    bool valid; // true if this is a fluid point
};

void convertDatToVTK(const string &inputFile, const string &outputFile)
{
    ifstream inFile(inputFile);
    if (!inFile.is_open())
    {
        cout << "ERROR: Could not open input file: " << inputFile << endl;
        return;
    }

    vector<GridPoint> points;
    string line;
    int nx = 0, ny = 0;

    // Skip header lines
    getline(inFile, line); // TITLE line
    getline(inFile, line); // VARIABLES line
    getline(inFile, line); // Zone line

    // Extract grid dimensions from Zone line
    // Format: "Zone T = Omega I = 120 J = 40"
    size_t pos_i = line.find("I = ");
    size_t pos_j = line.find("J = ");
    if (pos_i != string::npos && pos_j != string::npos)
    {
        string nx_str = line.substr(pos_i + 4);
        string ny_str = line.substr(pos_j + 4);
        nx = stoi(nx_str.substr(0, nx_str.find(' ')));
        ny = stoi(ny_str);
    }

    cout << "Grid dimensions detected: " << nx << " x " << ny << endl;

    // Read data points
    double x, y, rho, mach, pressure, u, v;
    while (inFile >> x >> y >> rho >> mach >> pressure >> u >> v)
    {
        GridPoint point;
        point.x = x;
        point.y = y;
        point.rho = rho;
        point.mach = mach;
        point.pressure = pressure;
        point.u = u;
        point.v = v;
        point.valid = (rho > 0.0); // Check if this is a valid fluid point
        points.push_back(point);
    }
    inFile.close();

    cout << "Read " << points.size() << " data points" << endl;

    // Write VTK file
    ofstream vtkFile(outputFile);
    if (!vtkFile.is_open())
    {
        cout << "ERROR: Could not create output file: " << outputFile << endl;
        return;
    }

    // VTK Header
    vtkFile << "# vtk DataFile Version 3.0" << endl;
    vtkFile << "2D CFD Solution - Converted from " << inputFile << endl;
    vtkFile << "ASCII" << endl;
    vtkFile << "DATASET STRUCTURED_GRID" << endl;

    // Grid dimensions
    vtkFile << "DIMENSIONS " << nx << " " << ny << " 1" << endl;

    // Points
    vtkFile << "POINTS " << points.size() << " double" << endl;
    for (const auto &point : points)
    {
        vtkFile << point.x << " " << point.y << " 0.0" << endl;
    }

    // Point data
    vtkFile << "POINT_DATA " << points.size() << endl;

    // Density
    vtkFile << "SCALARS Density double 1" << endl;
    vtkFile << "LOOKUP_TABLE default" << endl;
    for (const auto &point : points)
    {
        vtkFile << point.rho << endl;
    }

    // Pressure
    vtkFile << "SCALARS Pressure double 1" << endl;
    vtkFile << "LOOKUP_TABLE default" << endl;
    for (const auto &point : points)
    {
        vtkFile << point.pressure << endl;
    }

    // Mach number
    vtkFile << "SCALARS Mach_Number double 1" << endl;
    vtkFile << "LOOKUP_TABLE default" << endl;
    for (const auto &point : points)
    {
        vtkFile << point.mach << endl;
    }

    // Velocity magnitude
    vtkFile << "SCALARS Velocity_Magnitude double 1" << endl;
    vtkFile << "LOOKUP_TABLE default" << endl;
    for (const auto &point : points)
    {
        double vel_mag = sqrt(point.u * point.u + point.v * point.v);
        vtkFile << vel_mag << endl;
    }

    // Velocity vectors
    vtkFile << "VECTORS Velocity double" << endl;
    for (const auto &point : points)
    {
        vtkFile << point.u << " " << point.v << " 0.0" << endl;
    }

    // Fluid region indicator
    vtkFile << "SCALARS Fluid_Region int 1" << endl;
    vtkFile << "LOOKUP_TABLE default" << endl;
    for (const auto &point : points)
    {
        vtkFile << (point.valid ? 1 : 0) << endl;
    }

    vtkFile.close();

    cout << "VTK file created successfully: " << outputFile << endl;
    cout << "Open in ParaView or VisIt for visualization" << endl;
    cout << endl;
    cout << "Visualization Tips:" << endl;
    cout << "1. Load the .vtk file in ParaView" << endl;
    cout << "2. Use 'Density' or 'Mach_Number' for contour plots" << endl;
    cout << "3. Use 'Velocity' vectors for flow visualization" << endl;
    cout << "4. Use 'Fluid_Region' to highlight the computational domain" << endl;
}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        cout << "Usage: " << argv[0] << " <input_file.dat> <output_file.vtk>" << endl;
        cout << "Example: " << argv[0] << " Aerospike_IBM_MPI_Solution_790x1185_t1500.dat solution.vtk" << endl;
        return 1;
    }

    string inputFile = argv[1];
    string outputFile = argv[2];

    cout << "Converting " << inputFile << " to " << outputFile << endl;

    convertDatToVTK(inputFile, outputFile);

    return 0;
}
