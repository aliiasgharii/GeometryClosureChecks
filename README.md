# Geometry Closure (Watertightness) Checks
## Introduction 
The internal spatial consistency principles assure 3D geospatial data is correctly and accurately defined based on the standards in each application.
The spatial extent of legal space, to which rights, restrictions and responsibilities relate in a 3D digital cadastre needs to be accurately defined and geometrically closed; watertight. 
We have tried to engage several techniques to formulate watertight concept for 3D digital cadastre. The proposed method comprises of two level checks, namely Primitive and Advanced checks.

![Flowchart](https://github.com/aliiasgharii/GeometryClosureChecks/blob/master/Image/flowchart.jpg)
## Dataset
In this project, a real cadastral dataset is used as our case study. This dataset includes 50 legal spaces representing both private ownership and common property.To testify the algorithms, different scenarios were designed to showcase how the geometry closure of 3D digital data associated with various complex legal spaces of a real multi-storey building can be checked through the primitive and advanced checks.  
![Flowchart](https://github.com/aliiasgharii/GeometryClosureChecks/blob/master/Image/Scenarios.JPG)
## Requirements
This code was implemented with C++ programming language using Computational Geometrical Algorithms Library (CGAL). Other requirements are as follows: 
- Visual Studio 2019
- Vcpkg Library Manager
- QT (Required for drawing and visualisation)
- CMake

## Citation
```
Asghari A, Kalantari M, Rajabifard A. Advances in techniques to formulate the watertight concept for cadastre. Transactions in GIS. 2020;00:1–25. https://doi.org/10.1111/tgis.12695
```
## Author
Ali Asghari
- Email: <aasghari@student.unimelb.edu.au>
- Linkedin: [Here](https://www.linkedin.com/in/aliiasgharii/)

## Acknowledgments
This project incorporates code from the following repos:
- https://github.com/CGAL/cgal

