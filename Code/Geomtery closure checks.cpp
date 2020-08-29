// This code has been developed to check the geomtry closure of cadastral objects.
// The code includes primitive and advanced checks to check the spatial consistecny of 2-manifold cadastral objects in OFF format.
// Copyright (c) 2020 Ali Asghari

#include <iostream>
#include <string>
#include <CGAL/Cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/draw_polyhedron.h>
#include <CGAL/draw_Surface_mesh.h>
#include <vector>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <fstream>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <iostream>
#include <dirent.h>
#include <sys/types.h>
typedef double Real;
typedef CGAL::Cartesian < Real > Kernel0;
typedef CGAL::Filtered_kernel < Kernel0 > Kernel;
typedef CGAL::Polyhedron_3 < Kernel > Polyhedron;
typedef Kernel::Point_3 Point;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>             Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
//typedef boost::graph_traits<Mesh>::vertex_descriptor                 vertex_descriptor;
typedef Polyhedron::Halfedge_handle    Halfedge_handle;
typedef Polyhedron::Facet_handle       Facet_handle;
typedef Polyhedron::Vertex_handle      Vertex_handle;
//typedef boost::graph_traits<Mesh>::halfedge_descriptor               halfedge_descriptor;
//namespace NP = CGAL::parameters;
namespace PMP = CGAL::Polygon_mesh_processing;
using namespace std;

// A function for listing the absolute paths for each object in a folder
void list_filepaths(const char* path) {
    struct dirent* entry;
    DIR* dir = opendir(path);

    if (dir == NULL) {
        return;
    }
    while ((entry = readdir(dir)) != NULL) {
        string ss = entry->d_name;
        if (entry->d_type == DT_REG)
            if (ss[0] != '.') {

                //cout << "C:/dev/MyProg/Test/" << ss << endl;
                ofstream dataset_filepaths;
                dataset_filepaths.open("C:/dev/MyProg/Test/Result/dataset_filepaths.txt", std::ios_base::app);
                dataset_filepaths << "C:/dev/MyProg/Test/" << ss << endl;
            }


    }
    closedir(dir);
}
// A function for listing the names of objects in a folder
void list_filenames(const char* path) {
    struct dirent* entry;
    DIR* dir = opendir(path);
    if (dir == NULL) {
        return;
    }
    while ((entry = readdir(dir)) != NULL) {
        string ss = entry->d_name;
        if (entry->d_type == DT_REG)
            if (ss[0] != '.') {
                // cout << ss << endl;
                ofstream dataset_filenames;
                dataset_filenames.open("C:/dev/MyProg/Test/Result/dataset_filenames.txt", std::ios_base::app);
                dataset_filenames << entry->d_name << endl;
            }
    } closedir(dir);
}

int main(int argc, char* argv[])
{
    list_filepaths("C:/dev/MyProg/Test");
    list_filenames("C:/dev/MyProg/Test");
    std::ifstream myFile;
    std::string line;
    int lines;
    // Reading 3D objects in a dataset ----------------------------------------------------------------------------------
    myFile.open("C:/dev/MyProg/Test/Result/dataset_filepaths.txt");
    for (lines = 1; std::getline(myFile, line); lines++) {

        std::cout << "The object is being assessed: " << line.substr(19) << '\n';
        Mesh mesh0;
        ifstream input((argc > 1) ? argv[1] : line);

        // Triangulate faces for checking self-intersection -----------------------------------------------------------------
        if (!input || !(input >> mesh0) || !CGAL::is_triangle_mesh(mesh0))
        {

        // Confirm that all faces are triangles.

            CGAL::Polygon_mesh_processing::triangulate_faces(mesh0);
            for (face_descriptor fit : faces(mesh0))
                if (next(next(halfedge(fit, mesh0), mesh0), mesh0)
                    != prev(halfedge(fit, mesh0), mesh0))
                {
                    cerr << "Error: non-triangular face left in mesh." << endl;
                    cerr << "Not a valid input file." << endl;
                    return 1;
                }
        }

        // Primitive checks: Checking self-intersection --------------------------------------------------------------------------------------
        bool intersecting = PMP::does_self_intersect(mesh0,
            PMP::parameters::vertex_point_map(get(CGAL::vertex_point, mesh0)));

        vector<pair<face_descriptor, face_descriptor> > intersected_tris;

        PMP::self_intersections(mesh0, back_inserter(intersected_tris));

        // Primitive checks: Checking collinearity (First method) ----------------------------------------------------------------------------
        int notcollinear = 0;
        for (face_descriptor fit : faces(mesh0))
        {
            //cout << PMP::face_area(fit, mesh0) << endl;
            if (PMP::face_area(fit, mesh0) <= 0) {

                ++notcollinear;
            }

        }
        // Checking orientatin ---------------------------------------------------------------------------------------------

        //if (CGAL::is_closed(mesh0))
        //bool orinted = PMP::is_outward_oriented(mesh0);
        //---------------------------------------------------------------------------------------------

        Polyhedron mesh;
        Real InDegreeSum = 0;
        Real numFaces = 0;
        Real faceSum = 0;
        int num_border_edge = 0;
        int min_InDegree = -1;
        int max_InDegree = -1;
        string f;
        f = "The 3D object is Watertight.";

        // Read each object in our dataset in OFF format ----------------------------------------------------------------------------
        ifstream inp((argc > 1) ? argv[1] : line);
        if (!inp || !(inp >> mesh) || mesh.empty()) {
            std::cerr << "Not a valid dataset." << std::endl;
            return 1;
        }
        inp >> mesh;

        // Loop over all of the vertices in the mesh
        // Compute the minimum, maximum, and average InDegrees of the vertices.
        // For each vertex in the mesh...
        for (Polyhedron::Vertex_const_iterator vertexIter = mesh.vertices_begin();
            vertexIter != mesh.vertices_end(); ++vertexIter) {

            // Get the in-degree of the current vertex (i.e. The total number of incident edges of the current vertex).
            int InDegree = vertexIter->degree();

            // Update the minimum InDegree value.
            if (min_InDegree < 0 || InDegree < min_InDegree) {
                min_InDegree = InDegree;
            }
            // Update the maximum InDegree value.
            if (max_InDegree < 0 || InDegree > max_InDegree) {
                max_InDegree = InDegree;
            }
            InDegreeSum += InDegree;

            // Using a circulator that can be used to visit all of the
            // halfedges around the vertex in CW order.
            Polyhedron::Vertex::Halfedge_around_vertex_const_circulator
                halfEdgeCirc0 = vertexIter->vertex_begin();

            // Advanced checks: Checking the geometry closure (watertightness) of the legal spaces --------------------------------------------

            // ##################################################  Node_based method  #######################################################
            Real faceSum = 0;
            for (int k = 1; k <= InDegree; ++k) {
                const bool& c = halfEdgeCirc0->is_border();
                if (c == 0) {
                    numFaces = 1;
                    faceSum += numFaces;
                }

                ++halfEdgeCirc0;
            }
            faceSum;
            if (faceSum != InDegree)
            {
                f = "The 3D object is not Watertight!";
                //cout << "Number of incident faces for each vertex: " << faceSum << endl;
                //cout << "Number of incident edges for each vertex: " << InDegree << endl;
            } 
            // ##################################################   Edge_based method    ####################################################
           /*  int sum_num_border = 0;
            //int num_border = 0;
            for (int k = 1; k <= InDegree; ++k) {
                const bool c = halfEdgeCirc0->is_border();
                if (c == 1) {
                    // num_border = 1;
                    sum_num_border += 1;
                }
                ++halfEdgeCirc0;
            }
            sum_num_border;
            if (sum_num_border >= 1)
            {
                f = "The 3D object is not Watertight!";
            }*/
            //###############################################################################################################################
        }
        //Real meanInDegree = InDegreeSum / mesh.size_of_vertices();

        // Calculating the number of holes on the surface of the 3D object and repair the surface by filling the holes ----------------------------------------------------

        Polyhedron poly;
        poly = mesh;
        unsigned int nb_holes = 0;
        for (Halfedge_handle h : halfedges(poly))
        {
            if (h->is_border())
            {
                vector<Facet_handle>  patch_facets;
                vector<Vertex_handle> patch_vertices;
                bool success = get<0>(
                    PMP::triangulate_refine_and_fair_hole(
                        poly,
                        h,
                        back_inserter(patch_facets),
                        back_inserter(patch_vertices),
                        PMP::parameters::vertex_point_map(get(CGAL::vertex_point, poly)).
                        geom_traits(Kernel())));
                //std::cout << " Number of facets in constructed patch: " << patch_facets.size() << std::endl;
                //std::cout << " Number of vertices in constructed patch: " << patch_vertices.size() << std::endl;
                //std::cout << " Fairing : " << (success ? "succeeded" : "failed") << std::endl;
                ++nb_holes;
            }
        }
        /*ofstream out("C:/dev/MyProg/Test/repaired.off");
        out.precision(17);
        out << poly << endl;*/
        //-------------------------------------------------------------------------------------------------------------------------------------
        // Primitive checks: Check for collinearity and planarity ------------------------------------------------------------------------------------------------
        int numNonplanarFaces = 0;
        Real degreeSum = 0;
        int minDegree = -1;
        int maxDegree = -1;
        int triNum = 0;
        int triSum = 0;
        int quadNum = 0;
        int quadSum = 0;

        // For each face in the mesh...
        for (Polyhedron::Facet_const_iterator faceIter = mesh.facets_begin();
            faceIter != mesh.facets_end(); ++faceIter) {

            // Get the degree of the face.
            int degree = faceIter->facet_degree();

            // Using a circulator that can be used to visit all of the
            // halfedges around the face in CCW order (where the halfedge
            // has the face on its left side).
            Polyhedron::Facet::Halfedge_around_facet_const_circulator
                halfEdgeCirc = faceIter->facet_begin();

            // Get the first vertex of the face.
            const Point& v0 = halfEdgeCirc->vertex()->point();
            ++halfEdgeCirc;

            // Get the second vertex of the face.
            const Point& v1 = halfEdgeCirc->vertex()->point();
            ++halfEdgeCirc;

            // Get the third vertex of the face.
            const Point& v2 = halfEdgeCirc->vertex()->point();
            ++halfEdgeCirc;

            // Check collinearity of verticies for triangle faces (Second method). 
            if (degree == 3) {
                if (CGAL::collinear(v0, v1, v2))
                {
                    cout << " Collinearity detected : " << "\n"
                        << "(" << v0.x() << " ," << v0.y() << " ," << v0.z() << ") " << "\n"
                        << "(" << v1.x() << " ," << v1.y() << " ," << v1.z() << ") " << "\n"
                        << "(" << v2.x() << " ," << v2.y() << " ," << v2.z() << ") " << "\n";
                }
                //cout << (CGAL::collinear(v0, v1, v2) ? "No collinearity detected.\n" : "Collinearity detected.\n");
            }

            // Primitive checks: Checking planarity 
            // First method
            // The face can only be nonplanar if its degree exceeds three.
            if (degree >= 4) {
                int notCoplanar = 0;
                // Check that each remaining vertex is coplanar with the first three vertices.
                for (int i = 3; i < degree; ++i) {
                    // Get the next vertex of the face.
                    const Point& v = halfEdgeCirc->vertex()->point();
                    ++halfEdgeCirc;
                    // Check if the vertex is coplanar with the first three.
                    if (!CGAL::coplanar(v0, v1, v2, v)) {
                        ++notCoplanar;
                        if (notCoplanar == 1) {
                            cout << " nonplanar face detected : " << "\n"
                                << "(" << v0.x() << " ," << v0.y() << " ," << v0.z() << ") " << "\n"
                                << "(" << v1.x() << " ," << v1.y() << " ," << v1.z() << ") " << "\n"
                                << "(" << v2.x() << " ," << v2.y() << " ," << v2.z() << ") " << "\n";
                        }
                        cout
                            << "(" << v.x() << " ," << v.y() << " ," << v.z() << ") " << "\n";
                    }
                }
                if (notCoplanar >= 1) {
                    ++numNonplanarFaces;
                }
            }

            // Update the minimum degree.
            if (minDegree < 0 || degree < minDegree) {
                minDegree = degree;
            }
            // Update the maximum degree.
            if (maxDegree < 0 || degree > maxDegree) {
                maxDegree = degree;
            }
            degreeSum += degree;

            // Calculate the number of triangles, quads and faces
            bool tri = faceIter->is_triangle();
            bool qua = faceIter->is_quad();
            if (tri == 1) {
                triNum = 1;
                triSum += triNum;
            }
            if (qua == 1) {
                quadNum = 1;
                quadSum += quadNum;
            }
        }

        int E = mesh.size_of_halfedges() / 2;
        int F = mesh.size_of_facets();
        int T = triSum;
        int V = mesh.size_of_vertices();
        int Q = quadSum;
        int arbitraryFaces = F - T - Q;

        // Determine mesh type-----------------------------------
        string meshType;
        if (mesh.is_pure_triangle()) {
            meshType = string("triangles.");
        }
        else if (mesh.is_pure_quad()) {
            meshType = string("quads.");
        }
        else if (T >= 1 && Q >= 1) {
            meshType = string("triangles and quads.");
        }
        else {
            meshType = string("arbitrary convex faces.");
        }
        // ########################################  Determine mesh type and applying Euler Formula    #########################################
       /* string meshType;
        if (mesh.is_pure_triangle()) {
            meshType = std::string("triangles.");
            if (2 * E == 3 * T && T == 2 * V - 4 && E == 3 * V - 6)
            {
                f = "The 3D object is Watertight!";
            }
            else { f = "The 3D object is not Watertight!"; }
        }

        else if (mesh.is_pure_quad()) {
            meshType = std::string("quads.");
            if (2 * E == 4 * Q && Q == V - 2 && E == 2 * V - 4)
            {
                f = "The 3D object is Watertight.";
            }
            else { f = "The 3D object is not Watertight!"; }
        }

        else if (T >= 1 && Q >= 1) {
            meshType = std::string("triangles and quads.");
            if ((2 * E) == (3 * T) + (4 * Q) && T == (2 * V) - (2 * Q) - 4)
            {
                f = "The 3D object is Watertight!";
            }
            else { f = "The 3D object is not Watertight!"; }
        }

        else {
            meshType = std::string("arbitrary convex faces.");
            if (2 * E >= 3 * F && F <= 2 * V - 4 && E <= 3 * V - 6)
            {
                f = "The special 3D object is Watertight!";
            }
            else { f = "The special 3D object is not Watertight!"; }
        }*/

        // ######################################################################################################################################
        Real meanDegree = degreeSum / mesh.size_of_facets();
        //-------------------------------------------------------------------
        mesh.normalize_border();
        //---------------------------------------------------------------
        // Output the mesh information.
        ofstream file;
        file.open("C:/dev/MyProg/Test/Result/Validation Report.txt", std::ios_base::app);
        if (!file)
        {
            cout << "Error in creating file!!!";

        }
        file
            
            << " Internal Spatial Consistency Report... " << "\n" << "\n"
            << " Data Name: " << line.substr(19) << "\n"
            << " Data Type:" << " 3D Polyhedral model consisted of " << meshType << "\n" << "\n"
            
            << " Checking vertices... "  "\n"
            << " Number of vertices: " << mesh.size_of_vertices() << "\n"
            << " Minimum and Maximum vertex in-degree: (" << min_InDegree << "," << max_InDegree << ")" << "\n" << "\n"

            << " Checking edges... "  "\n"
            << " Number of edges and halfedges: (" << mesh.size_of_halfedges() / 2 << "," << mesh.size_of_halfedges() << ")" << "\n"
            << " Number of border edges: " << mesh.size_of_border_edges() << "\n" << "\n"

            << " Checking faces... "  "\n"
            << " A total number of "
            << F
            << " faces including: "   "\n"
            <<" "<< T << " trinagle(s), "
            << Q << " quadrat(s), and "
            << arbitraryFaces << " arbitrary face(s)." << "\n"
            << " Minimum and Maximum face degree: (" << minDegree << "," << maxDegree << ")" << "\n" << "\n"

            << " Primitive Checks... "  "\n"
            << " Number of nonplanar faces: " << numNonplanarFaces << "\n"
            << " Self-intersection status:" << (intersecting ? " The 3D object is self-intersected." : " The 3D object is not self-intersected.") << "\n"
            << " Pairs of faces intersected:" << intersected_tris.size() << "\n"
            << " Number of collinearity detected: " << notcollinear << "\n"
            << " Orientation status: " << (PMP::is_outward_oriented(mesh0) ? " The 3D object is correctly outward oriented." : "The 3D object is not correctly oriented.") << "\n" << "\n"

            << " Advanced Check...  "  "\n"
            << " Closure status: " << f << "\n"
            << " Number of holes detected: " << nb_holes << "\n"

            << " ---------------------------------------------------------------- " << "\n" 
            //<< " The Number of holes have been filled: " << nb_holes << "\n" << "\n"
            ;

        //CGAL::draw(mesh0);
        CGAL::draw(mesh);
    };
    return 0;
}
