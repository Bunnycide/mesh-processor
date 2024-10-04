#include <pmp/surface_mesh.h>
#include <pmp/visualization/mesh_viewer.h>
#include "pmp/io/io.h"

#include <math.h>

#include "pmp/algorithms/differential_geometry.h"

#define ANGLE_360_DEGREES 360.0f

using namespace pmp;

// Energy function to compute vertex importance
// Here the energy is relevant to angular distortion and visibility
// TODO: Implement this
Scalar get_edge_energy()
{
  return  0;
}

// Calculate angular distortion at a vertex
std::pair<Vertex, float> calculate_angular_distortion_at_vertex( const SurfaceMesh mesh, const Vertex& vertex)
{
  return {};
}

// Return the angle between two edges as a pair of (angle in radians, angle in degrees)
std::pair<float, float> get_angle_between_two_halfedges(const SurfaceMesh& mesh, const Halfedge& h1, const Halfedge& h2)
{
  // Get the first halfedge's start and end vertices
  auto start_vertex_curr = mesh.from_vertex(h1);
  auto end_vertex_curr = mesh.to_vertex(h1);

  // Get the second halfedge's start and end vertices
  auto start_vertex_next = mesh.from_vertex(h2);
  auto end_vertex_next = mesh.to_vertex(h2);

  // Get the vectors representing both the halfedges
  Point vector_curr = mesh.position(end_vertex_curr) - mesh.position(start_vertex_curr);
  Point vector_next = mesh.position(end_vertex_next) - mesh.position(start_vertex_next);

  // Compute the dot product between the current and next halfedge vectors
  float dot_product = dot(vector_curr, vector_next);

  // Compute the magnitudes of the vectors
  float magnitude_curr = norm(vector_curr);
  float magnitude_next = norm(vector_next);

  // Compute the cosine of the angle
  float cos_theta = dot_product / (magnitude_curr * magnitude_next);

  // Clamp cos_theta to [-1, 1] to avoid floating point errors
  cos_theta = std::max(-1.0f, std::min(1.0f, cos_theta));

  // Compute the angle in radians
  float angle_radians = std::acos(cos_theta);

  // Convert the angle to degrees
  float angle_degrees = angle_radians * (180.0f / M_PI);

  return {angle_radians, angle_degrees};
}

// Returns a pair of (veretx id, angular distortion)
std::vector<std::pair<Vertex, float>> get_angular_distortion(const SurfaceMesh& mesh)
{
  std::vector<std::pair<Vertex, float>> distortion_values_per_vertex;

  float angle_sum = 0.0f;

  for(auto v : mesh.vertices())
  {
    angle_sum = 0.0f;
    for(auto he : mesh.halfedges(v))
    {
      // Get the next halfedge in the face
      auto next_he = mesh.next_halfedge(he);

      auto _ = get_angle_between_two_halfedges(mesh, he, next_he);

      // Sum up the angles in degrees
      angle_sum += _.second;

      // If angle sum is over 360 mesh is not locally developable and distortion has to be calculated
      if( angle_sum > ANGLE_360_DEGREES)
      {
        // TODO : Find the angular distortion of these vertices as in the technical spec doc
        std::cout << "sub mesh not locally developable " << angle_sum << std::endl;

      } else
      {
        distortion_values_per_vertex.push_back({
          v, 0.0f
        });
      }
    }
  }

  return {};
}

int main(){
  SurfaceMesh mesh;
  read(mesh, "/home/jaeb/Workspace/geom-virt-test/common-3d-test-models/data/beast.obj");

  auto _ = get_angular_distortion(mesh);

  return 0;
}
