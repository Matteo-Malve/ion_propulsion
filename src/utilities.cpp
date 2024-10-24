#include "utilities.h"

std::string extract_mesh_name() {
  // Find the last '/' character to isolate the filename
  size_t lastSlash = PATH_TO_MESH.rfind('/');
  // Find the last '.' character to remove the file extension
  size_t lastDot = PATH_TO_MESH.rfind('.');
  // Extract the substring between the last '/' and the last '.'
  std::string meshName = PATH_TO_MESH.substr(lastSlash + 1, lastDot - lastSlash - 1);
  return meshName;
}