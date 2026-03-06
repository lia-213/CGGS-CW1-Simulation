#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>
#include "readOFF.h"
#include "scene.h"
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <set>
#include <array>
#include <chrono>

using namespace Eigen;
using namespace std;

bool isAnimating = false;

polyscope::SurfaceMesh* pMesh;
polyscope::CurveNetwork* pConstraints;

double currTime = 0;
double timeStep = 0.02;
double CRCoeff = 1.0;
double tolerance = 1e-3;
int maxIterations = 10000;
int screenshotCounter = 0;
bool pauseOnFirstContact = true;
bool firstContactReached = false;
bool showTrails = true;

// COM position history for trajectory trails
vector<vector<RowVector3d>> comTrails;
polyscope::CurveNetwork* pTrails = nullptr;

Scene scene;

void callback_function() {
  ImGui::PushItemWidth(50);
  ImGui::TextUnformatted("Animation Parameters");
  ImGui::Separator();
  ImGui::Checkbox("isAnimating", &isAnimating);
  ImGui::Checkbox("Naive O(n^2) broad-phase", &scene.useNaive);
  ImGui::PopItemWidth();
  ImGui::SetNextItemWidth(200);
  ImGui::InputInt("maxIterations", &maxIterations);
  ImGui::SetNextItemWidth(200);
  ImGui::SliderFloat("timeStep", (float*)&timeStep, 0.001f, 0.5f);

  if (ImGui::Button("Screenshot")) {
    string fname = "screenshot_" + to_string(screenshotCounter++) + ".png";
    polyscope::screenshot(fname, true);
    cout << "Saved " << fname << "\n";
  }
  ImGui::Checkbox("Pause on first contact", &pauseOnFirstContact);
  ImGui::Checkbox("Show trails", &showTrails);
  if (showTrails && ImGui::Button("Clear trails")) {
    for (auto& t : comTrails) t.clear();
    firstContactReached = false;
  }
  if (firstContactReached)
    ImGui::TextColored(ImVec4(1,0.3f,0.3f,1), "FIRST CONTACT");

  if (!isAnimating)
    return;

  scene.update_scene(timeStep, CRCoeff, maxIterations, tolerance);
  pMesh->updateVertexPositions(scene.currV);
  pConstraints->updateNodePositions(scene.currConstVertices);

  // Detect first floor contact: colour touching meshes red
  Eigen::MatrixXd vertexColors = Eigen::MatrixXd::Ones(scene.currV.rows(), 3) * 0.7;
  int offset = 0;
  for (int i = 0; i < (int)scene.meshes.size(); i++) {
    double minY = scene.meshes[i].currV.col(1).minCoeff();
    bool touching = minY <= 0.01;
    if (touching && !firstContactReached) {
      firstContactReached = true;
      if (pauseOnFirstContact) isAnimating = false;
      cout << "First floor contact: mesh " << i << " at t=" << scene.currTime << "\n";
    }
    Eigen::RowVector3d col = touching ? Eigen::RowVector3d(1.0, 0.2, 0.2)
                                      : Eigen::RowVector3d(0.6, 0.8, 1.0);
    for (int j = 0; j < scene.meshes[i].currV.rows(); j++)
      vertexColors.row(offset + j) = col;
    offset += scene.meshes[i].currV.rows();
  }
  pMesh->addVertexColorQuantity("contact", vertexColors)->setEnabled(true);

  // Record trails and redraw
  if (showTrails) {
    for (int i = 0; i < (int)scene.meshes.size(); i++)
      comTrails[i].push_back(scene.meshes[i].COM);

    // Build a single curve network from all trails
    int totalPts = 0;
    for (auto& t : comTrails) totalPts += (int)t.size();

    if (totalPts > (int)comTrails.size()) { // at least 2 points somewhere
      Eigen::MatrixXd trailV(totalPts, 3);
      vector<std::array<int,2>> trailE;
      int idx = 0;
      for (auto& t : comTrails) {
        for (int k = 0; k < (int)t.size(); k++) {
          trailV.row(idx + k) = t[k];
          if (k > 0) trailE.push_back({idx + k - 1, idx + k});
        }
        idx += (int)t.size();
      }
      Eigen::MatrixXi trailEMat(trailE.size(), 2);
      for (int k = 0; k < (int)trailE.size(); k++)
        trailEMat.row(k) << trailE[k][0], trailE[k][1];
      pTrails = polyscope::registerCurveNetwork("Trails", trailV, trailEMat);
      pTrails->setRadius(0.002);
      pTrails->setColor({1.0, 0.8, 0.0}); // yellow trails
    }
  }
}


int main()
{
  
  scene.load_scene("demo-spin-scene.txt","no-constraints.txt");

  // Give each object a different lateral velocity so paths trace visible parabolas - specifically for fig1 with mixed-scene.txt
  // for (int i = 0; i < (int)scene.meshes.size(); i++) {
  //   double vx = (i % 2 == 0 ? 1.0 : -1.0) * (3.0 + i * 1.5);
  //   double vz = (i % 3 == 0 ? 1.0 : -1.0) * 2.0;
  //   scene.meshes[i].comVelocity = RowVector3d(vx, 0.0, vz);
  // }
  comTrails.resize(scene.meshes.size());

  scene.meshes[0].comVelocity = RowVector3d(15.0, 0.0, 0.0);

  polyscope::init();
  scene.update_scene(0.0, CRCoeff, maxIterations, tolerance);

  // Visualization
  pMesh = polyscope::registerSurfaceMesh("Entire Scene", scene.currV, scene.allF);
  pConstraints = polyscope::registerCurveNetwork("Constraints", scene.currConstVertices, scene.constEdges);
  polyscope::options::groundPlaneHeightMode = polyscope::GroundPlaneHeightMode::Manual;
  polyscope::options::groundPlaneHeight = 0.;
  polyscope::state::userCallback = callback_function;

  polyscope::show();
  
}

