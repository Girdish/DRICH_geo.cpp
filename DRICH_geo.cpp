
// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022, 2023 Christopher Dilks, Junhuai Xu

//==========================================================================
//  dRICH: Dual Ring Imaging Cherenkov Detector
//--------------------------------------------------------------------------
//
// Author: Christopher Dilks (Duke University)
//
// - Design Adapted from Standalone Fun4all and GEMC implementations
//   [ Evaristo Cisbani, Cristiano Fanelli, Alessio Del Dotto, et al. ]
//
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"

#include <XML/Helper.h>

using namespace dd4hep;
using namespace dd4hep::rec;

// create the detector
static Ref_t createDetector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens) {

  xml::DetElement detElem       = handle;
  std::string detName           = detElem.nameStr();
  int detID                     = detElem.id();
  xml::Component dims           = detElem.dimensions();
  OpticalSurfaceManager surfMgr = desc.surfaceManager();
  DetElement det(detName, detID);
  sens.setType("tracker");

  // attributes, from compact file =============================================
  // - vessel
  auto vesselZmin      = dims.attr<double>(_Unicode(zmin));
  auto vesselLength    = dims.attr<double>(_Unicode(length));
  auto vesselRmin0     = dims.attr<double>(_Unicode(rmin0));
  auto vesselRmin1     = dims.attr<double>(_Unicode(rmin1));
  auto vesselRmax0     = dims.attr<double>(_Unicode(rmax0));
  auto vesselRmax1     = dims.attr<double>(_Unicode(rmax1));
  auto vesselRmax2     = dims.attr<double>(_Unicode(rmax2));
  auto snoutLength     = dims.attr<double>(_Unicode(snout_length));
  auto nSectors        = dims.attr<int>(_Unicode(nsectors));
  auto wallThickness   = dims.attr<double>(_Unicode(wall_thickness));
  auto windowThickness = dims.attr<double>(_Unicode(window_thickness));
  auto vesselMat       = desc.material(detElem.attr<std::string>(_Unicode(material)));
  auto gasvolMat       = desc.material(detElem.attr<std::string>(_Unicode(gas)));
  auto vesselVis       = desc.visAttributes(detElem.attr<std::string>(_Unicode(vis_vessel)));
  auto gasvolVis       = desc.visAttributes(detElem.attr<std::string>(_Unicode(vis_gas)));
  // - radiator (applies to aerogel)
  auto radiatorElem       = detElem.child(_Unicode(radiator));
  auto radiatorRmin       = radiatorElem.attr<double>(_Unicode(rmin));
  auto radiatorRmax       = radiatorElem.attr<double>(_Unicode(rmax));
  auto radiatorPitch      = radiatorElem.attr<double>(_Unicode(pitch));
  auto radiatorFrontplane = radiatorElem.attr<double>(_Unicode(frontplane));
  // - aerogel
  auto aerogelElem      = radiatorElem.child(_Unicode(aerogel));
  auto aerogelMat       = desc.material(aerogelElem.attr<std::string>(_Unicode(material)));
  auto aerogelVis       = desc.visAttributes(aerogelElem.attr<std::string>(_Unicode(vis)));
  auto aerogelThickness = aerogelElem.attr<double>(_Unicode(thickness));
  // - mirror
  auto mirrorElem      = detElem.child(_Unicode(mirror));
  auto mirrorMat       = desc.material(mirrorElem.attr<std::string>(_Unicode(material)));
  auto mirrorVis       = desc.visAttributes(mirrorElem.attr<std::string>(_Unicode(vis)));
  auto mirrorSurf      = surfMgr.opticalSurface(mirrorElem.attr<std::string>(_Unicode(surface)));
  auto mirrorBackplane = mirrorElem.attr<double>(_Unicode(backplane));
  auto mirrorThickness = mirrorElem.attr<double>(_Unicode(thickness));
  auto mirrorRmin      = mirrorElem.attr<double>(_Unicode(rmin));
  auto mirrorRmax      = mirrorElem.attr<double>(_Unicode(rmax));
  auto mirrorPhiw      = mirrorElem.attr<double>(_Unicode(phiw));
  auto focusTuneZ      = mirrorElem.attr<double>(_Unicode(focus_tune_z));
  auto focusTuneX      = mirrorElem.attr<double>(_Unicode(focus_tune_x));
  // - sensorboxes
  auto sensorboxLength = desc.constant<double>("DRICH_sensorbox_length");
  auto sensorboxRmin   = desc.constant<double>("DRICH_sensorbox_rmin");
  auto sensorboxRmax   = desc.constant<double>("DRICH_sensorbox_rmax");
  auto sensorboxDphi   = desc.constant<double>("DRICH_sensorbox_dphi");
  // - sensor photosensitive surface (pss)
  auto pssElem      = detElem.child(_Unicode(sensors)).child(_Unicode(pss));
  auto pssMat       = desc.material(pssElem.attr<std::string>(_Unicode(material)));
  auto pssVis       = desc.visAttributes(pssElem.attr<std::string>(_Unicode(vis)));
  auto pssSurf      = surfMgr.opticalSurface(pssElem.attr<std::string>(_Unicode(surface)));
  auto pssSide      = pssElem.attr<double>(_Unicode(side));
  auto pssThickness = pssElem.attr<double>(_Unicode(thickness));
  // - sensor resin
  auto resinElem      = detElem.child(_Unicode(sensors)).child(_Unicode(resin));
  auto resinMat       = desc.material(resinElem.attr<std::string>(_Unicode(material)));
  auto resinVis       = desc.visAttributes(resinElem.attr<std::string>(_Unicode(vis)));
  auto resinSide      = resinElem.attr<double>(_Unicode(side));
  auto resinThickness = resinElem.attr<double>(_Unicode(thickness));
  // - photodetector unit (PDU)
  auto pduElem       = detElem.child(_Unicode(sensors)).child(_Unicode(pdu));
  auto pduNumSensors = desc.constant<int>("DRICH_pdu_num_sensors");
  auto pduSensorGap  = desc.constant<double>("DRICH_pdu_sensor_gap");
  auto pduGap        = desc.constant<double>("DRICH_pdu_gap");
  // - sensor sphere
  auto sensorSphElem    = detElem.child(_Unicode(sensors)).child(_Unicode(sphere));
  auto sensorSphRadius  = sensorSphElem.attr<double>(_Unicode(radius));
  auto sensorSphCenterX = sensorSphElem.attr<double>(_Unicode(centerx));
  auto sensorSphCenterZ = sensorSphElem.attr<double>(_Unicode(centerz));
  // - sensor sphere patch cuts
  auto sensorSphPatchElem = detElem.child(_Unicode(sensors)).child(_Unicode(sphericalpatch));
  auto sensorSphPatchPhiw = sensorSphPatchElem.attr<double>(_Unicode(phiw));
  auto sensorSphPatchRmin = sensorSphPatchElem.attr<double>(_Unicode(rmin));
  auto sensorSphPatchRmax = sensorSphPatchElem.attr<double>(_Unicode(rmax));
  auto sensorSphPatchZmin = sensorSphPatchElem.attr<double>(_Unicode(zmin));
  // - sensor readout
  auto readoutName = detElem.attr<std::string>(_Unicode(readout));
  // - settings and switches
  auto debugOpticsMode = desc.constant<int>("DRICH_debug_optics");
  bool debugSector     = desc.constant<int>("DRICH_debug_sector") == 1;
  bool debugMirror     = desc.constant<int>("DRICH_debug_mirror") == 1;
  bool debugSensors    = desc.constant<int>("DRICH_debug_sensors") == 1;

  // if debugging optics, override some settings
  bool debugOptics = debugOpticsMode > 0;
  if (debugOptics) {
    printout(WARNING, "DRICH_geo", "DEBUGGING DRICH OPTICS");
    switch (debugOpticsMode) {
    case 1:
      vesselMat = aerogelMat = pssMat = gasvolMat = desc.material("VacuumOptical");
      break;
    case 2:
      vesselMat = aerogelMat = pssMat = desc.material("VacuumOptical");
      break;
    case 3:
      vesselMat = aerogelMat = gasvolMat = desc.material("VacuumOptical");
      break;
    default:
      printout(FATAL, "DRICH_geo", "UNKNOWN debugOpticsMode");
      return det;
    }
  }

  // if debugging anything, draw only one sector and adjust visibility
  if (debugOptics || debugMirror || debugSensors)
    debugSector = true;
  if (debugSector)
    gasvolVis = vesselVis = desc.invisible();

  // readout coder <-> unique sensor ID
  std::vector<std::string> sensorIDfields = {"pdu", "sipm", "sector"};
  const auto& readoutCoder                = *desc.readout(readoutName).idSpec().decoder();
  uint64_t cellMask = 0;
  for (const auto& idField : sensorIDfields)
    cellMask |= readoutCoder[idField].mask();
  desc.add(Constant("DRICH_cell_mask", std::to_string(cellMask)));
  auto encodeSensorID = [&readoutCoder](auto ids) {
    uint64_t enc = 0;
    for (const auto& [idField, idValue] : ids)
      enc |= uint64_t(idValue) << readoutCoder[idField].offset();
    return enc;
  };

  // BUILD VESSEL ====================================================================
  double tankLength = vesselLength - snoutLength;
  double vesselZmax = vesselZmin + vesselLength;

  // snout solids
  double boreDelta  = vesselRmin1 - vesselRmin0;
  double snoutDelta = vesselRmax1 - vesselRmax0;
  Cone vesselSnout(snoutLength / 2.0, vesselRmin0, vesselRmax0,
                   vesselRmin0 + boreDelta * snoutLength / vesselLength, vesselRmax1);
  Cone gasvolSnout(
      snoutLength / 2.0, vesselRmin0 + wallThickness, vesselRmax0 - wallThickness,
      vesselRmin0 + boreDelta * (snoutLength - windowThickness) / vesselLength + wallThickness,
      vesselRmax1 - wallThickness + windowThickness * (vesselRmax1 - vesselRmax0) / snoutLength);

  // tank solids
  Polycone vesselTank(
      0, 2 * M_PI,
      {vesselSnout.rMin2(),
       std::lerp(vesselSnout.rMin2(), vesselRmin1, (sensorboxLength - snoutLength) / tankLength),
       vesselRmin1},
      {vesselSnout.rMax2(), vesselRmax2, vesselRmax2},
      {-tankLength / 2.0, -tankLength / 2.0 + sensorboxLength - snoutLength, tankLength / 2.0});
  Polycone gasvolTank(
      0, 2 * M_PI,
      {gasvolSnout.rMin2(),
       std::lerp(gasvolSnout.rMin2(), vesselRmin1 + wallThickness,
                 (sensorboxLength - snoutLength) / tankLength),
       vesselRmin1 + wallThickness},
      {gasvolSnout.rMax2(), vesselRmax2 - wallThickness, vesselRmax2 - wallThickness},
      {-tankLength / 2.0 + windowThickness,
       -tankLength / 2.0 + windowThickness + sensorboxLength - snoutLength,
       tankLength / 2.0 - windowThickness});

  // sensorbox solids
  double dphi = atan2(wallThickness, sensorboxRmax);
  Tube vesselSensorboxTube(sensorboxRmin, sensorboxRmax, sensorboxLength / 2., -sensorboxDphi / 2.,
                           sensorboxDphi / 2.);
  Tube gasvolSensorboxTube(sensorboxRmin + wallThickness, sensorboxRmax - wallThickness,
                           sensorboxLength / 2., -sensorboxDphi / 2. + dphi,
                           sensorboxDphi / 2. - dphi);

  // union: snout + tank
  UnionSolid vesselUnion(vesselTank, vesselSnout, Position(0., 0., -vesselLength / 2.));
  UnionSolid gasvolUnion(gasvolTank, gasvolSnout,
                         Position(0., 0., -vesselLength / 2. + windowThickness));

  // union: add sensorboxes for all sectors
  for (int isec = 0; isec < nSectors; isec++) {
    RotationZ sectorRotation((isec + 0.5) * 2 * M_PI / nSectors);
    vesselUnion = UnionSolid(
        vesselUnion, vesselSensorboxTube,
        Transform3D(sectorRotation, Position(0., 0., -(snoutLength + sensorboxLength - 0.6) / 2.)));
    gasvolUnion = UnionSolid(
        gasvolUnion, gasvolSensorboxTube,
        Transform3D(sectorRotation,
                    Position(0., 0., -(snoutLength + sensorboxLength) / 2. + windowThickness)));
  }

  // extra solids for `debugOptics` only
  Box vesselBox(1001, 1001, 1001);
  Box gasvolBox(1000, 1000, 1000);

  // choose vessel and gasvol solids
  Solid vesselSolid, gasvolSolid;
  switch (debugOpticsMode) {
  case 0:
    vesselSolid = vesselUnion;
    gasvolSolid = gasvolUnion;
    break;
  case 1:
  case 3:
    vesselSolid = vesselBox;
    gasvolSolid = gasvolBox;
    break;
  case 2:
    vesselSolid = vesselBox;
    gasvolSolid = gasvolUnion;
    break;
  }

  // volumes
  Volume vesselVol(detName, vesselSolid, vesselMat);
  Volume gasvolVol(detName + "_gas", gasvolSolid, gasvolMat);
  vesselVol.setVisAttributes(vesselVis);
  gasvolVol.setVisAttributes(gasvolVis);

  // reference positions
  auto originFront = Position(0., 0., -tankLength / 2.0 - snoutLength);
  auto vesselPos = Position(0, 0, vesselZmin) - originFront;

  // place gas volume
  PlacedVolume gasvolPV = vesselVol.placeVolume(gasvolVol, Position(0, 0, 0));
  DetElement gasvolDE(det, "gasvol_de", 0);
  gasvolDE.setPlacement(gasvolPV);

  // place mother volume (vessel)
  Volume motherVol      = desc.pickMotherVolume(det);
  PlacedVolume vesselPV = motherVol.placeVolume(vesselVol, vesselPos);
  vesselPV.addPhysVolID("system", detID);
  det.setPlacement(vesselPV);
  
  
  
 ////////////////////////////////////////////////aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa


// BUILD RADIATOR ====================================================================
// solid and volume: create aerogel with rims
Cone aerogelSolid(aerogelThickness / 2, radiatorRmin, radiatorRmax,
                  radiatorRmin + boreDelta * aerogelThickness / vesselLength,
                  radiatorRmax + snoutDelta * aerogelThickness / snoutLength);

// rim structure parameters
float rimThickness = 0.5; // thickness of the rim in z-direction (unchanged)
float rimWidth = 5.0;     // increased radial width of the rim (adjustable)

// Define rim solids for inner (minimum) and outer (maximum) edges with increased width
Tube rimInnerSolid(radiatorRmin - rimWidth, radiatorRmin, aerogelThickness / 2. + rimThickness);
Tube rimOuterSolid(radiatorRmax, radiatorRmax + rimWidth, aerogelThickness / 2. + rimThickness);

// rib structure
float originalSideLength = 2 * (radiatorRmax + 0.5 * aerogelThickness);
int nTilesX = 10;
int nTilesY = 10;
float tileSize = originalSideLength / nTilesX;
float ribThickness = 0.5;
auto ribMat = desc.material("CarbonFiber_15percent");
//auto ribVis = desc.visAttributes("color: #FF0000; alpha: 0.8");
  auto filterElem      = radiatorElem.child(_Unicode(filter));
 //rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
 // auto filterMat       = desc.material(filterElem.attr<std::string>(_Unicode(material)));
  auto filterVis       = desc.visAttributes(filterElem.attr<std::string>(_Unicode(vis)));

// Define rib solids once, outside the loop
Box aerogelRibSolidX(ribThickness / 2, tileSize / 2, aerogelThickness / 2);
Box aerogelRibSolidY(tileSize / 2, ribThickness / 2, aerogelThickness / 2);

// Create a union of all ribs and rims
Solid ribUnionSolid = nullptr;
bool firstRib = true;

// Build the union of ribs
for (int i = 0; i < nTilesX; i++) {
    for (int j = 0; j < nTilesY; j++) {
        float xPos = -((nTilesX - 1) * tileSize / 2) + i * tileSize;
        float yPos = -((nTilesY - 1) * tileSize / 2) + j * tileSize;

        // Horizontal ribs - left edge
        auto aerogelRibTrXLeft = Transform3D(Translation3D(xPos - tileSize / 2, yPos, 0.) * RotationY(radiatorPitch));
        if (firstRib) {
            ribUnionSolid = aerogelRibSolidX;
            firstRib = false;
        } else {
            ribUnionSolid = UnionSolid(ribUnionSolid, aerogelRibSolidX, aerogelRibTrXLeft);
        }

        // Horizontal ribs - right edge
        auto aerogelRibTrXRight = Transform3D(Translation3D(xPos + tileSize / 2, yPos, 0.) * RotationY(radiatorPitch));
        ribUnionSolid = UnionSolid(ribUnionSolid, aerogelRibSolidX, aerogelRibTrXRight);

        // Vertical ribs - bottom edge
        auto aerogelRibTrYBottom = Transform3D(Translation3D(xPos, yPos - tileSize / 2, 0.) * RotationY(radiatorPitch));
        ribUnionSolid = UnionSolid(ribUnionSolid, aerogelRibSolidY, aerogelRibTrYBottom);

        // Vertical ribs - top edge
        auto aerogelRibTrYTop = Transform3D(Translation3D(xPos, yPos + tileSize / 2, 0.) * RotationY(radiatorPitch));
        ribUnionSolid = UnionSolid(ribUnionSolid, aerogelRibSolidY, aerogelRibTrYTop);
    }
}

// Add rims to the union
ribUnionSolid = UnionSolid(ribUnionSolid, rimInnerSolid, Transform3D(RotationY(radiatorPitch)));
ribUnionSolid = UnionSolid(ribUnionSolid, rimOuterSolid, Transform3D(RotationY(radiatorPitch)));

// Intersect the aerogel with the rib and rim union to keep only the overlapping parts
Solid intersectedAerogelSolid = IntersectionSolid(aerogelSolid, ribUnionSolid);

// Subtract the rib and rim union from the aerogel to get the remaining parts
Solid remainingAerogelSolid = SubtractionSolid(aerogelSolid, ribUnionSolid);

// Create volumes for both parts of the aerogel
Volume intersectedAerogelVol(detName + "_aerogel_intersected", intersectedAerogelSolid, ribMat);
intersectedAerogelVol.setVisAttributes(filterVis);

Volume remainingAerogelVol(detName + "_aerogel_remaining", remainingAerogelSolid, aerogelMat);
remainingAerogelVol.setVisAttributes(aerogelVis);

// aerogel placement (both parts)
auto radiatorPos = Position(0., 0., radiatorFrontplane + 0.5 * aerogelThickness) + originFront;
auto aerogelPlacement = Translation3D(radiatorPos) * RotationY(radiatorPitch);

// Place the intersected aerogel (ribs and rims)
auto intersectedAerogelPV = gasvolVol.placeVolume(intersectedAerogelVol, aerogelPlacement);
DetElement intersectedAerogelDE(det, "aerogel_intersected_de", 0);
intersectedAerogelDE.setPlacement(intersectedAerogelPV);

// Place the remaining aerogel
auto remainingAerogelPV = gasvolVol.placeVolume(remainingAerogelVol, aerogelPlacement);
DetElement remainingAerogelDE(det, "aerogel_remaining_de", 1);
remainingAerogelDE.setPlacement(remainingAerogelPV);

// radiator z-position (w.r.t. IP)
double aerogelZpos = vesselPos.z() + intersectedAerogelPV.position().z();
desc.add(Constant("DRICH_aerogel_zpos", std::to_string(aerogelZpos)));
desc.add(Constant("DRICH_aerogel_material", aerogelMat.ptr()->GetName(), "string"));
desc.add(Constant("DRICH_gasvol_material", gasvolVol.material().ptr()->GetName(), "string"));



////////////////////////////////////////////////aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa


  // SECTOR LOOP //////////////////////////////////////////////////////////////////////
  for (int isec = 0; isec < nSectors; isec++) {
    if (debugSector && isec != 0)
      continue;

    RotationZ sectorRotation((isec + 0.5) * 2 * M_PI / nSectors);
    std::string secName = "sec" + std::to_string(isec);

    // BUILD MIRRORS ====================================================================
    double zS = sensorSphCenterZ + vesselZmin;
    double xS = sensorSphCenterX;
    double b = vesselZmax - mirrorBackplane;
    double zF = zS + focusTuneZ;
    double xF = xS + focusTuneX;

    double mirrorCenterZ = b * zF / (2 * b - zF);
    double mirrorCenterX = b * xF / (2 * b - zF);
    double mirrorRadius  = b - mirrorCenterZ;
    mirrorCenterZ -= vesselZmin;

    double mirrorThetaRot = std::asin(mirrorCenterX / mirrorRadius);
    double mirrorTheta1   = mirrorThetaRot - std::asin((mirrorCenterX - mirrorRmin) / mirrorRadius);
    double mirrorTheta2   = mirrorThetaRot + std::asin((mirrorRmax - mirrorCenterX) / mirrorRadius);

    if (debugMirror) {
      mirrorTheta1 = 0;
      mirrorTheta2 = M_PI;
    }

    Sphere mirrorSolid1(mirrorRadius, mirrorRadius + mirrorThickness, mirrorTheta1, mirrorTheta2,
                        -40 * degree, 40 * degree);
    auto mirrorPos = Position(mirrorCenterX, 0., mirrorCenterZ) + originFront;
    auto mirrorPlacement = Translation3D(mirrorPos) * RotationY(-mirrorThetaRot);

    Tube pieSlice(0.01 * cm, vesselRmax2, tankLength / 2.0, -mirrorPhiw / 2.0, mirrorPhiw / 2.0);
    IntersectionSolid mirrorSolid2(pieSlice, mirrorSolid1, mirrorPlacement);

    Volume mirrorVol(detName + "_mirror_" + secName, mirrorSolid2, mirrorMat);
    mirrorVol.setVisAttributes(mirrorVis);
    auto mirrorSectorPlacement = Transform3D(sectorRotation);
    auto mirrorPV = gasvolVol.placeVolume(mirrorVol, mirrorSectorPlacement);

    DetElement mirrorDE(det, "mirror_de_" + secName, isec);
    mirrorDE.setPlacement(mirrorPV);
    SkinSurface mirrorSkin(desc, mirrorDE, "mirror_optical_surface_" + secName, mirrorSurf,
                           mirrorVol);
    mirrorSkin.isValid();

    auto mirrorFinalPlacement = mirrorSectorPlacement * mirrorPlacement;
    auto mirrorFinalCenter = vesselPos + mirrorFinalPlacement.Translation().Vect();
    desc.add(Constant("DRICH_mirror_center_x_" + secName, std::to_string(mirrorFinalCenter.x())));
    desc.add(Constant("DRICH_mirror_center_y_" + secName, std::to_string(mirrorFinalCenter.y())));
    desc.add(Constant("DRICH_mirror_center_z_" + secName, std::to_string(mirrorFinalCenter.z())));
    if (isec == 0)
      desc.add(Constant("DRICH_mirror_radius", std::to_string(mirrorRadius)));

    // BUILD SENSORS ====================================================================
    if (debugSensors) {
      pssSide = 2 * M_PI * sensorSphRadius / 64;
    }

    auto sensorSphPos = Position(sensorSphCenterX, 0., sensorSphCenterZ) + originFront;
    auto sensorSphFinalCenter = sectorRotation * Position(xS, 0.0, zS);
    desc.add(Constant("DRICH_sensor_sph_center_x_" + secName, std::to_string(sensorSphFinalCenter.x())));
    desc.add(Constant("DRICH_sensor_sph_center_y_" + secName, std::to_string(sensorSphFinalCenter.y())));
    desc.add(Constant("DRICH_sensor_sph_center_z_" + secName, std::to_string(sensorSphFinalCenter.z())));
    if (isec == 0)
      desc.add(Constant("DRICH_sensor_sph_radius", std::to_string(sensorSphRadius)));

    // SENSOR MODULE LOOP ------------------------
    int ipdu = 0;
    double pduPitch = pduNumSensors * resinSide + (pduNumSensors + 1) * pduSensorGap + pduGap;
    double nTheta = M_PI * sensorSphRadius / pduPitch;

    for (int t = 0; t < (int)(nTheta + 0.5); t++) {
      double thetaGen = t / ((double)nTheta) * M_PI;
      double nPhi = 2 * M_PI * sensorSphRadius * std::sin(thetaGen) / pduPitch;
      for (int p = 0; p < (int)(nPhi + 0.5); p++) {
        double phiGen = p / ((double)nPhi) * 2 * M_PI - M_PI;

        double xGen = sensorSphRadius * std::sin(thetaGen) * std::cos(phiGen);
        double yGen = sensorSphRadius * std::sin(thetaGen) * std::sin(phiGen);
        double zGen = sensorSphRadius * std::cos(thetaGen);
        double x = zGen;
        double y = xGen;
        double z = yGen;

        double zCheck = z + sensorSphCenterZ;
        double xCheck = x + sensorSphCenterX;
        double yCheck = y;
        double rCheck = std::hypot(xCheck, yCheck);
        double phiCheck = std::atan2(yCheck, xCheck);

        bool patchCut = std::fabs(phiCheck) < sensorSphPatchPhiw && zCheck > sensorSphPatchZmin &&
                        rCheck > sensorSphPatchRmin && rCheck < sensorSphPatchRmax;
        if (debugSensors)
          patchCut = std::fabs(phiCheck) < sensorSphPatchPhiw;
        if (patchCut) {
          Box pssSolid(pssSide / 2., pssSide / 2., pssThickness / 2.);
          Box resinSolid(resinSide / 2., resinSide / 2., resinThickness / 2.);
          SubtractionSolid resinSolidEmbedded(
              resinSolid, pssSolid,
              Transform3D(Translation3D(0., 0., (resinThickness - pssThickness) / 2.)));

          Volume pssVol(detName + "_pss_" + secName, pssSolid, pssMat);
          Volume resinVol(detName + "_resin_" + secName, resinSolidEmbedded, resinMat);
          pssVol.setVisAttributes(pssVis);
          resinVol.setVisAttributes(resinVis);

          if (!debugOptics || debugOpticsMode == 3)
            pssVol.setSensitiveDetector(sens);

          auto pduAssemblyPlacement =
            sectorRotation *
            Translation3D(sensorSphPos) *
            RotationX(phiGen) *
            RotationZ(thetaGen) *
            Translation3D(sensorSphRadius, 0., 0.) *
            RotationY(M_PI / 2) *
            RotationZ(-M_PI / 2);

          Assembly pduAssembly(detName + "_pdu_" + secName);
          double pduSensorPitch = resinSide + pduSensorGap;
          double pduSensorOffsetMax = pduSensorPitch * (pduNumSensors - 1) / 2.0;
          int isipm = 0;

          for (int sensorIx = 0; sensorIx < pduNumSensors; sensorIx++) {
            for (int sensorIy = 0; sensorIy < pduNumSensors; sensorIy++) {
              Assembly sensorAssembly(detName + "_sensor_" + secName);

              auto pssPlacement = Transform3D(Translation3D(0., 0., -pssThickness / 2.0));
              auto resinPlacement = Transform3D(Translation3D(0., 0., -resinThickness / 2.0));
              auto pduSensorOffsetX = sensorIx * pduSensorPitch - pduSensorOffsetMax;
              auto pduSensorOffsetY = sensorIy * pduSensorPitch - pduSensorOffsetMax;
              auto sensorAssemblyPlacement =
                  Transform3D(Translation3D(pduSensorOffsetX, pduSensorOffsetY, 0.0));

              auto pssPV = sensorAssembly.placeVolume(pssVol, pssPlacement);
              sensorAssembly.placeVolume(resinVol, resinPlacement);
              pduAssembly.placeVolume(sensorAssembly, sensorAssemblyPlacement);

              pssPV.addPhysVolID("sector", isec)
                   .addPhysVolID("pdu", ipdu)
                   .addPhysVolID("sipm", isipm);

              auto sensorID = encodeSensorID(pssPV.volIDs());
              std::string sensorIDname =
                  secName + "_pdu" + std::to_string(ipdu) + "_sipm" + std::to_string(isipm);
              DetElement pssDE(det, "sensor_de_" + sensorIDname, sensorID);
              pssDE.setPlacement(pssPV);

              if (!debugOptics || debugOpticsMode == 3) {
                SkinSurface pssSkin(desc, pssDE, "sensor_optical_surface_" + sensorIDname, pssSurf,
                                    pssVol);
                pssSkin.isValid();
              }

              auto pduOrigin = ROOT::Math::XYZPoint(0, 0, 0);
              auto sensorPos = Translation3D(vesselPos) * pduAssemblyPlacement * sensorAssemblyPlacement * pduOrigin;
              auto pduPos = Translation3D(vesselPos) * pduAssemblyPlacement * pduOrigin;

              auto normVector = [pduAssemblyPlacement](Direction n) {
                return pduAssemblyPlacement * n;
              };
              auto sensorNormX = normVector(Direction{1., 0., 0.});
              auto sensorNormY = normVector(Direction{0., 1., 0.});

              auto distActual = std::sqrt((sensorPos - sensorSphFinalCenter).Mag2());
              auto distExpected = std::hypot(pduSensorOffsetX, pduSensorOffsetY, sensorSphRadius);
              auto testOnSphere = distActual - distExpected;
              if (std::abs(testOnSphere) > 1e-6) {
                printout(ERROR, "DRICH_geo", "sensor %s failed on-sphere test; testOnSphere=%f",
                         sensorIDname.c_str(), testOnSphere);
                throw std::runtime_error("dRICH sensor position test failed");
              }

              Direction radialDir = Direction(pduPos) - sensorSphFinalCenter;
              auto sensorNormZ = sensorNormX.Cross(sensorNormY);
              auto testOrtho = sensorNormX.Dot(sensorNormY);
              auto testRadial = radialDir.Cross(sensorNormZ).Mag2();
              auto testDirection = radialDir.Dot(sensorNormZ);
              if (std::abs(testOrtho) > 1e-6 || std::abs(testRadial) > 1e-6 || testDirection <= 0) {
                printout(ERROR, "DRICH_geo", "sensor %s failed orientation test",
                         sensorIDname.c_str());
                printout(ERROR, "DRICH_geo", "  testOrtho     = %f; should be zero", testOrtho);
                printout(ERROR, "DRICH_geo", "  testRadial    = %f; should be zero", testRadial);
                printout(ERROR, "DRICH_geo", "  testDirection = %f; should be positive",
                         testDirection);
                throw std::runtime_error("dRICH sensor orientation test failed");
              }

              auto pssVarMap = pssDE.extension<VariantParameters>(false);
              if (pssVarMap == nullptr) {
                pssVarMap = new VariantParameters();
                pssDE.addExtension<VariantParameters>(pssVarMap);
              }
              auto addVecToMap = [pssVarMap](std::string key, auto vec) {
                pssVarMap->set<double>(key + "_x", vec.x());
                pssVarMap->set<double>(key + "_y", vec.y());
                pssVarMap->set<double>(key + "_z", vec.z());
              };
              addVecToMap("pos", sensorPos);
              addVecToMap("normX", sensorNormX);
              addVecToMap("normY", sensorNormY);
              printout(DEBUG, "DRICH_geo", "sensor %s:", sensorIDname.c_str());
              for (auto kv : pssVarMap->variantParameters)
                printout(DEBUG, "DRICH_geo", "    %s: %f", kv.first.c_str(),
                         pssVarMap->get<double>(kv.first));

              isipm++;
            }
          }

          Transform3D frontServiceTransformation = Transform3D(Translation3D(0., 0., -resinThickness));
          for (xml::Collection_t serviceElem(pduElem.child(_Unicode(frontservices)), _Unicode(service));
               serviceElem; ++serviceElem) {
            auto serviceName = serviceElem.attr<std::string>(_Unicode(name));
            auto serviceSide = serviceElem.attr<double>(_Unicode(side));
            auto serviceThickness = serviceElem.attr<double>(_Unicode(thickness));
            auto serviceMat = desc.material(serviceElem.attr<std::string>(_Unicode(material)));
            auto serviceVis = desc.visAttributes(serviceElem.attr<std::string>(_Unicode(vis)));
            Box serviceSolid(serviceSide / 2.0, serviceSide / 2.0, serviceThickness / 2.0);
            Volume serviceVol(detName + "_" + serviceName + "_" + secName, serviceSolid, serviceMat);
            serviceVol.setVisAttributes(serviceVis);
            frontServiceTransformation =
                Transform3D(Translation3D(0., 0., -serviceThickness / 2.0)) * frontServiceTransformation;
            pduAssembly.placeVolume(serviceVol, frontServiceTransformation);
            frontServiceTransformation =
                Transform3D(Translation3D(0., 0., -serviceThickness / 2.0)) * frontServiceTransformation;
          }

          auto boardsElem = pduElem.child(_Unicode(boards));
          auto boardsMat = desc.material(boardsElem.attr<std::string>(_Unicode(material)));
          auto boardsVis = desc.visAttributes(boardsElem.attr<std::string>(_Unicode(vis)));
          Transform3D backServiceTransformation;
          for (xml::Collection_t boardElem(boardsElem, _Unicode(board)); boardElem; ++boardElem) {
            auto boardName = boardElem.attr<std::string>(_Unicode(name));
            auto boardWidth = boardElem.attr<double>(_Unicode(width));
            auto boardLength = boardElem.attr<double>(_Unicode(length));
            auto boardThickness = boardElem.attr<double>(_Unicode(thickness));
            auto boardOffset = boardElem.attr<double>(_Unicode(offset));
            Box boardSolid(boardWidth / 2.0, boardThickness / 2.0, boardLength / 2.0);
            Volume boardVol(detName + "_" + boardName + "+" + secName, boardSolid, boardsMat);
            boardVol.setVisAttributes(boardsVis);
            auto boardTransformation =
                Translation3D(0., boardOffset, -boardLength / 2.0) * frontServiceTransformation;
            pduAssembly.placeVolume(boardVol, boardTransformation);
            if (boardName == "RDO")
              backServiceTransformation =
                  Translation3D(0., 0., -boardLength) * frontServiceTransformation;
          }

          for (xml::Collection_t serviceElem(pduElem.child(_Unicode(backservices)), _Unicode(service));
               serviceElem; ++serviceElem) {
            auto serviceName = serviceElem.attr<std::string>(_Unicode(name));
            auto serviceSide = serviceElem.attr<double>(_Unicode(side));
            auto serviceThickness = serviceElem.attr<double>(_Unicode(thickness));
            auto serviceMat = desc.material(serviceElem.attr<std::string>(_Unicode(material)));
            auto serviceVis = desc.visAttributes(serviceElem.attr<std::string>(_Unicode(vis)));
            Box serviceSolid(serviceSide / 2.0, serviceSide / 2.0, serviceThickness / 2.0);
            Volume serviceVol(detName + "_" + serviceName + "_" + secName, serviceSolid, serviceMat);
            serviceVol.setVisAttributes(serviceVis);
            backServiceTransformation =
                Transform3D(Translation3D(0., 0., -serviceThickness / 2.0)) * backServiceTransformation;
            pduAssembly.placeVolume(serviceVol, backServiceTransformation);
            backServiceTransformation =
                Transform3D(Translation3D(0., 0., -serviceThickness / 2.0)) * backServiceTransformation;
          }

          gasvolVol.placeVolume(pduAssembly, pduAssemblyPlacement);
          ipdu++;
        }
      }
    }

    if (isec == 0)
      desc.add(Constant("DRICH_num_pdus", std::to_string(ipdu)));
    else if (ipdu != desc.constant<int>("DRICH_num_pdus"))
      printout(WARNING, "DRICH_geo", "number of PDUs is not the same for each sector");
  }

  return det;
}

// clang-format off
DECLARE_DETELEMENT(epic_DRICH, createDetector)

