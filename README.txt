/** \mainpage
@brief
G4QSim --- Geant4 Quick Simulation
\verbatim

G4QSim stands for ``Geant4 Quick Simulation''.  It accepts text files
(gdml and text) for geometry description, and .mac files for gun,
source (/gps/) and run control.  The basic idea is not to modify the
source code while being able to simulate small setups as much and as
quickly as possible.

By setting `/G4QSim/event/storeTrajectory true', the program reuses
the same ROOT file structure but saves trajectory information into the
data file.  See `AnalysisMessenger' for details.

Haven't figured out ``TrackingAction'' from ``rdecay'' example yet, so
it is practically disabled (causes crash if enabled).  He6 and Sr90
decays seem to work without it.

################################################################################
See macro/dbdXe136.mac and geometry/Xe136HPG.gdml for an example to
generate Xe-136 double-beta decay events:

./G4QSim -g ../geometry/Xe136HPG.gdml -u 1 -m ../macro/dbdXe136.mac

\endverbatim
*/
