# SCRIPT TO RUN SABINO DOMAIN WITH TERRAIN-FOLLOWING GRID AND VARIABLE DZ
# DETAILS:
# Arugments are 1) runname 2) year

# Import the ParFlow TCL package
lappend   auto_path $env(PARFLOW_DIR)/bin
package   require parflow
namespace import Parflow::*

pfset     FileVersion    4

#-----------------------------------------------------------------------------
# Set Processor topology 
#-----------------------------------------------------------------------------
pfset Process.Topology.P 1
pfset Process.Topology.Q 1
pfset Process.Topology.R 1

#-----------------------------------------------------------------------------
# Make a directory for the simulation and copy inputs into it
#-----------------------------------------------------------------------------
exec mkdir "Outputs_spinup_CLMoff_OFoff"
cd "./Outputs_spinup_CLMoff_OFoff"

# ParFlow Inputs
file copy -force "../../parflow_inputs/tucson.slopex.pfb" .
file copy -force "../../parflow_inputs/tucson.slopey.pfb" .
file copy -force "../../parflow_inputs/geology_indicator.pfb"   .
file copy -force "../../parflow_inputs/eff_recharge_0013.pfb"  .
# file copy -force "../../parflow_input/press.init.nc"  .
# 
# CLM Inputs
# file copy -force "../../clm_input/drv_clmin.dat" .
# file copy -force "../../clm_input/drv_vegp.dat"  .
# file copy -force "../../clm_input/drv_vegm.dat"  . 
# file copy -force "../../clm_input/metForcing.nc"  . 

 puts "Files Copied"

#-----------------------------------------------------------------------------
# Computational Grid
#-----------------------------------------------------------------------------
pfset ComputationalGrid.Lower.X           0.0
pfset ComputationalGrid.Lower.Y           0.0
pfset ComputationalGrid.Lower.Z           0.0 

pfset ComputationalGrid.DX                90.0
pfset ComputationalGrid.DY                90.0
pfset ComputationalGrid.DZ                100.0

pfset ComputationalGrid.NX                246 
pfset ComputationalGrid.NY                178 
pfset ComputationalGrid.NZ                17  


#-----------------------------------------------------------------------------
# Names of the GeomInputs
#-----------------------------------------------------------------------------
pfset GeomInput.Names                     "box_input indi_input"

#-----------------------------------------------------------------------------
# Domain Geometry Input
#-----------------------------------------------------------------------------
pfset GeomInput.box_input.InputType      Box
pfset GeomInput.box_input.GeomName      domain

#-----------------------------------------------------------------------------
# Domain Geometry 
#-----------------------------------------------------------------------------
pfset Geom.domain.Lower.X                        0.0
pfset Geom.domain.Lower.Y                        0.0
pfset Geom.domain.Lower.Z                        0.0
 
pfset Geom.domain.Upper.X                        22140.0
pfset Geom.domain.Upper.Y                        16020.0
pfset Geom.domain.Upper.Z                         1700.0
pfset Geom.domain.Patches             "x-lower x-upper y-lower y-upper z-lower z-upper"

#-----------------------------------------------------------------------------
# Indicator Geometry Input
#-----------------------------------------------------------------------------
 pfset GeomInput.indi_input.InputType      IndicatorField
 pfset GeomInput.indi_input.GeomNames      "s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12 s13 s14 s15 s16 g1 g2 g3 g4"
 pfset Geom.indi_input.FileName            "geology_indicator.pfb"
 
 pfset GeomInput.s1.Value                3
 pfset GeomInput.s2.Value                5
 pfset GeomInput.s3.Value                13
 pfset GeomInput.s4.Value                30
 pfset GeomInput.s5.Value                35
 pfset GeomInput.s6.Value                36
 pfset GeomInput.s7.Value                40
 pfset GeomInput.s8.Value                47
 pfset GeomInput.s9.Value                54
 pfset GeomInput.s10.Value               55
 pfset GeomInput.s11.Value               58
 pfset GeomInput.s12.Value               60
 pfset GeomInput.s13.Value               62
 pfset GeomInput.s14.Value               72
 pfset GeomInput.s15.Value               81
 pfset GeomInput.s16.Value               86
 pfset GeomInput.g1.Value                1
 pfset GeomInput.g2.Value                2
 pfset GeomInput.g3.Value                4
 pfset GeomInput.g4.Value                6
 

#-----------------------------------------------------------------------------
# hydraulic conductivity (values in m/hr)
#-----------------------------------------------------------------------------
pfset Geom.Perm.Names                     "domain s1 s2 s4 s7 g1 g2 g3 g4"

pfset Geom.domain.Perm.Type           Constant
pfset Geom.domain.Perm.Value          0.126

pfset Geom.s1.Perm.Type               Constant
pfset Geom.s1.Perm.Value              0.36
 
pfset Geom.s2.Perm.Type               Constant
pfset Geom.s2.Perm.Value              1.8
 
pfset Geom.s4.Perm.Type               Constant
pfset Geom.s4.Perm.Value              0.0036
 
pfset Geom.s7.Perm.Type               Constant
pfset Geom.s7.Perm.Value              0.36
 
pfset Geom.g1.Perm.Type               Constant
pfset Geom.g1.Perm.Value              0.000036
 
pfset Geom.g2.Perm.Type               Constant
pfset Geom.g2.Perm.Value              0.0126
 
pfset Geom.g3.Perm.Type               Constant
pfset Geom.g3.Perm.Value              0.00036
 
pfset Geom.g4.Perm.Type               Constant
pfset Geom.g4.Perm.Value              0.0036

pfset Perm.TensorType                     TensorByGeom
pfset Geom.Perm.TensorByGeom.Names        "domain"
pfset Geom.domain.Perm.TensorValX         1.0d0
pfset Geom.domain.Perm.TensorValY         1.0d0
pfset Geom.domain.Perm.TensorValZ         1.0d0

#-----------------------------------------------------------------------------
# Specific Storage
#-----------------------------------------------------------------------------
pfset SpecificStorage.Type                Constant
pfset SpecificStorage.GeomNames           "domain"
pfset Geom.domain.SpecificStorage.Value   1.0e-5

#-----------------------------------------------------------------------------
# Phases
#-----------------------------------------------------------------------------
pfset Phase.Names                         "water"
pfset Phase.water.Density.Type            Constant
pfset Phase.water.Density.Value           1.0
pfset Phase.water.Viscosity.Type          Constant
pfset Phase.water.Viscosity.Value         1.0

#-----------------------------------------------------------------------------
# Contaminants
#-----------------------------------------------------------------------------
pfset Contaminants.Names                  ""

#-----------------------------------------------------------------------------
# Gravity
#-----------------------------------------------------------------------------
pfset Gravity                             1.0

#-----------------------------------------------------------------------------
# Timing (time units is set by units of permeability)
#-----------------------------------------------------------------------------
pfset TimingInfo.BaseUnit        1.0
pfset TimingInfo.StartCount      0.0
pfset TimingInfo.StartTime       0.0
pfset TimingInfo.StopTime        1000000.0
pfset TimingInfo.DumpInterval    1000.0

#pfset TimeStep.Type              Constant
#pfset TimeStep.Value             1.0

pfset TimeStep.Type              Growth
pfset TimeStep.InitialStep       0.0001
pfset TimeStep.GrowthFactor      1.4
pfset TimeStep.MaxStep           1000
pfset TimeStep.MinStep           0.0001

#-----------------------------------------------------------------------------
# Porosity
#-----------------------------------------------------------------------------
pfset Geom.Porosity.GeomNames             "domain s2 g1 g2 g3 g4"

pfset Geom.domain.Porosity.Type         Constant
pfset Geom.domain.Porosity.Value        0.3
 
pfset Geom.s2.Porosity.Type             Constant
pfset Geom.s2.Porosity.Value            0.4

pfset Geom.g1.Porosity.Type             Constant
pfset Geom.g1.Porosity.Value            0.01

pfset Geom.g2.Porosity.Type             Constant
pfset Geom.g2.Porosity.Value            0.2

pfset Geom.g3.Porosity.Type             Constant
pfset Geom.g3.Porosity.Value            0.05

pfset Geom.g4.Porosity.Type             Constant
pfset Geom.g4.Porosity.Value            0.2


#-----------------------------------------------------------------------------
# Domain
#-----------------------------------------------------------------------------
pfset Domain.GeomName                     "domain"

#----------------------------------------------------------------------------
# Mobility
#----------------------------------------------------------------------------
pfset Phase.water.Mobility.Type        Constant
pfset Phase.water.Mobility.Value       1.0

#-----------------------------------------------------------------------------
# Wells
#-----------------------------------------------------------------------------
pfset Wells.Names                         ""

#-----------------------------------------------------------------------------
# Time Cycles
#-----------------------------------------------------------------------------
pfset Cycle.Names "constant"
pfset Cycle.constant.Names              "alltime"
pfset Cycle.constant.alltime.Length      1
pfset Cycle.constant.Repeat             -1

#-----------------------------------------------------------------------------
# Boundary Conditions
#-----------------------------------------------------------------------------
pfset BCPressure.PatchNames                   [pfget Geom.domain.Patches]

pfset Patch.x-lower.BCPressure.Type		      FluxConst
pfset Patch.x-lower.BCPressure.Cycle		      "constant"
pfset Patch.x-lower.BCPressure.alltime.Value	      0.0

pfset Patch.y-lower.BCPressure.Type		      FluxConst
pfset Patch.y-lower.BCPressure.Cycle		      "constant"
pfset Patch.y-lower.BCPressure.alltime.Value	      0.0

pfset Patch.z-lower.BCPressure.Type		      FluxConst
pfset Patch.z-lower.BCPressure.Cycle		      "constant"
pfset Patch.z-lower.BCPressure.alltime.Value	      0.0

pfset Patch.x-upper.BCPressure.Type		      FluxConst
pfset Patch.x-upper.BCPressure.Cycle		      "constant"
pfset Patch.x-upper.BCPressure.alltime.Value	      0.0

pfset Patch.y-upper.BCPressure.Type		      FluxConst
pfset Patch.y-upper.BCPressure.Cycle		      "constant"
pfset Patch.y-upper.BCPressure.alltime.Value	      0.0

## overland flow boundary condition with spatially-distributed recharge from USGS report
pfset Patch.z-upper.BCPressure.Type		             OverlandFlow
pfset Patch.z-upper.BCPressure.Cycle		            "constant"
pfset Patch.z-upper.BCPressure.alltime.Value	      0.0

#-----------------------------------------------------------------------------
# Topo slopes in x-direction
#-----------------------------------------------------------------------------
pfset TopoSlopesX.Type                                "PFBFile"
pfset TopoSlopesX.GeomNames                           "domain"
pfset TopoSlopesX.FileName                            "tucson.slopex.pfb"

#-----------------------------------------------------------------------------
# Topo slopes in y-direction
#-----------------------------------------------------------------------------
pfset TopoSlopesY.Type                                "PFBFile"
pfset TopoSlopesY.GeomNames                           "domain"
pfset TopoSlopesY.FileName                            "tucson.slopey.pfb"

#-----------------------------------------------------------------------------
# Mannings coefficient
#-----------------------------------------------------------------------------
pfset Mannings.Type                                   "Constant"
pfset Mannings.GeomNames                              "domain"
pfset Mannings.Geom.domain.Value                       5.52e-6

#-----------------------------------------------------------------------------
# Relative Permeability
#-----------------------------------------------------------------------------
pfset Phase.RelPerm.Type                  VanGenuchten
pfset Phase.RelPerm.GeomNames             "domain s2 s3 s4 s7 s9 s14 g1 g2 g3 g4"

pfset Geom.domain.RelPerm.Alpha           2.691
pfset Geom.domain.RelPerm.N               2.0

pfset Geom.s2.RelPerm.Alpha         2.691
pfset Geom.s2.RelPerm.N             3

pfset Geom.s3.RelPerm.Alpha         3
pfset Geom.s3.RelPerm.N             2

pfset Geom.s4.RelPerm.Alpha         0.436
pfset Geom.s4.RelPerm.N             2
 
pfset Geom.s7.RelPerm.Alpha         3
pfset Geom.s7.RelPerm.N             2
 
pfset Geom.s9.RelPerm.Alpha         1.585
pfset Geom.s9.RelPerm.N             2

pfset Geom.s14.RelPerm.Alpha        1.585
pfset Geom.s14.RelPerm.N             2

pfset Geom.g1.RelPerm.Alpha         3
pfset Geom.g1.RelPerm.N             2

pfset Geom.g2.RelPerm.Alpha         3
pfset Geom.g2.RelPerm.N             2

pfset Geom.g3.RelPerm.Alpha         3
pfset Geom.g3.RelPerm.N             2

pfset Geom.g4.RelPerm.Alpha         3
pfset Geom.g4.RelPerm.N             2

#-----------------------------------------------------------------------------
# Saturation
#-----------------------------------------------------------------------------
pfset Phase.Saturation.Type               VanGenuchten
pfset Phase.Saturation.GeomNames          "domain s2 s3 s4 s7 s9 s14 g1 g2 g3 g4"

pfset Geom.domain.Saturation.Alpha        2.691
pfset Geom.domain.Saturation.N            2
pfset Geom.domain.Saturation.SRes         0.001
pfset Geom.domain.Saturation.SSat         1.0

pfset Geom.s2.Saturation.Alpha        2.691
pfset Geom.s2.Saturation.N            3
pfset Geom.s2.Saturation.SRes         0.001
pfset Geom.s2.Saturation.SSat         1.0
 
pfset Geom.s3.Saturation.Alpha        3
pfset Geom.s3.Saturation.N            2
pfset Geom.s3.Saturation.SRes         0.001
pfset Geom.s3.Saturation.SSat         1.0
 
pfset Geom.s4.Saturation.Alpha        0.436
pfset Geom.s4.Saturation.N            2
pfset Geom.s4.Saturation.SRes         0.001
pfset Geom.s4.Saturation.SSat         1.0
 
pfset Geom.s7.Saturation.Alpha        3
pfset Geom.s7.Saturation.N            2
pfset Geom.s7.Saturation.SRes         0.001
pfset Geom.s7.Saturation.SSat         1.0
 
pfset Geom.s9.Saturation.Alpha        1.585
pfset Geom.s9.Saturation.N            2
pfset Geom.s9.Saturation.SRes         0.001
pfset Geom.s9.Saturation.SSat         1.0

pfset Geom.s14.Saturation.Alpha       1.585
pfset Geom.s14.Saturation.N            2
pfset Geom.s14.Saturation.SRes         0.001
pfset Geom.s14.Saturation.SSat         1.0

pfset Geom.g1.Saturation.Alpha        3
pfset Geom.g1.Saturation.N            2
pfset Geom.g1.Saturation.SRes         0.001
pfset Geom.g1.Saturation.SSat         1.0

pfset Geom.g2.Saturation.Alpha        3
pfset Geom.g2.Saturation.N            2
pfset Geom.g2.Saturation.SRes         0.001
pfset Geom.g2.Saturation.SSat         1.0

pfset Geom.g3.Saturation.Alpha        3
pfset Geom.g3.Saturation.N            2
pfset Geom.g3.Saturation.SRes         0.001
pfset Geom.g3.Saturation.SSat         1.0

pfset Geom.g4.Saturation.Alpha        3
pfset Geom.g4.Saturation.N            2
pfset Geom.g4.Saturation.SRes         0.001
pfset Geom.g4.Saturation.SSat         1.0

#-----------------------------------------------------------------------------
# Phase sources:
#-----------------------------------------------------------------------------
pfset PhaseSources.water.Type                         "Constant"
pfset PhaseSources.water.GeomNames                    "domain"
pfset PhaseSources.water.Geom.domain.Value            0.0

#----------------------------------------------------------------
# CLM Settings:
# ------------------------------------------------------------
#for spin-up runs, CLM is initially turned off
pfset Solver.LSM                                        none
pfset Solver.CLM.CLMFileDir                           "clm_output/"
pfset Solver.CLM.Print1dOut                           False
pfset Solver.BinaryOutDir                             False
pfset Solver.CLM.CLMDumpInterval                      1
 
pfset Solver.CLM.MetForcing                           NC
pfset Solver.CLM.MetFileName                          "met_forcing_2000.nc"
pfset Solver.CLM.MetFilePath                          "../clm_input/"
pfset Solver.CLM.MetFileNT                            1
pfset Solver.CLM.IstepStart                           1
 
pfset Solver.CLM.EvapBeta                             Linear
pfset Solver.CLM.VegWaterStress                       Saturation
pfset Solver.CLM.ResSat                               0.1
pfset Solver.CLM.WiltingPoint                         0.12
pfset Solver.CLM.FieldCapacity                        0.98
pfset Solver.CLM.IrrigationType                       none

pfset Solver.EvapTransFile                            True
pfset Solver.EvapTrans.FileName                       "eff_recharge_0013.pfb"

#---------------------------------------------------------
# Initial conditions: water pressure
#---------------------------------------------------------
# pfset ICPressure.Type                                 NCFile
# pfset ICPressure.GeomNames                            domain
# pfset Geom.domain.ICPressure.RefPatch                 z-upper
# pfset Geom.domain.ICPressure.FileName                 press.init.nc

#pfset ICPressure.Type	 	                               PFBFile
#pfset ICPressure.GeomNames                             domain
#pfset Geom.domain.ICPressure.RefPatch                  z-upper
#pfset Geom.domain.ICPressure.FileName                  tucson.out.press.00030.pfb

pfset ICPressure.Type                                   HydroStaticPatch
pfset ICPressure.GeomNames                              domain
pfset Geom.domain.ICPressure.Value                      0.0
pfset Geom.domain.ICPressure.RefGeom                    domain
pfset Geom.domain.ICPressure.RefPatch                   z-lower

#----------------------------------------------------------------
# Outputs
# ------------------------------------------------------------
#Writing output (all pfb):
pfset Solver.PrintSubsurfData                         True
pfset Solver.PrintPressure                            True
pfset Solver.PrintSaturation                          True
pfset Solver.PrintMask                                True

pfset Solver.WriteCLMBinary                           False
pfset Solver.PrintCLM                                 False
pfset Solver.WriteSiloSpecificStorage                 False
pfset Solver.WriteSiloMannings                        False
pfset Solver.WriteSiloMask                            False
pfset Solver.WriteSiloSlopes                          False
pfset Solver.WriteSiloSubsurfData                     False
pfset Solver.WriteSiloPressure                        False
pfset Solver.WriteSiloSaturation                      False
pfset Solver.WriteSiloEvapTrans                       False
pfset Solver.WriteSiloEvapTransSum                    False
pfset Solver.WriteSiloOverlandSum                     False
pfset Solver.WriteSiloCLM                             False


#-----------------------------------------------------------------------------
# Exact solution specification for error calculations
#-----------------------------------------------------------------------------
pfset KnownSolution                                   NoKnownSolution

#-----------------------------------------------------------------------------
# Set solver parameters
#-----------------------------------------------------------------------------
# ParFlow Solution
pfset Solver                                          Richards
pfset Solver.TerrainFollowingGrid                     True
pfset Solver.Nonlinear.VariableDz                     True
pfset dzScale.GeomNames                               domain
pfset dzScale.Type                                    nzList
pfset dzScale.nzListNumber                            17
pfset Cell.0.dzScale.Value                             2.0
pfset Cell.1.dzScale.Value                             1.0
pfset Cell.2.dzScale.Value                             0.5
pfset Cell.3.dzScale.Value                             0.5
pfset Cell.4.dzScale.Value                             0.2
pfset Cell.5.dzScale.Value                             0.2
pfset Cell.6.dzScale.Value                             0.2
pfset Cell.7.dzScale.Value                             0.2
pfset Cell.8.dzScale.Value                             0.1
pfset Cell.9.dzScale.Value                             0.02
pfset Cell.10.dzScale.Value                            0.02
pfset Cell.11.dzScale.Value                            0.02
pfset Cell.12.dzScale.Value                            0.02
pfset Cell.13.dzScale.Value                            0.01
pfset Cell.14.dzScale.Value                            0.006
pfset Cell.15.dzScale.Value                            0.003
pfset Cell.16.dzScale.Value                            0.001

pfset Solver.MaxIter                                  25000
pfset Solver.Drop                                     1E-20
pfset Solver.AbsTol                                   1E-8
pfset Solver.MaxConvergenceFailures                   8
pfset Solver.Nonlinear.MaxIter                        80
pfset Solver.Nonlinear.ResidualTol                    1e-6

## new solver settings for Terrain Following Grid
pfset Solver.Nonlinear.EtaChoice                         EtaConstant
pfset Solver.Nonlinear.EtaValue                          0.001
pfset Solver.Nonlinear.UseJacobian                       True 
pfset Solver.Nonlinear.DerivativeEpsilon                 1e-16
pfset Solver.Nonlinear.StepTol	                         1e-30
pfset Solver.Nonlinear.Globalization                     LineSearch
pfset Solver.Linear.KrylovDimension                      70
pfset Solver.Linear.MaxRestarts                          2

pfset Solver.Linear.Preconditioner                       PFMG
pfset Solver.Linear.Preconditioner.PCMatrixType          FullJacobian

#keys for first round of spin-up
pfset OverlandFlowSpinUp                             1
pfset OverlandSpinupDampP1                           10.0
pfset OverlandSpinupDampP2                           0.1

 pfset NetCDF.NumStepsPerFile	                         24
 pfset NetCDF.CLMNumStepsPerFile                       24
 pfset NetCDF.WritePressure	                           False
 pfset NetCDF.WriteSaturation	                         False
 pfset NetCDF.WriteMannings                            False
 pfset NetCDF.WriteSubsurface	                         False
 pfset NetCDF.WriteSlopes	                             False
 pfset NetCDF.WriteMask                                False
 pfset NetCDF.WriteDZMultiplier                        False
 pfset NetCDF.WriteEvapTrans	                          False
 pfset NetCDF.WriteEvapTransSum                        False
 pfset NetCDF.WriteOverlandSum	                        False
 pfset NetCDF.WriteOverlandBCFlux	                     False
 pfset NetCDF.WriteCLM		                               False

#-----------------------------------------------------------------------------
# Distribute inputs
#-----------------------------------------------------------------------------
pfset ComputationalGrid.NX                246 
pfset ComputationalGrid.NY                178 
pfset ComputationalGrid.NZ                1
pfdist tucson.slopex.pfb
pfdist tucson.slopey.pfb

pfset ComputationalGrid.NX                246 
pfset ComputationalGrid.NY                178 
pfset ComputationalGrid.NZ                17
pfdist geology_indicator.pfb
pfdist eff_recharge_0013.pfb

#-----------------------------------------------------------------------------
# Run Simulation 
#-----------------------------------------------------------------------------
set runname "sabino"
puts $runname
#pfwritedb $runname
pfrun    $runname

##-----------------------------------------------------------------------------
## Undistribute outputs
##-----------------------------------------------------------------------------
pfundist $runname
#pfundist press.init.pfb
pfundist tucson.slopex.pfb
pfundist tucson.slopey.pfb
pfundist geology_indicator.pfb
pfundist eff_recharge_0013.pfb

puts "ParFlow run Complete"

#
