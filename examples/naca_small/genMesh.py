#------------------------------------------------------------------------------#
import pyCAPS
import os

#------------------------------------------------------------------------------#

meshStructure = True

def run_pointwise(pointwise):
    # Run AIM pre-analysis
    pointwise.preAnalysis()

    ####### Run pointwise #################
    currentDirectory = os.getcwd() # Get current working directory
    os.chdir(pointwise.analysisDir)    # Move into test directory

    CAPS_GLYPH = os.environ["CAPS_GLYPH"]
    os.system("pointwise -b " + CAPS_GLYPH + "/GeomToMesh.glf caps.egads capsUserDefaults.glf")

    os.chdir(currentDirectory)     # Move back to top directory
    #######################################

    # Run AIM post-analysis
    pointwise.postAnalysis()

#------------------------------------------------------------------------------#
# Load CSM file
#os.system("serveCSM naca_small.csm")


if (meshFluid):
    filename = os.path.join("cfd.egads")
    caps = pyCAPS.Problem(problemName = "myCAPS",
                        capsFile = filename,
                        outLevel = 1)
    # Create pointwise aim
    pointwise = caps.analysis.create(aim = "pointwiseAIM",
                                        name = "pointwise")

    # Dump VTK files for visualization
    pointwise.input.Proj_Name   = "TransportWing"
    pointwise.input.Mesh_Format = "VTK"

    # Connector level
    pointwise.input.Connector_Turn_Angle       = 10
    pointwise.input.Connector_Prox_Growth_Rate = 1.2
    pointwise.input.Connector_Source_Spacing   = True

    # Domain level
    pointwise.input.Domain_Algorithm    = "AdvancingFront"
    pointwise.input.Domain_Max_Layers   = 15
    pointwise.input.Domain_Growth_Rate  = 1.25
    pointwise.input.Domain_TRex_ARLimit = 40.0
    pointwise.input.Domain_Decay        = 0.8
    pointwise.input.Domain_Iso_Type = "Triangle"

    # Block level
    pointwise.input.Block_Boundary_Decay       = 0.8
    pointwise.input.Block_Collision_Buffer     = 1.0
    pointwise.input.Block_Max_Skew_Angle       = 160.0
    pointwise.input.Block_Edge_Max_Growth_Rate = 1.5
    pointwise.input.Block_Full_Layers          = 1
    pointwise.input.Block_Max_Layers           = 100
    pointwise.input.Block_TRexType = "TetPyramid"
    #T-Rex cell type (TetPyramid, TetPyramidPrismHex, AllAndConvertWallDoms).

    # Set wall spacing for capsMesh == leftWing and capsMesh == riteWing
    viscousWall  = {"boundaryLayerSpacing" : 0.001}
    pointwise.input.Mesh_Sizing = {"OML": viscousWall,
            "Farfield": {"bcType":"farfield"}}

    # Execute pointwise
    run_pointwise(pointwise)

    # View the surface tessellation
    #pointwise.geometry.view()
    #os.system("pointwise ")

##-------------------------------------------------------##
if (meshStructure):
    filename = os.path.join("wing.egads")
    caps = pyCAPS.Problem(problemName = "myCAPS",
                        capsFile = filename,
                        outLevel = 1)

    #initialize egads Aim
    egadsAim = caps.analysis.create(aim="egadsTessAIM")

    #initialize tacs Aim
    tacsAim = caps.analysis.create(aim = "tacsAIM", name = "tacs")


    ##-----------PyCAPS setup function-------##
    #setup function for naca_small.csm
        
    #Egads Aim section, for mesh
    egadsAim.input.Edge_Point_Min = 5
    egadsAim.input.Edge_Point_Max = 10

    egadsAim.input.Mesh_Elements = "Quad"

    egadsAim.input.Tess_Params = [.25,.01,15]

    #increase the precision in the BDF file
    tacsAim.input.File_Format = "Large"
    tacsAim.input.Mesh_File_Format = "Large"

    # Link the mesh
    tacsAim.input["Mesh"].link(egadsAim.output["Surface_Mesh"])

    # Set analysis type
    tacsAim.input.Analysis_Type = "Static"

    #materials section    
    madeupium    = {"materialType" : "isotropic",
                    "youngModulus" : 72.0E9 ,
                    "poissonRatio": 0.33,
                    "density" : 2.8E3,
                    "tensionAllow" :  20.0e7}

    tacsAim.input.Material = {"madeupium": madeupium}

    # Material properties section
    OMLshell = {"propertyType" : "Shell",
                "membraneThickness" : 0.01,
                "material"        : "madeupium",
                "bendingInertiaRatio" : 1.0, # Default
                "shearMembraneRatio"  : 5.0/6.0} # Default
    ribshell  = {"propertyType" : "Shell",
                "membraneThickness" : 0.02,
                "material"        : "madeupium",
                "bendingInertiaRatio" : 1.0, # Default
                "shearMembraneRatio"  : 5.0/6.0} # Default

    sparshell = {"propertyType" : "Shell",
                "membraneThickness" : 0.05,
                "material"        : "madeupium",
                "bendingInertiaRatio" : 1.0, # Default
                "shearMembraneRatio"  : 5.0/6.0} # Default

    tacsAim.input.Property = {"rib": ribshell,
    "spar" : sparshell,
    "OML" : OMLshell}

    # constraint section
    constraint1 = {"groupName" : "wingRoot",
                    "dofConstraint" : 123456}

    tacsAim.input.Constraint = {"fixRoot": constraint1}

    # Set load
    liftload = {"groupName"         : "bottomWing",
        "loadType"          : "GridForce",
        "forceScaleFactor"  : 1.0e2,
        "directionVector"   : [0.0, 1.0, 0.0]}

    # Set loads
    tacsAim.input.Load = {"lift": liftload }

    #return the capsGroups you want to have thickDVs in that order
    #[thick1, thick2, thick3]
    capsDVgroups = ["rib", "spar", "OML"]

    ##-----------setup DVs and DVRs-----------##
    #setup the Design_Variable and Design_Variable_Relations

    #list of design variables, with thick1-thickN for thickness DVs
    desvars = ["area","aspect","taper","ctwist","lesweep","dihedral","thick1", "thick2", "thick3"]
    nvar = len(desvars)

    #where the thickDVs have thick1 is for "rib", thick2 for "spar" etc.
    capsGroups = ["rib","spar","OML"]
    thickness = [0.01, 0.02, 0.03]

    def makeThicknessDV(capsGroup, thickness):
        #thick DV dictionary for Design_Variable Dict
        desvar    = {"groupName" : capsGroup,
                "initialValue" : thickness,
                "lowerBound" : thickness*0.5,
                "upperBound" : thickness*1.5,
                "maxDelta"   : thickness*0.1}
        return desvar
        
    def makeThicknessDVR(DVname):
        #thick DV dictionary for Design_Variable_Relation Dict
        DVR = {"variableType": "Property",
        "fieldName" : "T",
        "constantCoeff" : 0.0,
        "groupName" : DVname,
        "linearCoeff" : 1.0}
        return DVR

    #make initial DV and DVRdict
    DVdict = {}
    DVRdict = {}

    #add thickDVs and geomDVs to caps
    thickCt = 0
    for desvar in desvars:
        if ("thick" in desvar):
            #add thickDV entry into DV_Relations and DV Dicts
            DVRdict[desvar] = makeThicknessDVR(desvar)
            DVdict[desvar] = makeThicknessDV(capsGroups[thickCt],thickness[thickCt])
            thickCt += 1
        else: #geomDV, add empty entry into DV dicts
            DVdict[desvar] = {}

    #input DVdict and DVRdict into tacsAim
    tacsAim.input.Design_Variable = DVdict
    tacsAim.input.Design_Variable_Relation = DVRdict

    ##---------PYTACS to run TACS-------------##
    #first run the preanalysis to prepare to run tacs through pytacs
    tacsAim.preAnalysis()

    #data file
    datFile = os.path.join(tacsAim.analysisDir, tacsAim.input.Proj_Name + '.dat')

    #