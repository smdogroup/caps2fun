"""
Sean Engelstad, Brian Burke, Graeme J Kennedy
Georgia Tech SMDO 2022

Goal is to generate ESP/CAPS geometry files
"""

from re import T
from caps2fun.caps2fun.caps2fun import Caps2Fun
import pyCAPS
import mpi4py
from typing import TYPE_CHECKING, List, TypedDict

class CAPS:
    """
    wrapper function class for Caps functions
    TODO : add wrapper function @CAPS for that only lets you run it for root procs or no procs
    """
    def __init__(self, func, obj):
        self._func = func
        self._obj = obj # something like this or func.__self__ is the obj

    def __call__(self):
        if (self._obj.can_use_caps):
            # call function
            pass

# TODO : alternative wrapper classes CAPS_Fluid, CAPS_Struct that check not None

class ESPGeometry:
    def __init__(self, csm_file : str, comm):
        assert(".csm" in csm_file)

        self.csm_file = csm_file
        self.comm = comm

        self.no_comm : bool = self.comm is None
        self.root_proc : bool = self.comm.Get_rank() == 0
        
        self.caps_struct = None
        self.caps_fluid = None
        self.tacs_aim = None
        self.egads_struct_aim = None
        self.egads_fluid_aim = None
        self.pointwise_aim = None
        self.fun3d_aim = None

        # booleans for preparing struct analysis from tacs AIM
        self.set_struct_mesh_settings = False
        self.has_set_struct_analysis = False
        self.has_struct_materials = False
        self.has_struct_properties = False
        self.has_struct_constraints = False
        self.has_struct_DVs = False

        # full dictionaries to be sent to tacsAIMs
        self._constraints = {}
        self._materials = {}
        self._properties = {}
        self._DVs = {}
        self._DVRs = {}

    
    @property
    def root_proc(self) -> bool:
        return self.comm.Get_rank() == 0

    @property
    def can_use_caps(self) -> bool:
        return self.comm is None or self.root_proc

    @CAPS
    def start_caps_struct(self, start_aims : bool = True) -> None:
        """
        start a CAPS_Struct problem from the struct mode of the .csm file
        """
        self.caps_struct = pyCAPS.Problem(problemName = "CAPS_struct",
                    capsFile = self.csm_file,
                    outLevel = 1)
        assert("structOn" in self.capsFluid.geometry.cfgpmtr.keys())
        # maybe add a similar attribute StructOn

        if (start_aims): self.start_struct_aims()

    @CAPS
    def start_caps_fluid(self, start_aims:bool = True) -> None:
        """
        start a CAPS_Fluid problem from the fluid mode of the .csm file
        """
        self.caps_fluid = pyCAPS.Problem(problemName = "CAPS_fluid",
                        capsFile = self.config["csmFile"],
                        outLevel = 1)
        
        assert("cfdOn" in self.capsFluid.geometry.cfgpmtr.keys())
        self.caps_fluid.geometry.cfgpmtr["cfdOn"].value = 1
        # maybe add a FluidOn setting in each

        if (start_aims): self.start_fluid_aims()

    def start_struct_aims(self) -> None:
        """
        initialize each of the CAPS_struct aims
        """
        self.start_egads_struct_aim()
        self.start_tacs_aim()

    def start_fluid_aims(self) -> None:
        """
        initialize each of the CAPS_fluid aims
        """
        self.start_egads_fluid_aim()
        self.start_pointwise_aim()
        self.start_fun3d_aim()

    @CAPS
    def start_tacs_aim(self) -> None:
        """
        start the tacs AIM
        requires capsAIM attribute tacsAIM in list
        """
        assert(self.caps_struct is not None)
        self.caps_struct.analysis.create(aim = "tacsAIM", name = "tacs")

    @CAPS
    def start_egads_struct_aim(self) -> None:
        """
        start the EGADS struct AIM for structure meshing
        requires capsAIM attribute egadsAIM in list and applied to structure
        """
        assert(self.caps_struct is not None)
        self.egads_struct_aim = self.caps_struct.analysis.create(aim="egadsTessAIM")

    @CAPS
    def start_egads_fluid_aim(self) -> None:
        """
        start the egads fluid AIM
        requires capsAIM attribute egadsAIM in list and applied to fluid
        """
        assert(self.caps_fluid is not None)
        self.egads_struct_aim = self.caps_struct.analysis.create(aim="egadsTessAIM")

    @CAPS
    def start_pointwise_aim(self) -> None:
        """
        start the pointwise meshing aim
        requires capsAIM attribute pointwiseAIM in list and applied to fluid
        """
        assert(self.caps_fluid is not None)
        self.pointwise_aim = self.caps_fluid.analysis.create(aim = "pointwiseAIM",
                                            name = "pointwise")

    @CAPS
    def start_fun3d_aim(self) -> None:
        """
        start the fun3d AIM
        requires capsAIM attribute fun3dAIM in list and applied to fluid
        """
        assert(self.caps_fluid is not None)
        self.fun3d_aim = self.caps_fluid.analysis.create(aim = "fun3dAIM",
                                    name = "fun3d")

    def struct_mesh_settings(
        self, mesh_style : str, 
        edge_pt_min : float = 5,
        edge_pt_max : float = 100,
        tess_params : List[float] = [0.03, 0.01, 15],
    ) -> None:

        self.set_struct_mesh_settings = True

        #names the bdf and dat files as pointwise.ext or tetgen.ext
        self.tacs_aim.input.Proj_Name = mesh_style
        
        #Egads Aim section, for mesh
        self.egads_struct_aim.input.Edge_Point_Min = edge_pt_min
        self.egads_struct_aim.input.Edge_Point_Max = edge_pt_max

        self.egads_struct_aim.input.Mesh_Elements = "Quad"

        assert(len(tess_params) == 3)
        self.egads_struct_aim.input.Tess_Params = tess_params

        #increase the precision in the BDF file
        self.tacs_aim.input.File_Format = "Large"
        self.tacs_aim.input.Mesh_File_Format = "Large"

        # Link the mesh
        self.tacs_aim.input["Mesh"].link(self.egadsAim.output["Surface_Mesh"])

    def add_material(self, name : str, young_modulus : float, poisson_ratio : float, 
                        density : float, tension_allowed : float, type : str = "isotropic") -> None:
        self.has_struct_materials = True

        material_dict = {"materialType" : type,
                            "youngModulus" : young_modulus ,
                            "poissonRatio": poisson_ratio,
                            "density" : density,
                            "tensionAllow" :  tension_allowed}
        self._materials[name] = material_dict

    @property
    def ready_for_struct_analysis(self):
        return self.has_struct_constraints and self.has_struct_materials and self.has_struct_properties and \
            self.set_struct_mesh_settings and self.has_set_struct_analysis and self.has_struct_DVs

    def _add_struct_config_settings(self):
        """
        add the struct config settings to the tacs AIM
        """

        if (self.ready_for_struct_analysis):
            self.tacs_aim.input.Material = self._materials
            self.tacs_aim.input.Constraint = self._constraints
            self.tacs_aim.input.Property = self._properties
            self.tacs_aim.input.Design_Variable = self._DVs
            self.tacs_aim.input.Design_Variable_Relation = self._DVRs
        else:
            raise AssertionError("Haven't set up the struct problem correctly: need materials, constraints, properties, analysis, and mesh settings.")

    def get_shell_property(self, thickness : float, material : str, boost_factor : float = 1.0) -> dict:
        """
        get a shell property dict
        """
        return {"propertyType" : "Shell",
                "membraneThickness" : thickness,
                "material"        : "aluminum",
                "bendingInertiaRatio" : 1.0 * boost_factor, # Default, TODO : double-check which one boost_factor should be on or both
                "shearMembraneRatio"  : 5.0/6.0} # Default

    def add_struct_constraint(self, name : str, caps_constraint : str, dof_values : List[int]) -> None:
        self.has_struct_constraints = True
        self._constraints[name] = {"groupName" : caps_constraint,
                "dofConstraint" : dof_values}

    def set_struct_analysis(self, analysis_type = "Static"):
        self.has_set_struct_analysis = True
        self.tacsAim.input.Analysis_Type = analysis_type

    def add_thickness_DV(self, name : str, caps_group : str, thickness : float, add_DVR=True):
        #thick DV dictionary for Design_Variable Dict
        self.has_struct_DVs = True
        DVdict = {"groupName" : caps_group,
                "initialValue" : thickness,
                "lowerBound" : thickness*0.5,
                "upperBound" : thickness*1.5,
                "maxDelta"   : thickness*0.1}

        DVRdict = {"variableType": "Property",
                    "fieldName" : "T",
                    "constantCoeff" : 0.0,
                    "groupName" : name,
                    "linearCoeff" : 1.0}

        self._DVs[name] = DVdict
        self._DVRs[name] = DVRdict

    def fluid_mesh_settings(self):
        #wall bc settings (wall is the OML)
            if (self.config["fun3d_analysis_type"] == "inviscid"):
                self.wallBC = {"bcType" : "inviscid"}
                wallSpacing = 0.1
            elif (self.config["fun3d_analysis_type"] == "laminar"):
                wallSpacing = 0.01
                self.wallBC = {"bcType" : "viscous",
                "boundaryLayerSpacing" : wallSpacing}
            elif (self.config["fun3d_analysis_type"] == "turbulent"):
                wallSpacing = 0.01
                self.wallBC = {"bcType" : "viscous",
                "boundaryLayerSpacing" : wallSpacing}

            # Dump VTK files for visualization
            self.pointwiseAim.input.Proj_Name   = "TransportWing"
            self.pointwiseAim.input.Mesh_Format = "VTK"

            # Connector level
            self.pointwiseAim.input.Connector_Turn_Angle       = 1
            self.pointwiseAim.input.Connector_Prox_Growth_Rate = 1.2
            self.pointwiseAim.input.Connector_Source_Spacing   = True

            # Domain level
            self.pointwiseAim.input.Domain_Algorithm    = "AdvancingFront"
            self.pointwiseAim.input.Domain_Max_Layers   = self.config["domain_max_layers"]
            self.pointwiseAim.input.Domain_Growth_Rate  = self.config["domain_growth_rate"]
            self.pointwiseAim.input.Domain_TRex_ARLimit = 40.0 #def 40.0, lower inc mesh size
            self.pointwiseAim.input.Domain_Decay        = 0.5
            self.pointwiseAim.input.Domain_Iso_Type = "Triangle" #"TriangleQuad"
            self.pointwiseAim.input.Domain_Wall_Spacing = wallSpacing

            # Block level
            self.pointwiseAim.input.Block_Boundary_Decay       = self.config["block_boundary_decay"]
            self.pointwiseAim.input.Block_Collision_Buffer     = 1.0
            self.pointwiseAim.input.Block_Max_Skew_Angle       = self.config["block_max_skew_angle"]
            self.pointwiseAim.input.Block_Edge_Max_Growth_Rate = 2.0
            self.pointwiseAim.input.Block_Full_Layers          = self.config["block_full_layers"]
            self.pointwiseAim.input.Block_Max_Layers           = self.config["block_max_layers"]
            self.pointwiseAim.input.Block_TRexType = "TetPyramid"
            #T-Rex cell type (TetPyramid, TetPyramidPrismHex, AllAndConvertWallDoms)        

        
            self.pointwiseAim.input.Mesh_Sizing = {"wall": self.wallBC,
                "Farfield": {"bcType":"Farfield"}}
