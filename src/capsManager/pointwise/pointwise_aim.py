
__all__ = ["PointwiseAim"]

from typing import TYPE_CHECKING
import pyCAPS
import os

class PointwiseAim:
    def __init__(self, caps_problem:pyCAPS.Problem):
        """
        Pointwise AIM wrapper object which builds automatic fluid meshes through pointwise
            use run_pointwise() method to 
        """

        # initialize the pointwise aim
        self._aim = caps_problem.analysis.create(aim = "pointwiseAIM", name = "pointwise")

        self._is_setup = False
        self._built_mesh = False

        self._wall_bc = None
        self._boundary_condition = None

    @property
    def is_setup(self) -> bool:
        return self._is_setup

    def set_mesh(self,
        wall_spacing:float=0.01, 
        domain_max_layers:int=15,
        domain_growth_rate:float=1.75,
        block_boundary_decay:float=0.5,
        block_max_skew_angle:float=170.0,
        block_full_layers:int=1,
        block_max_layers:int=100,
        inviscid:bool=True
    ):
        """
        TODO : add parameters maybe here which set the mesh or auto-generate?
        """

        # Dump VTK files for visualization
        self.aim.input.Proj_Name   = "capsFluid"
        self.aim.input.Mesh_Format = "VTK"

        # Connector level
        self.aim.input.Connector_Turn_Angle       = 1
        self.aim.input.Connector_Prox_Growth_Rate = 1.2
        self.aim.input.Connector_Source_Spacing   = True

        # Domain level
        self.aim.input.Domain_Algorithm    = "AdvancingFront"
        self.aim.input.Domain_Max_Layers   = domain_max_layers
        self.aim.input.Domain_Growth_Rate  = domain_growth_rate
        self.aim.input.Domain_TRex_ARLimit = 40.0 #def 40.0, lower inc mesh size
        self.aim.input.Domain_Decay        = 0.5
        self.aim.input.Domain_Iso_Type = "Triangle" #"TriangleQuad"
        self.aim.input.Domain_Wall_Spacing = wall_spacing
        # Block level
        self.aim.input.Block_Boundary_Decay       = block_boundary_decay
        self.aim.input.Block_Collision_Buffer     = 1.0
        self.aim.input.Block_Max_Skew_Angle       = block_max_skew_angle
        self.aim.input.Block_Edge_Max_Growth_Rate = 2.0
        self.aim.input.Block_Full_Layers          = block_full_layers
        self.aim.input.Block_Max_Layers           = block_max_layers
        self.aim.input.Block_TRexType = "TetPyramid"
        #T-Rex cell type (TetPyramid, TetPyramidPrismHex, AllAndConvertWallDoms)        

        # maybe change the constraint names or something or have a class for this?
        self.aim.input.Mesh_Sizing = self.boundary_condition

        if inviscid:
            self._wall_bc = {"bcType" : "inviscid"}
        else:
            self._wall_bc = {"bcType" : "viscous",
                "boundaryLayerSpacing" : wall_spacing}

        self._is_setup = True
        self._built_mesh = False

    @property
    def aim(self):
        return self._aim

    @property
    def built_mesh(self) -> bool:
        return self._built_mesh

    @property
    def analysis_dir(self) -> str:
        return self.aim.analysisDir

    def run_pointwise(self) -> bool:
        """
        run automatic pointwise meshing through ESP/CAPS pointwise AIM
        """

        #run AIM pre-analysis
        self.aim.preAnalysis()

        try:
            CAPS_GLYPH = os.environ["CAPS_GLYPH"]
            #for i in range(1): #can run extra times if having license issues
            os.system(f"pointwise -b {CAPS_GLYPH}/GeomToMesh.glf {self.analysis_dir}/caps.egads {self.analysis_dir}/capsUserDefaults.glf")
                #ranPointwise = os.path.isfile('caps.GeomToMesh.gma') and os.path.isfile('caps.GeomToMesh.ugrid')
                #if ranPointwise: break
            ran_pointwise = True
        except:
            ran_pointwise = False
            raise AssertionError("Pointwise mesh generation failed to run...")

        if ran_pointwise:
            #run AIM postanalysis, files in self.aim.analysisDir
            self.aim.postAnalysis() 

        self._built_mesh = ran_pointwise

        return ran_pointwise


    @property
    def boundary_condition(self) -> dict:
        if self._boundary_condition is None:
            return {"wall": self._wall_bc,
            "Farfield": {"bcType":"Farfield"}}
        else:
            return self._boundary_condition

    @boundary_condition.setter
    def boundary_condition(self, new_bc:dict):
        self._boundary_condition = new_bc