
__all__ = ["PointwiseAim"]

from typing import TYPE_CHECKING
import pyCAPS
import os

class PointwiseAim:
    def __init__(self, 
        caps_problem:pyCAPS.Problem, 
        wall_spacing:float=0.01, 
        domain_max_layers:int=15,
        domain_growth_rate:float=1.75,
        block_boundary_decay:float=0.5,
        block_max_skew_angle:float=170.0,
        block_full_layers:int=1,
        block_max_layers:int=100.0,
        inviscid:bool=True
        ):

        # initialize the pointwise aim
        self._aim = caps_problem.analysis.create(aim = "pointwiseAIM",
                                            name = "pointwise")

        self._wall_spacing = wall_spacing
        self._domain_max_layers = domain_max_layers
        self._domain_growth_rate = domain_growth_rate
        self._block_boundary_decay = block_boundary_decay
        self._block_max_skew_angle = block_max_skew_angle
        self._block_full_layers = block_full_layers
        self._block_max_layers = block_max_layers
        self._inviscid = inviscid

        if self._inviscid:
            self._wall_bc = {"bcType" : "inviscid"}
        else:
            self._wall_bc = {"bcType" : "viscous",
                "boundaryLayerSpacing" : self._wall_spacing}

    def set_mesh(self):
        # Dump VTK files for visualization
        self.pointwiseAim.input.Proj_Name   = "capsFluid"
        self.pointwiseAim.input.Mesh_Format = "VTK"

        # Connector level
        self.pointwiseAim.input.Connector_Turn_Angle       = 1
        self.pointwiseAim.input.Connector_Prox_Growth_Rate = 1.2
        self.pointwiseAim.input.Connector_Source_Spacing   = True

        # Domain level
        self.pointwiseAim.input.Domain_Algorithm    = "AdvancingFront"
        self.pointwiseAim.input.Domain_Max_Layers   = self._domain_max_layers
        self.pointwiseAim.input.Domain_Growth_Rate  = self._domain_growth_rate
        self.pointwiseAim.input.Domain_TRex_ARLimit = 40.0 #def 40.0, lower inc mesh size
        self.pointwiseAim.input.Domain_Decay        = 0.5
        self.pointwiseAim.input.Domain_Iso_Type = "Triangle" #"TriangleQuad"
        self.pointwiseAim.input.Domain_Wall_Spacing = self._wall_spacing
        # Block level
        self.pointwiseAim.input.Block_Boundary_Decay       = self._block_boundary_decay
        self.pointwiseAim.input.Block_Collision_Buffer     = 1.0
        self.pointwiseAim.input.Block_Max_Skew_Angle       = self._block_max_skew_angle
        self.pointwiseAim.input.Block_Edge_Max_Growth_Rate = 2.0
        self.pointwiseAim.input.Block_Full_Layers          = self._block_full_layers
        self.pointwiseAim.input.Block_Max_Layers           = self._block_max_layers
        self.pointwiseAim.input.Block_TRexType = "TetPyramid"
        #T-Rex cell type (TetPyramid, TetPyramidPrismHex, AllAndConvertWallDoms)        

        # maybe change the constraint names or something or have a class for this?
        self.pointwiseAim.input.Mesh_Sizing = self.mesh_sizing


    @property
    def aim(self):
        return self._aim

    @property
    def analysis_dir(self) -> str:
        return self.aim.analysisDir

    def run_pointwise(self) -> bool:
        
        ran_pointwise = False

        #run AIM pre-analysis
        self.aim.preAnalysis()

        try:
            CAPS_GLYPH = os.environ["CAPS_GLYPH"]
            #for i in range(1): #can run extra times if having license issues
            os.system(f"pointwise -b {CAPS_GLYPH} {self.aim.analysis_dir}/GeomToMesh.glf caps.egads capsUserDefaults.glf")
                #ranPointwise = os.path.isfile('caps.GeomToMesh.gma') and os.path.isfile('caps.GeomToMesh.ugrid')
                #if ranPointwise: break
            ran_pointwise = True
        except:
            pass

        #run AIM postanalysis, files in self.pointwiseAim.analysisDir
        self.aim.postAnalysis()   

        return ran_pointwise


    @property
    def mesh_sizing(self) -> dict:
        return {"wall": self._wall_bc,
            "Farfield": {"bcType":"Farfield"}}