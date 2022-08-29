
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
        block_max_layers:int=100,
        inviscid:bool=True
        ):
        """
        Pointwise AIM wrapper object which builds automatic fluid meshes through pointwise
            use run_pointwise() method to 
        """

        # initialize the pointwise aim
        self._aim = caps_problem.analysis.create(aim = "pointwiseAIM", name = "pointwise")

        self._setup_settings = False
        self._setup_mesh = False

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

        # auto set mesh
        self.set_mesh()

    def set_mesh(self):
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
        self.aim.input.Domain_Max_Layers   = self._domain_max_layers
        self.aim.input.Domain_Growth_Rate  = self._domain_growth_rate
        self.aim.input.Domain_TRex_ARLimit = 40.0 #def 40.0, lower inc mesh size
        self.aim.input.Domain_Decay        = 0.5
        self.aim.input.Domain_Iso_Type = "Triangle" #"TriangleQuad"
        self.aim.input.Domain_Wall_Spacing = self._wall_spacing
        # Block level
        self.aim.input.Block_Boundary_Decay       = self._block_boundary_decay
        self.aim.input.Block_Collision_Buffer     = 1.0
        self.aim.input.Block_Max_Skew_Angle       = self._block_max_skew_angle
        self.aim.input.Block_Edge_Max_Growth_Rate = 2.0
        self.aim.input.Block_Full_Layers          = self._block_full_layers
        self.aim.input.Block_Max_Layers           = self._block_max_layers
        self.aim.input.Block_TRexType = "TetPyramid"
        #T-Rex cell type (TetPyramid, TetPyramidPrismHex, AllAndConvertWallDoms)        

        # maybe change the constraint names or something or have a class for this?
        self.aim.input.Mesh_Sizing = self.mesh_sizing

        self._setup_settings = True

    @property
    def aim(self):
        return self._aim

    @property
    def is_setup(self) -> bool:
        return self._setup_settings and self._setup_mesh

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

        #run AIM postanalysis, files in self.aim.analysisDir
        self.aim.postAnalysis() 

        self._setup_mesh = ran_pointwise

        return ran_pointwise


    @property
    def mesh_sizing(self) -> dict:
        return {"wall": self._wall_bc,
            "Farfield": {"bcType":"Farfield"}}