
__all__ = ["TacsMaterial"]

import numpy as np
from tacs.constitutive import MaterialProperties

class TacsMaterial:
    def __init__(self, material_card):
        self._card = material_card

    @property
    def card(self):
        return self._card
    
    @property
    def tacs_material(self):
        """
        getter method for tacs material objects
        """
        assert(self.material_type in ["MAT1", "MAT2", "MAT8"])
        if self.material_type == 'MAT1':
            return self.isotropic_material
        elif self.material_type == 'MAT8':
            return self.orthotropic_material
        elif self.material_type == 'MAT2':
            return self.anisotropic_material
        else:
            AssertionError("Material must be in [MAT1, MAT2, MAT8]")

    @property
    def material_type(self) -> str:
        return self.card.type

    @property
    def isotropic_material(self):
        """
        generate tacs isotropic material from nastran isotropic material card
        """
        return MaterialProperties(
                rho=self.card.rho, 
                E=self.card.e,
                nu=self.card.nu,
                ys=self.card.St,
                alpha=self.card.a
                )

    @property
    def orthotropic_material(self):
        """
        generate tacs orthotropic material from nastran orthotropic material card
        """
        E1 = self.card.e11
        E2 = self.card.e22
        nu12 = self.card.nu12
        G12 = self.card.g12
        G13 = self.card.g1z
        G23 = self.card.g2z
        # If out-of-plane shear values are 0, Nastran defaults them to the in-plane
        if G13 == 0.0:
            G13 = G12
        if G23 == 0.0:
            G23 = G12
        rho = self.card.rho
        Xt = self.card.Xt
        Xc = self.card.Xc
        Yt = self.card.Yt
        Yc = self.card.Yc
        S12 = self.card.S
        # TODO: add alpha
        return MaterialProperties(
            rho=rho, 
            E1=E1, 
            E2=E2, 
            nu12=nu12, 
            G12=G12, 
            G13=G13, 
            G23=G23,
            Xt=Xt, 
            Xc=Xc, 
            Yt=Yt, 
            Yc=Yc, 
            S12=S12)

    @property
    def anisotropic_material(self):
        C11 = self.card.G11
        C12 = self.card.G12
        C22 = self.card.G22
        C13 = self.card.G13
        C23 = self.card.G23
        C33 = self.card.G33
        rho = self.card.rho
        # See if this card features anisotropic coupling terms (which we don't support yet)
        if np.abs(C13)/(C11+C22) >= 1e-8 or np.abs(C23)/(C11+C22) >= 1e-8:
            self._TACSWarning(f"MAT2 card {self.card.mid} has anisotropic stiffness components that are not currently supported. "
                                "These terms will be dropped and the material treated as orthotropic. "
                                "Result accuracy may be affected.")
        nu12 = C12 / C22
        nu21 = C12 / C11
        E1 = C11 * (1 - nu12 * nu21)
        E2 = C22 * (1 - nu12 * nu21)
        G12 = G13 = G23 = C33
        return MaterialProperties(
            rho=rho, 
            E1=E1, 
            E2=E2, 
            nu12=nu12, 
            G12=G12, 
            G13=G13,
            G23=G23
            )