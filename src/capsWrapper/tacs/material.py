
__all__ = ["Material", "Isotropic", "Orthotropic"]

from typing import TYPE_CHECKING

class Material:
    def __init__(self, 
        name:str, material_type:str, young_modulus:float, poisson_ratio:float, density:float, tension_allow:float,
        compression_allow:float=None, shear_allow:float=None, yield_allow:float=None, thermExpCoeff:float=None
    ):
        assert(material_type in ["Isotropic", "Anisothotropic", "Orthotropic", "Anisotropic"])
        self._name = name
        self._material_type = material_type
        self._young_modulus = young_modulus
        self._poisson_ratio = poisson_ratio
        self._density = density
        self._tension_allow = tension_allow
        self._compression_allow = compression_allow
        self._shear_allow = shear_allow
        self._yield_allow = yield_allow 
        self._thermExpCoeff = thermExpCoeff

    @property
    def name(self) -> str:
        return self._name

    @property
    def dictionary(self) -> dict:
        """
        return dictionary of material settings
        """
        m_dict = {}
        m_dict["materialType"] = self._material_type
        m_dict["youngModulus"] = self._young_modulus
        m_dict["poissonRatio"] = self._poisson_ratio
        m_dict["density"] = self._density
        m_dict["thermalExpCoeff"] = self._thermExpCoeff
        m_dict["tensionAllow"] = self._tension_allow
        m_dict["compressionAllow"] = self._compression_allow
        m_dict["shearAllow"] = self._shear_allow
        m_dict["yieldAllow"] = self._yield_allow

        # return all items that are not None
        return {k: v for k, v in m_dict.items() if v is not None}
        

class Isotropic(Material):
    def __init__(self, 
        name:str, young_modulus:float, poisson_ratio:float, density:float, tension_allow:float,
        compression_allow:float=None, shear_allow:float=None, yield_allow:float=None, thermExpCoeff:float=None
    ):
        super(Isotropic,self).__init__( 
        name = name,
        material_type = "Isotropic",
        young_modulus = young_modulus,
        poisson_ratio = poisson_ratio,
        density = density,
        tension_allow = tension_allow,
        compression_allow = compression_allow,
        shear_allow = shear_allow,
        yield_allow = yield_allow,
        thermExpCoeff = thermExpCoeff
        )

    @classmethod
    def madeupium(cls, young_modulus=72.0E9, poisson_ratio=0.33, density=2.8E3, tension_allow=20.0E7):
        return cls(name="Madeupium", young_modulus=young_modulus,poisson_ratio=poisson_ratio, density=density, tension_allow=tension_allow)


class Orthotropic(Material):
    # TBD
    pass