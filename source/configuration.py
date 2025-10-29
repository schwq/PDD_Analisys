import json
from photon_DiffMV import PhotonEnergy

class PhotonConfiguration:
    
    def __init__(self, energies: tuple[PhotonEnergy] , numBeams: int, beamWidth: float, angles: tuple[float], photonEnergy: PhotonEnergy):
        self.energies = energies
        self.numBeams = numBeams
        self.beamWidth = beamWidth
        self.angles = angles 
        self.photonEnergy = photonEnergy

class ProtonConfiguration:

    def __init__(self, da : float, db: float, d0: float, numBeams: int, beamWidth: float, angles: tuple[float]):
        self.da = da 
        self.db = db 
        self.d0 = d0
        self.numBeams = numBeams
        self.beamWidth = beamWidth
        self.angles = angles

class SimulationConfiguration:

    def __init__(self, image : str, mask : str, alpha: float, gridSize: int, imageSize: int):
        self.image = image
        self.mask = mask
        self.alpha = alpha
        self.gridSize = gridSize
        self.imageSize = imageSize

def _get_member_enum_from_value(enumclass, value):
    for member in enumclass:
        if member.value == value: return member 
    return None 

class ProgramConfiguration:
    def __init__(self, path : str):
        self.path = path 
        with open(path, 'r') as file: data = json.load(file)
        self.photon = PhotonConfiguration([], 0, 0, [], PhotonEnergy.MV6)
        self.photon.energies = [ _get_member_enum_from_value(PhotonEnergy, int(n)) for n in data['Photon']['PlotEnergies']]        
        self.photon.numBeams = int(data['Photon']['NumBeams'])
        self.photon.beamWidth = float(data['Photon']['BeamWidth'])
        self.photon.photonEnergy = _get_member_enum_from_value(PhotonEnergy, int(data['Photon']['PhotonEnergyPlot']))
        self.photon.angles = [float(n) for n in data['Photon']['Angles']]
        self.proton = ProtonConfiguration(0, 0, 0, 0, 0, [])
        self.proton.da = float(data['Proton']['Da'])
        self.proton.db = float(data['Proton']['Db'])
        self.proton.d0 = float(data['Proton']['D0'])
        self.proton.numBeams = int(data['Proton']['NumBeams'])
        self.proton.beamWidth = float(data['Proton']['BeamWidth'])
        self.proton.angles = [float(n) for n in data['Proton']['Angles']]
        self.sim = SimulationConfiguration("", "", 0, 0, 0)
        self.sim.image = data['Simulation']['ImagePath']
        self.sim.mask = data['Simulation']['MaskPath']
        self.sim.alpha = float(data['Simulation']['ImageAlpha'])
        self.sim.gridSize = int(data['Simulation']['GridSize'])
        self.sim.imageSize = int(data['Simulation']['ImageSize'])

        