from proton_SOBP import Plot_Proton_SOBP
from photon_DiffMV import PhotonEnergy, Plot_PhotonBeam
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.ndimage import gaussian_filter
from photon_DiffMV import _Photon_ExpFunc
from skimage.transform import resize
import matplotlib.image as mpimg

class Simulation:

    @staticmethod
    def overlay_on_mri(dose_map, mri_path, alpha=0.4):
        """
        Overlay a 2D dose map on top of an MRI image.
        """
        try:
            mri = mpimg.imread(mri_path)
        except FileNotFoundError:
            print(f"[!] MRI image not found: {mri_path}")
            return

        # Convert to grayscale if needed
        if mri.ndim == 3:
            mri_gray = np.mean(mri, axis=2)
        else:
            mri_gray = mri

        # Normalize and resize MRI to match dose shape
        mri_gray = mri_gray / np.max(mri_gray)
        mri_resized = resize(mri_gray, dose_map.shape, mode='reflect', anti_aliasing=True)

        # Overlay plot
        plt.figure(figsize=(9, 9))
        plt.imshow(mri_resized, cmap='gray', extent=[-12.5, 12.5, -12.5, 12.5])
        plt.imshow(dose_map, cmap='inferno', alpha=alpha, extent=[-12.5, 12.5, -12.5, 12.5])
        plt.colorbar(label="Relative Dose (%)")
        plt.title("Dose Overlay on MRI (Proton or Photon)")
        plt.axis('off')

    @staticmethod
    def simulate_photon_beams_2d(grid_size=300, tumor_radius=3, num_beams=4, beam_width=4, angles=[]):
        # Grid parameters (cm)
        size = 25  # 25x25 cm area
        x = np.linspace(-size/2, size/2, grid_size)
        y = np.linspace(-size/2, size/2, grid_size)
        X, Y = np.meshgrid(x, y)
        
        # Tumor mask
        tumor_mask = X**2 + Y**2 <= tumor_radius**2
        
        # Photon beam model (example: 6 MV)
        a, b = 2.2, 0.06  
        z = np.linspace(0, size, grid_size)
        D, D_norm, _ = _Photon_ExpFunc(a, b, z)

        total_dose = np.zeros_like(X)

        for i in range(num_beams):
            # Beam rotation
            angle = angles[i]
            theta = np.deg2rad(angle)

            # Rotate coordinates for beam direction
            Xr = X*np.cos(theta) - Y*np.sin(theta)
            Yr = X*np.sin(theta) + Y*np.cos(theta)

            # Depth (distance along beam)
            depth = np.clip(Yr + size/2, 0, size)

            # Dose by depth
            dose_depth = np.interp(depth, z, D_norm)

            # Apply Gaussian lateral falloff (beam width control)
            lateral_profile = np.exp(-(Xr**2) / (2 * (beam_width/2)**2))
            beam_dose = dose_depth * lateral_profile

            total_dose += beam_dose

        # Normalize and limit
        total_dose = 100 * total_dose / total_dose.max()

        # Visualization
        plt.figure(figsize=(9, 9))
        plt.imshow(total_dose, extent=[-size/2, size/2, -size/2, size/2],
                origin='lower', cmap='inferno', vmin=0, vmax=100)
        plt.colorbar(label="Relative Dose (%)")

        # Tumor contour
        plt.contour(X, Y, tumor_mask, colors='red', linewidths=2)
        plt.title(f"Photon Beam Dose Distribution ({num_beams} beams, width={beam_width} cm)")
        plt.xlabel("x (cm)")
        plt.ylabel("y (cm)")

    @staticmethod
    def simulate_proton_beams_2d(grid_size=400, tumor_radius=3, num_beams=4, beam_width=2, da=10, db=15, d0=1, angles=[]):
        # Grid setup (cm)
        size = 25
        x = np.linspace(-size/2, size/2, grid_size)
        y = np.linspace(-size/2, size/2, grid_size) 
        X, Y = np.meshgrid(x, y)

        # Tumor mask (center)
        tumor_mask = X**2 + Y**2 <= tumor_radius**2

        # Proton SOBP curve (depth-dose)
        depth, _, _, dose_profile = Plot_Proton_SOBP(da, db, d0)
        dose_profile /= dose_profile.max()  # Normalize to 1

        total_dose = np.zeros_like(X)

        for i in range(num_beams):
            # Beam direction (rotating around tumor)
            angle = angles[i]
            theta = np.deg2rad(angle)

            # Rotate coordinates so beam goes along +Y' axis
            Xr = X*np.cos(theta) - Y*np.sin(theta)
            Yr = X*np.sin(theta) + Y*np.cos(theta)

            # Depth is the distance along beam direction (clipped to 0–max)
            depth_along_beam = np.clip(Yr + size/2, 0, size)

            # Interpolate proton dose at this depth
            dose_depth = np.interp(depth_along_beam, depth, dose_profile)

            # Gaussian lateral falloff (narrower than photons)
            lateral_profile = np.exp(-(Xr**2) / (2 * (beam_width/2)**2))

            # Beam dose
            beam_dose = dose_depth * lateral_profile

            total_dose += beam_dose

        # Normalize to 100%
        total_dose = gaussian_filter(total_dose, sigma=3)
        total_dose = 100 * total_dose / total_dose.max()

        # --- Visualization ---
        plt.figure(figsize=(9, 9))
        plt.imshow(total_dose, extent=[-size/2, size/2, -size/2, size/2],
                origin='lower', cmap='inferno', vmin=0, vmax=100)
        plt.colorbar(label="Relative Dose (%)")

        # Tumor outline
        plt.contour(X, Y, tumor_mask, colors='red', linewidths=2)
        plt.title(f"Proton Beam Dose Distribution ({num_beams} beams, width={beam_width} cm)")
        plt.xlabel("x (cm)")
        plt.ylabel("y (cm)")

        return total_dose

    @staticmethod
    def start():
        photon_dose = Simulation.simulate_photon_beams_2d(grid_size=400, tumor_radius=2, num_beams=4, beam_width=3, angles=[0, 30, 60, 90])
        proton_dose = Simulation.simulate_proton_beams_2d(grid_size=400, tumor_radius=2, num_beams=4, beam_width=3, angles=[0, 30, 60, 90])
        Simulation.overlay_on_mri(proton_dose, "./data/mri.png", alpha=0.8)
        plt.show()
