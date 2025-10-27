from proton_SOBP import Plot_Proton_SOBP, get_proton_SOBP_data
from photon_DiffMV import PhotonEnergy, Plot_PhotonBeam
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.ndimage import gaussian_filter
from photon_DiffMV import EMPDDMV_PhotonBeam, PhotonsEnergiesParameters
from skimage.transform import resize
import matplotlib.image as mpimg
from skimage import io, color, filters, measure, morphology
import random

class Image:
    def __init__(self, data, alpha, mask):
        self.data = data
        self.alpha = alpha
        self.mask = mask 

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

import matplotlib.pyplot as plt
import numpy as np

def select_cancer_region(image_path, size):

    img = mpimg.imread(image_path)
    if img.ndim == 3:
        img_gray = np.mean(img, axis=2)
    else:
        img_gray = img

    fig, ax = plt.subplots(figsize=(9, 9))
    img_gray = np.rot90(img_gray, 2)

    ax.imshow(img_gray, cmap='gray', origin='lower', extent=[-size/2, size/2, -size/2, size/2])
    ax.set_title("Click tumor center, then click tumor border")

    pts = []
    circle_artist = [None]

    def onclick(event):
        if event.inaxes != ax:
            return
        if len(pts) == 0:
            pts.append((event.xdata, event.ydata))
            ax.plot(event.xdata, event.ydata, 'ro')
            fig.canvas.draw()
        elif len(pts) == 1:
            # Second click defines border
            pts.append((event.xdata, event.ydata))
            if circle_artist[0]:
                circle_artist[0].remove()
            dx = pts[1][0] - pts[0][0]
            dy = pts[1][1] - pts[0][1]
            r = np.hypot(dx, dy)
            circle = plt.Circle(pts[0], r, color='r', fill=False, lw=2)
            ax.add_patch(circle)
            circle_artist[0] = circle
            fig.canvas.draw()
            plt.close(fig)

    def onmove(event):
        if len(pts) == 1 and event.inaxes == ax:
            if circle_artist[0]:
                circle_artist[0].remove()
            dx = event.xdata - pts[0][0]
            dy = event.ydata - pts[0][1]
            r = np.hypot(dx, dy)
            circle = plt.Circle(pts[0], r, color='r', fill=False, lw=1)
            ax.add_patch(circle)
            circle_artist[0] = circle
            fig.canvas.draw_idle()

    fig.canvas.mpl_connect('button_press_event', onclick)
    fig.canvas.mpl_connect('motion_notify_event', onmove)
    plt.show()

    if len(pts) < 2:
        print("Canceled or incomplete selection.")
        return None, None, None

    cx, cy = pts[0]
    bx, by = pts[1]
    r = np.hypot(bx - cx, by - cy)
    return cx, cy, r


def add_contour_levels(x, y ,z, start=10, end=101, step=10, color='white', linewidth=1.5, fontsize=8):
    levels = np.arange(start, end, step)
    contours = plt.contour(x, y, z, levels=levels, colors=color, linewidths=linewidth)
    plt.clabel(contours, inline=True, fontsize=fontsize, fmt='%d%%')

def info_mean_values(total_dose, tumor_mask, image):
    mean_dose_tumor = np.mean(total_dose[tumor_mask]) 
    mean_dose_all = np.mean(total_dose)
    mean_dose_not_tumor = np.mean(total_dose[image.mask & ~tumor_mask])
    print(f"Mean dose in tumor tissue = {mean_dose_tumor:.2f}%")
    print(f"Mean dose in all tissue = {mean_dose_all:.2f}%")
    print(f"Mean dose in not tumor tissue = {mean_dose_not_tumor:.2f}%")

def add_contour_tumor(x, y, tx, ty, tumor_mask, radius, color='cyan', linewidth=2):
    plt.contour(x, y, tumor_mask, colors=color, linewidths=linewidth)
    plt.text(tx + radius + 0.2, ty, "Tumor tissue", color=color)

def create_mri_image(path, alpha=0.4):
    try:
        mri = mpimg.imread(path)
    except FileNotFoundError:
        print(f"[!] MRI image not found: {path}")
        return
    # Convert to grayscale if needed
    if mri.ndim == 3:
        mri_gray = np.mean(mri, axis=2)
    else:
        mri_gray = mri
    # Normalize and resize MRI to match dose shape
    mri_gray = mri_gray / np.max(mri_gray)
    return Image(mri_gray, alpha, None)

def apply_red_mask(image : Image, path):
    mask_image = mpimg.imread(path)
    rc = mask_image[:, :, 0]
    gc = mask_image[:, :, 1]
    bc = mask_image[:, :, 2]
    mask = (rc > 0.3) & (gc < 0.2) & (bc < 0.2)
    image.mask = mask
    return image

def create_meshgrid(size=20, grid_size=300):  
    x = np.linspace(-size/2, size/2, grid_size)
    y = np.linspace(-size/2, size/2, grid_size)
    X, Y = np.meshgrid(x, y)
    return x, y, X, Y, size

def create_tumor_mask(X, Y, tx=0, ty=0, radius=3):
    tumor_mask = (X - tx)**2 + (Y - ty)**2 <= radius**2
    not_tumor_mask = ~tumor_mask
    return tumor_mask, not_tumor_mask

def get_photon_data(X, size, grid_size):
    n, u = PhotonsEnergiesParameters[PhotonEnergy.MV6]
    z = np.linspace(0, size, grid_size)
    D, D_norm, _ = EMPDDMV_PhotonBeam(z, n, u)
    total_dose = np.zeros_like(X)
    return z, D, D_norm, total_dose

def compute_tissue_depth(Xr, Yr, mask, size, grid_size):
    """
    Compute effective tissue depth along the beam direction
    (accumulated only inside tissue regions).
    """
    tissue_depth = np.zeros_like(Xr)

    # Rescale mask to match beam grid
    mask_resized = resize(mask.astype(float), Xr.shape, mode='reflect', anti_aliasing=False) > 0.5
    
    # Sort coordinates along the beam axis (Yr direction)
    # We assume Yr increases in beam propagation direction
    sort_idx = np.argsort(Yr, axis=0)
    Yr_sorted = np.take_along_axis(Yr, sort_idx, axis=0)
    mask_sorted = np.take_along_axis(mask_resized, sort_idx, axis=0)
    
    # Increment depth only inside tissue
    dy = size / grid_size
    depth_sorted = np.cumsum(mask_sorted * dy, axis=0)
    
    # Unsort back to original order
    tissue_depth = np.take_along_axis(depth_sorted, np.argsort(sort_idx, axis=0), axis=0)
    
    return tissue_depth


def add_beams(num_beams, angles, X, Y, tx, ty, z, D_norm, beam_width, size, total_dose, image):
    
    for i in range(num_beams):
        
        # Beam rotation
        angle = angles[i]
        theta = np.deg2rad(angle)
        
        xs = X - tx 
        ys = Y - ty
        
        # Rotate coordinates for beam direction
        Xr = xs *np.cos(theta) -  ys *np.sin(theta)
        Yr = xs *np.sin(theta) + ys *np.cos(theta)
        
        if hasattr(image, "mask") and image.mask is not None:
            # Compute tissue depth (only inside head)
            tissue_depth = compute_tissue_depth(Xr, Yr, image.mask, size, X.shape[0])
            depth = np.clip(tissue_depth, 0, size)
        else:
            # fallback if no mask
            depth = np.clip(Yr + size/2, 0, size)

        dose_depth = np.interp(depth, z, D_norm)

        # Apply Gaussian lateral falloff (beam width control)
        lateral_profile = np.exp(-(Xr**2) / (2 * (beam_width/2)**2))
        beam_dose = dose_depth * lateral_profile
        total_dose += beam_dose
    return total_dose

def apply_mask(image : Image, total_dose):
    if image.mask is not None:
        mask_resized = resize(image.mask.astype(float), total_dose.shape, mode='reflect', anti_aliasing=False) > 0.5
        image.mask = mask_resized
        total_dose = total_dose * mask_resized  # zero dose outside mask
    return total_dose, image

def adjust_total_dose(total_dose):
    total_dose = gaussian_filter(total_dose, sigma=10)
    total_dose = 100 * total_dose / total_dose.max()
    return total_dose
    
def show_image_mri(image : Image, total_dose, size):
    print(f"Total dose shape is {total_dose.shape}")
    if image is not None:
        mri_resized = resize(image.data, total_dose.shape, mode='reflect', anti_aliasing=True)
        plt.imshow(mri_resized, cmap='gray', extent=[-size/2, size/2, -size/2, size/2], alpha=image.alpha)

def plot_total_dose(title, total_dose, size, beam_width, num_beams):
    plt.imshow(total_dose, extent=[-size/2, size/2, -size/2, size/2],
            origin='lower', cmap='gnuplot2', vmin=0, vmax=100, alpha=0.5)
    plt.colorbar(label="Relative Dose (%)")
    
    plt.title(f"{title} Dose Distribution ({num_beams} beams, width={beam_width} cm)")
    plt.xlabel("x (cm)")
    plt.ylabel("y (cm)")

def simulate_photon_beams_2d(grid_size=300, tumor_radius=3, num_beams=4, beam_width=4, angles=[], image=None, tx=0, ty=0):
    
    x, y, X, Y, size = create_meshgrid(20, grid_size)
    
    tumor_mask, not_tumor_mask = create_tumor_mask(X, Y, tx, ty, tumor_radius)

    z, D, D_norm, total_dose = get_photon_data(X, size, grid_size)    

    total_dose = add_beams(num_beams, angles, X, Y, tx, ty, z, D_norm, beam_width, size, total_dose, image)

    total_dose, image = apply_mask(image, total_dose)

    total_dose = adjust_total_dose(total_dose)

    info_mean_values(total_dose, tumor_mask, image)
    
    plt.figure(figsize=(9, 9))
    
    show_image_mri(image, total_dose, size)

    add_contour_levels(x, y, total_dose)        
    add_contour_tumor(x, y, tx, ty, tumor_mask, tumor_radius)
    
    plot_total_dose("Photon beam", total_dose, size, beam_width, num_beams)
    return total_dose
    

def get_proton_data(da, db, d0, X):
    depth, _, _, dose_profile = get_proton_SOBP_data(da, db, d0)
    dose_profile /= dose_profile.max()  
    total_dose = np.zeros_like(X)
    return total_dose, dose_profile, depth
    

def simulate_proton_beams_2d(grid_size=400, tumor_radius=3, num_beams=4, beam_width=2, da=3, db=6, d0=1, angles=[], image=None, tx=0, ty=0):
    
    x, y, X, Y, size = create_meshgrid(20, grid_size)

    tumor_mask, not_tumor_mask = create_tumor_mask(X, Y, tx, ty, tumor_radius)

    total_dose, dose_profile, depth = get_proton_data(da, db, d0, X)

    total_dose = add_beams(num_beams, angles, X, Y, tx, ty, depth, dose_profile, beam_width, size, total_dose, image)    
    
    total_dose, image = apply_mask(image, total_dose)
    
    total_dose = adjust_total_dose(total_dose)

    info_mean_values(total_dose, tumor_mask, image)
    
    plt.figure(figsize=(9, 9))
    
    show_image_mri(image, total_dose, size)

    add_contour_levels(X, Y, total_dose)
    add_contour_tumor(x, y, tx, ty, tumor_mask, tumor_radius)
    
    plot_total_dose("Proton beam", total_dose, size, beam_width, num_beams)
    
    return total_dose

def _op_proton_beam_2d(grid_size=400, tumor_radius=3, num_beams=4, beam_width=2, da=3, db=6, d0=1, angles=[], image=None, tx=0, ty=0):
    x, y, X, Y, size = create_meshgrid(20, grid_size)

    tumor_mask, not_tumor_mask = create_tumor_mask(X, Y, tx, ty, tumor_radius)

    total_dose, dose_profile, depth = get_proton_data(da, db, d0, X)

    total_dose = add_beams(num_beams, angles, X, Y, tx, ty, depth, dose_profile, beam_width, size, total_dose, image)    
    
    total_dose, image = apply_mask(image, total_dose)
    
    total_dose = adjust_total_dose(total_dose)

    return total_dose, tumor_mask, image

def _op_photon_beam_2d(grid_size=300, tumor_radius=3, num_beams=4, beam_width=4, angles=[], image=None, tx=0, ty=0):
    
    x, y, X, Y, size = create_meshgrid(20, grid_size)
    
    tumor_mask, not_tumor_mask = create_tumor_mask(X, Y, tx, ty, tumor_radius)

    z, D, D_norm, total_dose = get_photon_data(X, size, grid_size)    

    total_dose = add_beams(num_beams, angles, X, Y, tx, ty, z, D_norm, beam_width, size, total_dose, image)

    total_dose, image = apply_mask(image, total_dose)

    total_dose = adjust_total_dose(total_dose)

    return total_dose, tumor_mask, image
"""
def score_beam(image, num_beams, angles, beam_width, da, db,
                penalty, grid_size, tumor_radius, cx, cy):
    total_dose_photon, tumor_mask_photon,  image = _op_photon_beam_2d(grid_size=grid_size, tumor_radius=tumor_radius, num_beams=num_beams, beam_width=beam_width, angles=angles, image=image, tx=cx, ty=cy)
    total_dose_proton, tumor_mask_proton,  image = _op_proton_beam_2d(grid_size=grid_size, da=da, db=db, d0=1, tumor_radius=tumor_radius, num_beams=num_beams, beam_width=beam_width, angles=angles, image=image, tx=cx, ty=cy)

    not_tumor_mask =  image.mask & ~tumor_mask_proton

    mean_tumor_photon = np.mean(total_dose_photon[tumor_mask_photon])
    mean_not_tumor_photon = np.mean(total_dose_photon[not_tumor_mask])

    mean_tumor_proton = np.mean(total_dose_proton[tumor_mask_proton])
    mean_not_tumor_proton = np.mean(total_dose_proton[not_tumor_mask])

    score_photon = mean_tumor_photon - penalty * mean_not_tumor_photon
    score_proton = mean_tumor_proton - penalty * mean_not_tumor_proton
    return score_photon, score_proton


def optimize_beams(image,grid_size, tumor_radius, cx, cy, generations=200, population_size=10, penalty=1.0):
    best_photon = -np.inf
    best_params_photon = None
    best_proton = -np.inf
    best_params_proton = None 

    for gen in range(generations):
        print(f"Generation {gen+1}/{generations}")
        for i in range(population_size):
            num_beams = random.randint(1, 20)
            beam_width = np.random.uniform(0.1, 8.0)
            da = np.random.uniform(0.0, 20.0)
            db = np.random.uniform(da, 25.0)
            angles = sorted(random.sample(range(0, 180, 5), num_beams))
            photon_score, proton_score = score_beam(image,
                num_beams, angles, beam_width, da, db, 
                penalty, grid_size, tumor_radius, cx, cy)
            if proton_score > best_proton:
                best_proton = proton_score
                best_params_proton = dict(num_beams=num_beams, beam_width=beam_width, da=da, db=db, angles=angles, score=proton_score)
            if photon_score > best_photon:
                best_photon = photon_score
                best_params_photon = dict(num_beams=num_beams, beam_width=beam_width, angles=angles, score=photon_score)
    
    return best_params_photon, best_params_proton
"""

def start():
    
    cx, cy, radius = select_cancer_region("./assets/medulloblastoma.jpg", size=20)
    image = create_mri_image("./assets/medulloblastoma.jpg",alpha=1)
    image = apply_red_mask(image, "./assets/mask.png")
    #best_photon, best_proton = optimize_beams(image, 400, radius, cx, cy)
 
    #simulate_photon_beams_2d(grid_size=400, tumor_radius=radius, num_beams=best_photon["num_beams"], beam_width=best_photon["beam_width"], angles=best_photon["angles"], image=image, tx=cx, ty=cy)
    #simulate_proton_beams_2d(grid_size=400, tumor_radius=radius, num_beams=best_proton["num_beams"], beam_width=best_proton["beam_width"], angles=best_proton["angles"], da=best_proton["da"], db=best_proton["db"], image=image, tx=cx, ty=cy)
    
    simulate_photon_beams_2d(grid_size=400, tumor_radius=radius, num_beams=10, beam_width=2, angles=[0, 10, 20, 30, 40, 50, 60, 70, 80, 90], image=image, tx=cx, ty=cy)
    simulate_proton_beams_2d(grid_size=400, tumor_radius=radius, num_beams=8, beam_width=2, angles=[0, 5, 10, 15, 20, 25, 30, 35], da=3, db=7, image=image, tx=cx, ty=cy)
    
    plt.show()