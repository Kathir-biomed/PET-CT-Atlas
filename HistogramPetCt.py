import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from datetime import timedelta
import os

# Constants and paths
ct_nifti_path = "C:/Personal/ABX/patient_1/NIFTI/PETCT_0117d7f11f/09-13-2001-NA-PET-CT Ganzkoerper  primaer mit KM-68547/CTres.nii.gz"
pet_nifti_path = "C:/Personal/ABX/patient_1/NIFTI/PETCT_0117d7f11f/09-13-2001-NA-PET-CT Ganzkoerper  primaer mit KM-68547/PET.nii.gz"
segmented_organs_folder = "C:/Personal/ABX/patient_1/Segmented_organs"  

# Hard-coded metadata for patient 1 (can be automated later)
radionuclide_total_dose = 313000000.0  # hard-coded from DICOM
rescale_slope_pet = 1.77409  # Standard_assumption for PET
rescale_intercept_pet = 0.0  # Standard_assumption for PET
rescale_slope_ct = 1  # Standard_assumption for CT
rescale_intercept_ct = -1024  # Standard_assumption for CT
patient_weight = 65.0  # hard-coded from DICOM
radionuclide_half_life = 6586.2  # actual half-life of F-18 # hard-coded from DICOM
acquisition_time = timedelta(seconds=84202.000013)  # hard-coded from DICOM
series_time = timedelta(seconds=84202.000000)  # hard-coded from DICOM
#series_date = 20010913  # %Y%m%d
#acquisition_date = 20010913  # %Y%m%d
patient_id = "PETCT_0117d7f11f"


def plot_histograms(ct_nifti_path, pet_nifti_path, segmented_nifti_path):
    # Load CT and PET NIfTI files
    ct_img = nib.load(ct_nifti_path)
    ct_data = ct_img.get_fdata()
    pet_img = nib.load(pet_nifti_path)
    pet_data = pet_img.get_fdata()

    # Load segmented NIfTI file
    segmented_img = nib.load(segmented_nifti_path)
    segmented_data = segmented_img.get_fdata()

    # CT Histogram
    hu_values = (ct_data * segmented_data * rescale_slope_ct) + rescale_intercept_ct
    valid_ct_data_points = hu_values[segmented_data > 0]
    plt.figure(figsize=(10, 6))
    plt.hist(valid_ct_data_points.flatten(), bins=50, color='blue', alpha=0.7)
    plt.title(f'CT Histogram for {os.path.basename(segmented_nifti_path)}')
    plt.xlabel('CT Number (HU)')
    plt.ylabel('Number of Pixels')
    plt.grid(True)
    plt.show()

    # PET Histogram
    decay_time = acquisition_time - series_time
    decay_time_seconds = decay_time.total_seconds()  # Convert timedelta to seconds
    decayed_dose = radionuclide_total_dose * (2 ** (-decay_time_seconds / radionuclide_half_life))
    suvbw_scale_factor = (patient_weight * 1000) / decayed_dose
    suvbw_values = (pet_data * segmented_data * rescale_slope_pet + rescale_intercept_pet) * suvbw_scale_factor
    valid_pet_data_points = suvbw_values[segmented_data > 0]
    plt.figure(figsize=(10, 6))
    plt.hist(valid_pet_data_points.flatten(), bins=50, color='green', alpha=0.7)
    plt.title(f'PET Histogram for {os.path.basename(segmented_nifti_path)}')
    plt.xlabel('SUVbw')
    plt.ylabel('Number of Pixels')
    plt.grid(True)
    plt.show()

# Iterating over segmented organs and plotting histograms
segmented_nifti_paths = [os.path.join(segmented_organs_folder, f) for f in os.listdir(segmented_organs_folder) if f.endswith('.nii.gz')]
for segmented_nifti_path in segmented_nifti_paths:
    plot_histograms(ct_nifti_path, pet_nifti_path, segmented_nifti_path)
