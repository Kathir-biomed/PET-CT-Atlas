import nibabel as nib
import numpy as np
import pandas as pd
import os
import re
from datetime import timedelta

# Hard-coded metadata for patient 1 (can be automated later)
radionuclide_total_dose = 313000000.0  # hard-coded from DICOM
rescale_slope_pet = 1.77409  # Standard_assumption for PET
rescale_intercept_pet = 0.0  # Standard_assumption for PET
rescale_slope_ct = 1  # Standard_assumption for CT
rescale_intercept_ct = -1024  # Standard_assumption for CT
weight = 65.0  # hard-coded from DICOM
radionuclide_half_life = 6586.2  # actual half-life of F-18 # hard-coded from DICOM
acquisition_time = timedelta(seconds=84202.000013)  # hard-coded from DICOM
series_time = timedelta(seconds=84202.000000)  # hard-coded from DICOM
series_date = 20010913  # %Y%m%d
acquisition_date = 20010913  # %Y%m%d
patient_id = "PETCT_0117d7f11f"

def calculate_hu_metrics(ct_data, segmented_data):
    # Calculate HU values for the organ region
    hu_values_organ = (ct_data * segmented_data * rescale_slope_ct) + rescale_intercept_ct

    # Check if there are valid data points in the segmented region
    valid_data_points = hu_values_organ[segmented_data > 0]

    if valid_data_points.size > 0:
        # Calculate HU mean and standard deviation for the organ
        hu_mean_organ = valid_data_points.mean()
        hu_std_organ = valid_data_points.std()

        print("HU Mean:")
        print(hu_mean_organ)

        print("HU Standard Deviation:")
        print(hu_std_organ)

    else:
        print("No valid data points in the segmented organ region. Unable to calculate HU metrics.")
        hu_mean_organ, hu_std_organ = np.nan, np.nan

    return hu_mean_organ, hu_std_organ

def calculate_suv_metrics(pet_data, segmented_data):
    # Calculate decayed dose
    decay_time = (acquisition_time - series_time).total_seconds()
    decayed_dose = radionuclide_total_dose * (2 ** (-decay_time / radionuclide_half_life))

    # Calculate SUVbw scale factor
    suvbw_scale_factor = (weight * 1000) / decayed_dose

    # Calculate SUVbw only for the organ region
    suvbw_organ = (pet_data * segmented_data * rescale_slope_pet + rescale_intercept_pet) * suvbw_scale_factor

    # Check if there are valid data points in the segmented region
    valid_data_points = suvbw_organ[segmented_data > 0]

    if valid_data_points.size > 0:
        # Calculate the mean SUVbw value for the organ
        mean_suvbw_organ = valid_data_points.mean()
        print("Mean SUVbw:")
        print(mean_suvbw_organ)

        # Calculate SUV Max for the organ
        suv_max_organ = valid_data_points.max()
        print("SUV Max:")
        print(suv_max_organ)

        # Find the spatial coordinates (indices) of the maximum SUV voxel
        max_suv_index = np.unravel_index(np.argmax(suvbw_organ), suvbw_organ.shape)

        # Define a small region (1 cm^3 sphere) around the maximum SUV voxel
        sphere_radius = int(np.ceil(1 / pet_img.header.get_zooms()[0]))  # in voxel units
        region_around_max_suv = suvbw_organ[
            max(0, max_suv_index[0] - sphere_radius):min(suvbw_organ.shape[0], max_suv_index[0] + sphere_radius + 1),
            max(0, max_suv_index[1] - sphere_radius):min(suvbw_organ.shape[1], max_suv_index[1] + sphere_radius + 1),
            max(0, max_suv_index[2] - sphere_radius):min(suvbw_organ.shape[2], max_suv_index[2] + sphere_radius + 1)
        ]

        # Calculate SUV Peak for the organ (average value within the 1-cm^3 sphere)
        suv_peak_organ = region_around_max_suv.mean()
        print("SUV Peak:")
        print(suv_peak_organ)

    else:
        print("No valid data points in the segmented organ region. Unable to calculate SUV metrics.")
        mean_suvbw_organ, suv_peak_organ, suv_max_organ = np.nan, np.nan, np.nan

    return mean_suvbw_organ, suv_peak_organ, suv_max_organ

def calculate_hu_and_suv_metrics(nifti_file_path, segmented_nifti_paths):
    # Load PET/CT NIfTI file
    img = nib.load(nifti_file_path)
    data = img.get_fdata()

    # Iterate through segmented organ NIfTI files
    for segmented_nifti_path in segmented_nifti_paths:
        # Load segmented NIfTI file
        segmented_img = nib.load(segmented_nifti_path)
        segmented_data = segmented_img.get_fdata()

        # Extract organ name from the segmented nifti file name
        organ_name = os.path.splitext(os.path.splitext(os.path.basename(segmented_nifti_path))[0])[0]

        print(f"\nCalculating metrics for organ: {organ_name}")

        # Calculate HU metrics
        if "CT" in nifti_file_path:
            hu_mean, hu_std = calculate_hu_metrics(data, segmented_data)
        else:
            print("Invalid file type. HU metrics can only be calculated for CT images.")

        # Calculate SUV metrics
        if "PET" in nifti_file_path:
            suv_mean, suv_peak, suv_max = calculate_suv_metrics(data, segmented_data)
        else:
            print("Invalid file type. SUV metrics can only be calculated for PET images.")

# NIfTI file paths
nifti_file_path_pet = "C:/Personal/ABX/patient_1/NIFTI/PETCT_0117d7f11f/09-13-2001-NA-PET-CT Ganzkoerper  primaer mit KM-68547/PET.nii.gz"
nifti_file_path_ct = "C:/Personal/ABX/patient_1/NIFTI/PETCT_0117d7f11f/09-13-2001-NA-PET-CT Ganzkoerper  primaer mit KM-68547/CTres.nii.gz"

# ... (previous code)

# NIfTI file paths
nifti_file_path_pet = "C:/Personal/ABX/patient_1/NIFTI/PETCT_0117d7f11f/09-13-2001-NA-PET-CT Ganzkoerper  primaer mit KM-68547/PET.nii.gz"
nifti_file_path_ct = "C:/Personal/ABX/patient_1/NIFTI/PETCT_0117d7f11f/09-13-2001-NA-PET-CT Ganzkoerper  primaer mit KM-68547/CTres.nii.gz"

# Organ directory
organ_directory = "C:/Personal/ABX/patient_1/segmentations"

# Discover all organ NIfTI files in the directory
segmented_nifti_paths = [os.path.join(organ_directory, file) for file in os.listdir(organ_directory) if file.endswith(".nii.gz")]

# Calculate HU and SUV metrics for each organ
calculate_hu_and_suv_metrics(nifti_file_path_ct, segmented_nifti_paths)
calculate_hu_and_suv_metrics(nifti_file_path_pet, segmented_nifti_paths)
