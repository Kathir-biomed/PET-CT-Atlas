import nibabel as nib
from datetime import timedelta
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import re
'''
# Hard-coded metadata for patient 1 (can be automated later)
radionuclide_total_dose = 313000000.0  # hard-coded from DICOM
rescale_slope_pet = 1.77409  # Standard_assumption for PET
rescale_intercept_pet = 0.0  # Standard_assumption for PET
rescale_slope_ct = 0.9236  # Standard_assumption for CT
rescale_intercept_ct = -1.6565  # Standard_assumption for CT
weight = 65.0  # hard-coded from DICOM
radionuclide_half_life = 6586.2  # actual half-life of F-18 # hard-coded from DICOM
acquisition_time = timedelta(seconds=84202.000013)  # hard-coded from DICOM
series_time = timedelta(seconds=84202.000000)  # hard-coded from DICOM
series_date = 20010913  # %Y%m%d
acquisition_date = 20010913  # %Y%m%d
patient_id = "PETCT_0117d7f11f"
##########
(0008, 0021): SeriesDate
(0008, 0022): AcquisitionDate
(0008, 0031): SeriesTime
(0008, 0032): AcquisitionTime
(0008, 0060): PT
(0010, 1030): PatientWeight)
(0028, 0051): CorrectedImage (ATTN, DECY)
(0054, 0016): RadiopharmaceuticalInformationSequence
(0018, 1072): RadiopharmaceuticalStartTime
(0018, 1074): RadionuclideTotalDose
(0018, 1075): RadionuclideHalfLife
(0054, 1001): BQML
(0054, 1002): EMISSION
(0054, 1102): START
(0054, 1300): FrameReferenceTime
(0054, 1321): DecayFactor
(0010, 0040): (M,F)
##########
'''
# Hard-coded metadata for patient 2 (can be automated later)
#PET
patient_id = "PETCT_4ef69de4e1" # (0010, 0020) Patient ID PETCT_4ef69de4e1
series_date = 20021025  # %Y%m%d # (0008, 0021) Series Date 20021025
acquisition_date = 20021025  # %Y%m%d # (0008, 0022) Acquisition Date 20021025
series_time = timedelta(seconds=131136.000000)  # (0008, 0031) Series Time 131136.000000
acquisition_time = timedelta(seconds=131136.000000)  # (0008, 0032) Acquisition Time 131136.000000 
weight = 67  # (0010, 1030) Patient's Weight 67
radionuclide_total_dose = 329000000.0  # (0018, 1074) Radionuclide Total Dose             DS: '329000000.0'
radionuclide_half_life = 6586.2  # actual half-life of F-18 # (0018, 1075) Radionuclide Half Life              DS: '6586.2'
rescale_slope_pet = 0.878121  # (0028, 1053) Rescale Slope 0.878121
rescale_intercept_pet = 0  # (0028, 1052) Rescale Intercept 0
#CT
#rescale_slope_ct = 1  # (0028, 1053) Rescale Slope 1
#rescale_intercept_ct = -1024  # (0028, 1052) Rescale Intercept -1024
rescale_slope_ct = 0.9236  # (0019, 1092) [Osteo Regression Line Slope] 0.9236
rescale_intercept_ct = -1.6565  # (0019, 1093) [Osteo Regression Line Intercept] -1.6565

def calculate_hu_metrics(ct_nifti_path, segmented_nifti_path):
    # Load CT NIfTI file
    ct_img = nib.load(ct_nifti_path)
    ct_data = ct_img.get_fdata()

    # Load segmented NIfTI file for HU calculation and mask metrics
    segmented_img = nib.load(segmented_nifti_path)
    segmented_data = segmented_img.get_fdata()

    # Calculate HU values for the organ region
    hu_values_organ = (ct_data * segmented_data * rescale_slope_ct) + rescale_intercept_ct

    # Initialize variables with default values
    hu_mean_organ, hu_std_organ, mask_volume_ml = np.nan, np.nan, np.nan

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

    # Calculate mask volume in milliliters
    voxel_volume_mm3 = np.prod(ct_img.header.get_zooms())  # Voxel volume in mm^3
    mask_volume_ml = np.sum(segmented_data) * voxel_volume_mm3 / 1000.0  # Convert to milliliters

    print("Mask Volume (ml):")
    print(mask_volume_ml)

    return hu_mean_organ, hu_std_organ, mask_volume_ml
'''
def calculate_organ_suv(nifti_file_path, segmented_nifti_path, weight, radionuclide_half_life, series_time, acquisition_time):
    # Load PET NIfTI file
    pet_img = nib.load(nifti_file_path)
    pet_data = pet_img.get_fdata()

    # Load segmented NIfTI file
    segmented_img = nib.load(segmented_nifti_path)
    segmented_data = segmented_img.get_fdata()

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
    return mean_suvbw_organ, suv_peak_organ, suv_max_organ
'''
def calculate_organ_suv(nifti_file_path, segmented_nifti_path, weight, radionuclide_half_life, series_time, acquisition_time):
    # Load PET NIfTI file
    pet_img = nib.load(nifti_file_path)
    pet_data = pet_img.get_fdata()

    # Load segmented NIfTI file
    segmented_img = nib.load(segmented_nifti_path)
    segmented_data = segmented_img.get_fdata()

    # Calculate decayed dose
    decay_time = (acquisition_time - series_time).total_seconds()
    decayed_dose = radionuclide_total_dose * (2 ** (-decay_time / radionuclide_half_life))

    # Calculate SUVbw scale factor
    suvbw_scale_factor = (weight * 1000) / decayed_dose

    # Calculate SUVbw only for the organ region
    suvbw_organ = (pet_data * segmented_data * rescale_slope_pet + rescale_intercept_pet) * suvbw_scale_factor

    # Initialize variables with default values
    mean_suvbw_organ, suv_peak_organ, suv_max_organ = np.nan, np.nan, np.nan

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

    return mean_suvbw_organ, suv_peak_organ, suv_max_organ

# NIfTI file paths
nifti_file_path_pet = "C:/Personal/ABX/patient_2/NIFTI/PETCT_4ef69de4e1/10-25-2002-NA-PET-CT Teilkoerper  primaer mit KM-18049/PET.nii.gz"
nifti_file_path_ct = "C:/Personal/ABX/patient_2/NIFTI/PETCT_4ef69de4e1/10-25-2002-NA-PET-CT Teilkoerper  primaer mit KM-18049/CTres.nii.gz"
existing_excel_file_path = f"C:/Personal/ABX/patient_2/Results/{patient_id}_12_12_23_results.xlsx"

# Get a list of all NIfTI files in the "segmented_organs" folder
segmented_organs_folder = "C:/Personal/ABX/patient_2/segmentations"
segmented_nifti_paths = [os.path.join(segmented_organs_folder, file) for file in os.listdir(segmented_organs_folder) if file.endswith(".nii.gz")]

# Initialize an empty DataFrame
existing_df = pd.DataFrame(columns=[
    'Patient ID', 'Organ', 'HU Mean', 'HU Std Dev', 'SUV Mean', 'SUV Peak', 'SUV Max', 'Mask Volume (ml)'
])


# Iterate through segmented organ NIfTI files
for segmented_nifti_path in segmented_nifti_paths:
    # Calculate HU metrics
    hu_mean, hu_std, mask_volume_ml = calculate_hu_metrics(nifti_file_path_ct, segmented_nifti_path)

    # Calculate SUV metrics (using PET for SUV)
    pet_weight = 65.0  # Assuming a constant weight for simplicity
    pet_radionuclide_half_life = 6586.2  # Assuming a constant half-life for simplicity
    pet_mean_suvbw, pet_suv_peak, pet_suv_max = calculate_organ_suv(nifti_file_path_pet, segmented_nifti_path, pet_weight, pet_radionuclide_half_life, series_time, acquisition_time)

    # Extract organ name from the segmented nifti file name
    organ_name = os.path.splitext(os.path.splitext(os.path.basename(segmented_nifti_path))[0])[0]

    # Extract patient ID from nifti_file_path
    match = re.search(r'/NIFTI/(.*?)/', nifti_file_path_ct)
    if match:
        patient_id = match.group(1)
    else:
        print("Unable to extract Patient ID from nifti_file_path.")

    # Determine the index of the next empty row
    next_index = len(existing_df) + 1

    # Append the new results to the existing DataFrame with the calculated index
    new_data = {
        'Patient ID': [patient_id],
        'Organ': [organ_name],
        'HU Mean': [hu_mean],
        'HU Std Dev': [hu_std],
        'SUV Mean': [pet_mean_suvbw],
        'SUV Peak': [pet_suv_peak],
        'SUV Max': [pet_suv_max],
        'Mask Volume (ml)': [mask_volume_ml],
    }

    new_df = pd.DataFrame(new_data)
    new_df.index = [next_index]

    existing_df = pd.concat([existing_df, new_df], ignore_index=False)

# Save the updated DataFrame to the Excel file without the index column
existing_df.to_excel(existing_excel_file_path, index=False)

# Print a message indicating successful export
print(f"HU and SUV Results appended to {existing_excel_file_path}")

