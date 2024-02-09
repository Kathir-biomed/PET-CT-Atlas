import nibabel as nib
from datetime import timedelta
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import re

# hard-coded the necessary metadata for patient 1 from DICOM file as I'm working with NIfTI file. Later the hard-coded values can be automated.
radionuclide_total_dose = 313000000.0  # hard-coded from DICOM
rescale_slope = 1.0  # Standard_assumption
rescale_intercept = 0.0  # Standard_assumption
weight = 65.0  # hard-coded from DICOM
radionuclide_half_life = 6586.2  # actual half-life of F-18 # hard-coded from DICOM
acquisition_time = timedelta(seconds=84202.000013)  # hard-coded from DICOM
series_time = timedelta(seconds=84202.000000)  # hard-coded from DICOM
series_date = 20010913  # %Y%m%d
acquisition_date = 20010913  # %Y%m%d
patient_id= "PETCT_0117d7f11f"


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
    suvbw_organ = (pet_data * segmented_data * rescale_slope + rescale_intercept) * suvbw_scale_factor

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
    # ToPlott histogram of SUV values in the segmented organ region
    plt.hist(valid_data_points, bins=50, color='blue', alpha=0.7)
    plt.xlabel('SUV Value')
    plt.ylabel('Frequency')
    plt.title('Histogram of SUV Values in Segmented Organ Region')
    plt.show()
    '''

def calculate_hu_and_mask_metrics(segmented_nifti_path):
    # Load segmented NIfTI file for HU calculation and mask metrics
    ct_img = nib.load(segmented_nifti_path)  # Assuming segmented NIfTI is a CT dataset
    ct_data = ct_img.get_fdata()
    segmented_data = ct_data > 0  # Assuming values > 0 indicate the organ region

    # Calculate HU mean and standard deviation for the organ region
    hu_values_organ = ct_data[segmented_data]
    
    # Check if there are valid data points in the segmented region
    if hu_values_organ.size > 0:
        hu_mean_organ = hu_values_organ.mean()
        hu_std_organ = hu_values_organ.std()
        
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

nifti_file_path = "C:/Personal/ABX/patient_1/NIFTI/PETCT_0117d7f11f/09-13-2001-NA-PET-CT Ganzkoerper  primaer mit KM-68547/PET.nii.gz"
#segmented_nifti_path = "C:/Personal/ABX/patient_1/Segmented_organs/liver.nii.gz"
#segmented_nifti_path = "C:/Personal/ABX/patient_1/segmentations/brain.nii.gz"
#calculate_hu_and_mask_metrics(segmented_nifti_path)
# Load the existing DataFrame from the Excel file
existing_excel_file_path = f"C:/Personal/ABX/patient_1/Results/{patient_id}_results.xlsx"

if os.path.exists(existing_excel_file_path):
    existing_df = pd.read_excel(existing_excel_file_path)
else:
    # Tocreate an empty DataFrame if the file doesn't exist
    existing_df = pd.DataFrame()

# Determine the index of the next empty row
next_index = len(existing_df) + 1

# Calculate SUV and HU metrics for the new segmented NIfTI
new_segmented_nifti_path = "C:/Personal/ABX/patient_1/segmentations/vertebrae_T12.nii.gz"
new_mean_suvbw, new_suv_peak, new_suv_max = calculate_organ_suv(nifti_file_path, new_segmented_nifti_path, weight, radionuclide_half_life, series_time, acquisition_time)
new_hu_mean, new_hu_std, new_mask_volume = calculate_hu_and_mask_metrics(new_segmented_nifti_path)

# Extract organ name from the new segmented nifti file name
#new_organ_name = os.path.splitext(os.path.basename(new_segmented_nifti_path))[0].split("_organ")[0].split("_")[-1][:-4]
#new_organ_name = os.path.splitext(os.path.basename(new_segmented_nifti_path))[0]
new_organ_name = os.path.splitext(os.path.splitext(os.path.basename(new_segmented_nifti_path))[0])[0]

# Extract patient ID from nifti_file_path
match = re.search(r'/NIFTI/(.*?)/', nifti_file_path)
if match:
    patient_id = match.group(1)
else:
    print("Unable to extract Patient ID from nifti_file_path.")

# Append the new results to the existing DataFrame with the calculated index
new_data = {
    'Patient ID': [patient_id],
    'Organ': [new_organ_name],
    'SUV Mean': [new_mean_suvbw],
    'SUV Peak': [new_suv_peak],
    'SUV Max': [new_suv_max],
    'HU Mean': [new_hu_mean],
    'HU Std Dev': [new_hu_std],
    'Mask Volume (ml)': [new_mask_volume]
}

new_df = pd.DataFrame(new_data)
new_df.index = [next_index]

existing_df = pd.concat([existing_df, new_df], ignore_index=False)

# Save the updated DataFrame to the Excel file without the index column
existing_df.to_excel(existing_excel_file_path, index=False)

# Print a message indicating successful export
print(f"Results appended to {existing_excel_file_path}")


'''
# Save the updated DataFrame to the Excel file
existing_df.to_excel(existing_excel_file_path, index=True)

# Print a message indicating successful export
print(f"Results appended to {existing_excel_file_path}")
#############
'''
