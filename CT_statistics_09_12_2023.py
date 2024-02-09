import nibabel as nib
import numpy as np
import pandas as pd
import os
import re

# Hard-coded metadata for patient 1 (can be automated later)
#rescale_slope = 1.000000
#rescale_intercept = -1024.000000

rescale_slope = 0.9236
rescale_intercept = -1.6565
patient_id = "PETCT_0117d7f11f"

def calculate_hu_metrics(ct_nifti_path, segmented_nifti_path):
    # Load CT NIfTI file
    ct_img = nib.load(ct_nifti_path)
    ct_data = ct_img.get_fdata()

    # Load segmented NIfTI file for HU calculation and mask metrics
    segmented_img = nib.load(segmented_nifti_path)
    segmented_data = segmented_img.get_fdata()

    # Calculate HU values for the organ region
    hu_values_organ = (ct_data * segmented_data * rescale_slope) + rescale_intercept

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

# NIfTI file paths
nifti_file_path = "C:/Personal/ABX/patient_1/NIFTI/PETCT_0117d7f11f/09-13-2001-NA-PET-CT Ganzkoerper  primaer mit KM-68547/CTres.nii.gz"


# Get a list of all NIfTI files in the "segmented_organs" folder
segmented_organs_folder = "C:/Personal/ABX/patient_1/segmented_organs"
segmented_nifti_paths = [os.path.join(segmented_organs_folder, file) for file in os.listdir(segmented_organs_folder) if file.endswith(".nii.gz")]

# Iterate through segmented organ NIfTI files
for segmented_nifti_path in segmented_nifti_paths:
    # Calculate HU mean
    hu_mean = calculate_hu_metrics(nifti_file_path, segmented_nifti_path)

    # Extract organ name from the segmented nifti file name
    organ_name = os.path.splitext(os.path.splitext(os.path.basename(segmented_nifti_path))[0])[0]

    # Extract patient ID from nifti_file_path
    match = re.search(r'/NIFTI/(.*?)/', nifti_file_path)
    if match:
        patient_id = match.group(1)
    else:
        print("Unable to extract Patient ID from nifti_file_path.")

    # Print the organ name and HU mean
    print(f"Organ: {organ_name}, Mean HU: {hu_mean}")

'''
# Sample segmented organ NIfTI file paths
segmented_nifti_paths = [
    "C:/Personal/ABX/patient_1/segmentations/stomach.nii.gz",
    "C:/Personal/ABX/patient_1/segmentations/brain.nii.gz",
    "C:/Personal/ABX/patient_1/segmentations/femur_left.nii.gz",
    "C:/Personal/ABX/patient_1/segmentations/femur_right.nii.gz",
    # Add more paths for other segmented organs
]
'''
# Load the existing DataFrame from the Excel file
existing_excel_file_path = f"C:/Personal/ABX/patient_1/Results/{patient_id}_CTnewwresults.xlsx"

if os.path.exists(existing_excel_file_path):
    existing_df = pd.read_excel(existing_excel_file_path)
else:
    # Create an empty DataFrame if the file doesn't exist
    existing_df = pd.DataFrame()

# Iterate through segmented organ NIfTI files
for segmented_nifti_path in segmented_nifti_paths:
    # Calculate HU metrics
    hu_mean, hu_std = calculate_hu_metrics(nifti_file_path, segmented_nifti_path)

    # Extract organ name from the segmented nifti file name
    organ_name = os.path.splitext(os.path.splitext(os.path.basename(segmented_nifti_path))[0])[0]

    # Extract patient ID from nifti_file_path
    match = re.search(r'/NIFTI/(.*?)/', nifti_file_path)
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
    }

    new_df = pd.DataFrame(new_data)
    new_df.index = [next_index]

    existing_df = pd.concat([existing_df, new_df], ignore_index=False)

# Save the updated DataFrame to the Excel file without the index column
existing_df.to_excel(existing_excel_file_path, index=False)

# Print a message indicating successful export
print(f"HU Results appended to {existing_excel_file_path}")
