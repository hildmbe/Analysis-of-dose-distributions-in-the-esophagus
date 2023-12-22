from connect import *
from CreateROIs import createRoiFromDose, createIntersecROI
import os; import csv
 
def GetROIDoseVolumeAndIndices(plan, ROI):
   dose = plan.TreatmentCourse.TotalDose
   dose_matrix = dose.DoseValues.DoseData
   dose_matrix_1d = dose_matrix.flatten()
   ROI_dg = dose.GetDoseGridRoi(RoiName=ROI)
   voxel_indices = ROI_dg.RoiVolumeDistribution.VoxelIndices
   
   relative_volumes = ROI_dg.RoiVolumeDistribution.RelativeVolumes
   total_vol = ROI_dg.RoiVolumeDistribution.TotalVolume
   
   ROI_dose_matrix_1d = np.zeros(len(voxel_indices))
   volumeMatrix_1d = np.zeros(len(voxel_indices))

   for i, elem in enumerate(voxel_indices):
       ROI_dose_matrix_1d[i] = dose_matrix_1d[elem]
       volumeMatrix_1d[i] = relative_volumes[i] * total_vol

   return ROI_dose_matrix_1d, volumeMatrix_1d, voxel_indices
   

def getSlicethicknessAndDim(plan, ROI):
    dose = plan.TreatmentCourse.TotalDose
    dose_matrix = dose.DoseValues.DoseData
    
    dim = dose_matrix.shape
    ROI_dg = dose.GetDoseGridRoi(RoiName=ROI)
    
    slice_thickness = ROI_dg.InDoseGrid.VoxelSize['z']
    voxel_volume = ROI_dg.InDoseGrid.VoxelSize['z']*ROI_dg.InDoseGrid.VoxelSize['x']*ROI_dg.InDoseGrid.VoxelSize['y']

    return dim, slice_thickness, voxel_volume
    
outputdir = 'H:/DoseData'

patient_db = get_current('PatientDB')


for i, p in enumerate(pat_numbers):
    pat_info = patient_db.QueryPatientInfo(Filter={'LastName': p})
    patient = patient_db.LoadPatient(PatientInfo=pat_info[0])
    case = patient.Cases[0]
    case.SetCurrent()
    plan = patient.Cases[0].TreatmentPlans[plans[i]]
    plan.SetCurrent()
    examination = patient.Cases[0].Examinations['CT 1']
        
    with open(outputdir+'/'+'SliceThicknessAndDim.csv', 'a', newline='') as csvfile:

        csvwriter = csv.writer(csvfile)
        #csvwriter.writerow(['Patient, Matrix_dim[0], Matrix_dim[1], Matrix_dim[2], Slice thickness, Voxel Volume'])
        
        dim, slice_thickness, voxel_vol = getSlicethicknessAndDim(plan, 'Esophagus')

        csvwriter.writerow([p, dim[0], dim[1], dim[2], slice_thickness, voxel_vol])
     
    Esophagus_dose_matrix_1d, Esophagus_volume_matrix_1d, Esophagus_voxel_indices = GetROIDoseVolumeAndIndices(plan, 'Esophagus')
    
    os.makedirs('H:/DoseData/pasientnr' + str(p))
    
    Esophagus_dose_matrix_1d.tofile('H:/DoseData/pasientnr' + str(p) +'/dose_matrix_1d', sep = "\r\n")
    Esophagus_volume_matrix_1d.tofile('H:/DoseData/pasientnr' + str(p) +'/volume_matrix_1d', sep = "\r\n")
    Esophagus_voxel_indices.tofile('H:/DoseData/pasientnr' + str(p) +'/voxel_indices', sep = "\r\n")
    
    
    
    

