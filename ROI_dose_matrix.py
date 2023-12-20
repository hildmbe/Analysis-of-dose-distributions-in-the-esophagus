''' Funksjonen GetRoiDoseMatrix() tar inn pasientinfo og valgt ROI. Den returnerer en dose-matrise for valgt ROI.
Denne matrisen har samme dimensjoner som den totale dosematrisen, men forskjellen er at alle elementene utenfor
valgt ROI er satt til 0.
'''

from connect import *
 
def GetRoiDoseMatrix(plan, ROI):
   dose = plan.TreatmentCourse.TotalDose
   dose_matrix = dose.DoseValues.DoseData
   dose_matrix_1d = dose_matrix.flatten()
   dim = dose_matrix.shape

   ROI_dg = dose.GetDoseGridRoi(RoiName=ROI)
   voxel_indices = ROI_dg.RoiVolumeDistribution.VoxelIndices

   ROI_dose_matrix_1d = np.zeros(len(dose_matrix_1d))

   for elem in voxel_indices:
       ROI_dose_matrix_1d[elem] = dose_matrix_1d[elem]

   ROI_dose_matrix = ROI_dose_matrix_1d.reshape(dim)

   return ROI_dose_matrix

   



