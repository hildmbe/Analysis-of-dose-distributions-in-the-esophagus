from connect import *
import csv
from CreateROIs import createRoiFromDose, createIntersecROI
from VOL_matrix import GetVolumeDoseMatrix, relativeAreaCovered, GetLengthOfROI
from ROI_dose_matrix import GetRoiDoseMatrix

"This script calculates the relative area covered between the PTV and the Esophagus, as saves it to a CSV file"

outputdir = 'H:/NyttigData'

patient_db = get_current('PatientDB')

with open(outputdir+'/'+'AreaCoveredPTV.csv', 'a', newline='') as csvfile:

    csvwriter = csv.writer(csvfile)
    #csvwriter.writerow(['Pasient, Length of Esophagus, Slice thickness'])
    
    for i, p in enumerate(pat_numbers):
        pat_info = patient_db.QueryPatientInfo(Filter={'LastName': p})
        patient = patient_db.LoadPatient(PatientInfo=pat_info[0])
        case = patient.Cases[0]
        case.SetCurrent()
        plan = patient.Cases[0].TreatmentPlans[plans[i]]
        plan.SetCurrent()
        examination = patient.Cases[0].Examinations['CT 1']
        
        dose = plan.TreatmentCourse.TotalDose

        dg_esoph = dose.GetDoseGridRoi(RoiName='Esophagus')

        voxel_indices_esoph = dg_esoph.RoiVolumeDistribution.VoxelIndices
        
        dg_PTV = dose.GetDoseGridRoi(RoiName='PTV')

        voxel_indices_PTV = dg_PTV.RoiVolumeDistribution.VoxelIndices
        
        overlap_indices = 0
        
        for ind in voxel_indices_PTV:
            for i in voxel_indices_esoph:
                if (ind == i):
                    overlap_indices = 1
                    break
            if overlap_indices == 1:
                break
                    
        if overlap_indices == 1: 

            Volume_Matrix_ROI, slice_thick1 = GetVolumeDoseMatrix(plan, 'Esophagus')
        
        
            Volume_Matrix_ROI_intersection, slice_thick2 = GetVolumeDoseMatrix(plan, 'Esoph_intersec_PTV')
            
            RelativeAreaOverlap, LengthROI, LengthIntersec = relativeAreaCovered(Volume_Matrix_ROI, Volume_Matrix_ROI_intersection, slice_thick1)
        
        else: 
            RelativeAreaOverlap = []
            
            Volume_Matrix_ROI, slice_thick1 = GetVolumeDoseMatrix(plan, 'Esophagus')
            LengthROI = GetLengthOfROI(Volume_Matrix_ROI, slice_thick1)

        csvwriter.writerow([p] + [LengthROI, slice_thick1] + RelativeAreaOverlap)
