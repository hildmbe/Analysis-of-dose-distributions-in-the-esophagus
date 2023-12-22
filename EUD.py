from connect import *
import csv

outputdir = 'H:/NyttigData'

patient_db = get_current('PatientDB')

n = 1
alpha_beta = 10

def compute_eud(dose, roi_name, n):
    
    dose_values = dose.DoseValues.DoseData.flatten()  
    dgr = dose.GetDoseGridRoi(RoiName = roi_name) 
    indices = dgr.RoiVolumeDistribution.VoxelIndices 
    relative_volumes = dgr.RoiVolumeDistribution.RelativeVolumes 

    dose_sum = 0.0
    scalefactor = 1.0 / max(dose_values) # Scale factor to prevent overflow with small n-values

    ##### EUD calculation
    for i, v in zip(indices, relative_volumes):
        d = dose_values[i] * scalefactor
        dose_sum += v * (d ** (1 / n))
    eud = (1 / scalefactor) * dose_sum ** n
    return eud


def compute_eqd2_eud(dose, num_frac, roi_name, ab, n):
    
    dose_values = dose.DoseValues.DoseData.flatten() 
    dgr = dose.GetDoseGridRoi(RoiName=roi_name)
    indices = dgr.RoiVolumeDistribution.VoxelIndices
    relative_volumes = dgr.RoiVolumeDistribution.RelativeVolumes 

    ##### EQD2 calculation
    dose_eqd2 = dose_values * (((dose_values/num_frac) + ab*100.0)/((2.0*100.0)+(ab*100.0)))

    ##### EUD calculation
    dose_sum = 0.0
    scalefactor = 1 / max(dose_eqd2)  # Scale factor to prevent overflow with small n-values
    for i, v in zip(indices, relative_volumes):
        d = dose_eqd2[i] * scalefactor 
        dose_sum += v * (d ** (1 / n))
    eud = (1 / scalefactor) * dose_sum ** n
    return eud
    

with open(outputdir+'/'+'EUD_intersec40.csv', 'w', newline='') as csvfile:

    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(['Pasient, EUD_phys, EUD_EQD2'])


    for i, p in enumerate(pat_numbers):
        pat_info = patient_db.QueryPatientInfo(Filter={'LastName': p})
        patient = patient_db.LoadPatient(PatientInfo=pat_info[0])
        case = patient.Cases[0]
        case.SetCurrent()
        plan = patient.Cases[0].TreatmentPlans[plans[i]]
        plan.SetCurrent()
        examination = patient.Cases[0].Examinations['CT 1']
        
        dose = plan.TreatmentCourse.TotalDose
        num_frac = plan.TreatmentCourse.TreatmentFractions[0].BeamSet.FractionationPattern.NumberOfFractions
        
        eud_phys = compute_eud(dose, 'Esophagus', n)
        eud_eqd2 = compute_eqd2_eud(dose, num_frac, 'Esophagus', alpha_beta, n)

        csvwriter.writerow([p] + [eud_phys/100] + [eud_eqd2/100])
        
   

    
    
