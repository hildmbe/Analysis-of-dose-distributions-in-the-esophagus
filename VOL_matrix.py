''' Funksjonen GetVolumeDoseMatrix tar inn tar inn pasientinfo og ROI av interesse. Funksjonen returner en volum-matrise
der alle elementer som ligger utenfor valgt ROI er satt til 0, mens alle elementer som ligger innenfor valgt ROI er satt
til volumet av ROI-en i en gitt dose-grid voxel. Volum-matrisen har samme dimensjoner som den totale dosematrisen.

Ved å lage en volum-matrise for Øsofagus og en volum-matrise for Overlapp Øsofagus/40Gy, så kan den relative andelen
av >40Gy som dekker Øsofagus beregnes vha. funksjonen relativeAreaCovered'''

'''Volum-matriser er lagd vha funskjonen ROI_dg.RoiVolumeDistribution.RelativeVolumes. 
Definisjon på RelativeVolumes: the fraction of the total volume of the ROI in a given dose grid voxel.
relative_volumes blir altså en array med lik lengde som voxel_indices. Summen av alle elementene i relative_volumes 
blir 1. Må dermed gange med totalt volum av ROI for a finne volum av ROI in hver voxel. '''

from connect import *


def GetVolumeDoseMatrix(plan, ROI):
    dose = plan.TreatmentCourse.TotalDose

    dose_matrix = dose.DoseValues.DoseData
    dose_matrix_1d = dose_matrix.flatten()

    dim = dose_matrix.shape

    dose.UpdateDoseGridStructures()

    ROI_dg = dose.GetDoseGridRoi(RoiName=ROI)

    voxel_indices= ROI_dg.RoiVolumeDistribution.VoxelIndices
    relative_volumes= ROI_dg.RoiVolumeDistribution.RelativeVolumes
    total_vol = ROI_dg.RoiVolumeDistribution.TotalVolume
    slice_thickness = ROI_dg.InDoseGrid.VoxelSize['z']

    volumeMatrix_1d = np.zeros(len(dose_matrix_1d))

    for index, elem in enumerate(voxel_indices):
        volumeMatrix_1d[elem] = relative_volumes[index] * total_vol

    volumeMatrix = volumeMatrix_1d.reshape(dim)

    return volumeMatrix, slice_thickness

def relativeAreaCovered(matrix_ROI, matrix_intersection, slice_thickness):
    dim = matrix_ROI.shape

    relativeAreaCovered = []
    length_of_ROI = 0
    length_of_intersection = 0

    for k in range(dim[0]):
        area_ROI = 0
        area_intersection = 0
        for s in range(dim[1]):
            for r in range(dim[2]):
                area_ROI += matrix_ROI[k, s, r]
                area_intersection += matrix_intersection[k, s, r]
        if (area_ROI != 0):
            relativeAreaCovered.append(area_intersection/area_ROI)
            length_of_ROI += slice_thickness
            if (area_intersection / area_ROI > 0.90):
                length_of_intersection += slice_thickness

    return relativeAreaCovered, length_of_ROI, length_of_intersection

def GetLengthOfROI(matrix_ROI, slice_thickness):
    dim = matrix_ROI.shape
    length_of_ROI = 0 
    
    for k in range(dim[0]):
        area_ROI = 0
        for s in range(dim[1]):
            for r in range(dim[2]):
                area_ROI += matrix_ROI[k, s, r]
        if (area_ROI != 0):
            length_of_ROI += slice_thickness
     
    return length_of_ROI
    
    
    





