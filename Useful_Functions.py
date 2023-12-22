'''Useful functions'''

import numpy as np; import statistics; import csv

def createTotalROIDoseMatrix(Dose_matrix_1d, dim, indicies):
    Total_ROI_dose_matrix = np.zeros(dim[0]*dim[1]*dim[2])
    for index, dose in zip(indicies, Dose_matrix_1d):
        Total_ROI_dose_matrix[index] = dose
    
    return Total_ROI_dose_matrix.reshape(dim)

def createVolumeROIMatrix(Volume_matrix_1d, dim, indicies):
    Volume_ROI_matrix = np.zeros(dim[0]*dim[1]*dim[2])
    for index, vol in zip(indicies, Volume_matrix_1d):
        Volume_ROI_matrix[index] = vol
        
    return Volume_ROI_matrix.reshape(dim)
        
def collectData(filename):
    rows = []
    with open(filename, 'r') as file:
        csvreader = csv.reader(file)
        next(csvreader)
        for row in csvreader:
        
            y = [float(j) for j in row]
            rows.append(y)
    return rows

def ConvertDoseMatrixToEQD2(Dose_matrix_1d, num_frac, ab):
    doseMatrixEQD2 = np.zeros(len(Dose_matrix_1d))
    for i, dose in enumerate(Dose_matrix_1d):
        doseMatrixEQD2[i] = (dose/100)*(ab + (dose/(100*num_frac)))/(ab+2)
    
    return doseMatrixEQD2

def ConvertDoseMatrixToEQD2withTIME(Dose_matrix_1d, num_frac, ab, T, T_ref):
    doseMatrixEQD2 = np.zeros(len(Dose_matrix_1d))
    for i, dose in enumerate(Dose_matrix_1d):
         x = (dose/100)*(ab + (dose/(100*num_frac)))/(ab+2) - 0.8*(T-T_ref)
         if x>0:
             doseMatrixEQD2[i] = x
         else: 
             doseMatrixEQD2[i] = 0
    
    return doseMatrixEQD2

def getEUD(Dose_matrix_1D, Volume_matrix_1d):
    total_volume = 0 
    for vol in Volume_matrix_1d:
        total_volume += vol
    EUD = 0
    for dose, volume in zip(Dose_matrix_1D, Volume_matrix_1d):
        EUD += dose*(volume/total_volume)
    
    return EUD

def getVolume(Dose_matrix_1d, Vol_matrix_1d, N):
    doses = np.linspace(0, 60, N)
    volumes = np.zeros(N)
    for i, thr in enumerate(doses):
        for dose, vol in zip(Dose_matrix_1d, Vol_matrix_1d):
            if dose >= thr:
                volumes[i] += vol
    return volumes, doses

def getDifferentialVolume(Dose_matrix_1d, Vol_matrix_1d, bin_size, N):
    doses = np.linspace(0, 60, N)
    volumes = np.zeros(N)
    for i, thr in enumerate(doses):
        for dose, vol in zip(Dose_matrix_1d, Vol_matrix_1d):
            if dose < (thr + bin_size) and dose > (thr - bin_size):
                volumes[i] += vol
                
    return volumes, doses
    
def simpleAreaCoveredFAST(Vol_matrix_1d, Dose_matrix_1d, indicies, slice_thickness, dim, overlap_threshold, N):
    thresholds = np.linspace(0, 60, N)
    lengths = np.zeros(len(thresholds))
    for i, thr in enumerate(thresholds): 
        relevant_indices = []
        relevant_volumes = []
        
        for dose, ind in zip(Dose_matrix_1d, indicies):
            if dose > thr:
                relevant_indices.append(ind)
        
        for j in relevant_indices:
            index = np.where(indicies == j)[0][0]
            relevant_volumes.append(Vol_matrix_1d[index])
        
        area_ROI = np.zeros(dim[0])
        area_intersection = np.zeros(dim[0])
        relativeAreaCovered = []
        
        if len(relevant_indices) > 0:
        
            slice_num = indicies[0]//(dim[1]*dim[2])
            for k, elem in enumerate(indicies):
                if ((slice_num)*dim[1]*dim[2]-1) >= elem:
                    area_ROI[slice_num] += Vol_matrix_1d[k]
                else:
                    area_ROI[(elem//(dim[1]*dim[2]))+1] += Vol_matrix_1d[k]
                    slice_num += 1
        
            slice_num_1 = relevant_indices[0]//(dim[1]*dim[2])
            
            for z, elem_1 in enumerate(relevant_indices):
                if ((slice_num_1)*dim[1]*dim[2]-1) >= elem_1:
                    area_intersection[slice_num_1] += relevant_volumes[z]
                else:
                    area_intersection[(elem_1//(dim[1]*dim[2]))+1] += relevant_volumes[z]
                    slice_num_1 += 1
            
            for a_r, a_i in zip(area_ROI, area_intersection):
                if a_r != 0:
                    relativeAreaCovered.append(a_i/a_r)    
                    
        for elem in relativeAreaCovered:
            if elem > overlap_threshold:
                lengths[i] += slice_thickness
    
    return lengths, thresholds 

def getMaxPercentCircumferenceIrradiatedPerDoseLevel(Dose_matrix, Vol_matrix, dim, dose_level):
    Max_overlap = []
    
    for k in range(dim[0]):
        vol_roi = 0.0
        vol_intersection = 0.0
        for s in range(dim[1]):
            for r in range(dim[2]):
                vol_roi += (Vol_matrix[k, s, r])
                if Dose_matrix[k, s, r] > dose_level:
                    vol_intersection += (Vol_matrix[k, s, r])
        if vol_roi > 0:
            Max_overlap.append(vol_intersection/vol_roi)
    
    return max(Max_overlap)

def GetDoseStatisticsPerSliceNEW(Dose_matrix, Vol_matrix, indices, dim, thr, slice_thickness):
    
    Median_Dose_per_slice, EightyPercent_dose_per_slice, TwentyPercent_dose_per_slice  = [], [], []
    length_of_full_dose = 0
    
    for k in range(dim[0]):
        vol_in_slice = []
        doses_in_slice = []
        total_volume = 0
        for s in range(dim[1]):
            for r in range(dim[2]):
                if Dose_matrix[k, s, r] > 0:
                    doses_in_slice.append(Dose_matrix[k, s, r])
                    vol_in_slice.append(Vol_matrix[k, s, r])
                    
        if len(doses_in_slice) > 1:
            #print(len(doses_in_slice))
            
            y = np.argsort(doses_in_slice)
            
            rel_volumes_sorted = np.zeros(len(doses_in_slice))
            
            for i, elem in enumerate(y):
                rel_volumes_sorted[i] = vol_in_slice[elem]
                
            doses_in_slice.sort()
            
            total_volume = sum(rel_volumes_sorted) 
            
            sum_vol = 0
            
            for index, vol in enumerate(rel_volumes_sorted):
                sum_vol += vol
                if sum_vol >= 0.2*total_volume:
                    correct_index_80_prosent = index
                    break
                
            sum_vol_1 = 0
            
            for index_1, vol_1 in enumerate(rel_volumes_sorted):
                sum_vol_1 += vol_1
                if sum_vol_1 >= 0.8*total_volume:
                    if index_1 < len(rel_volumes_sorted)-1:
                        correct_index_20_prosent = index_1 + 1
                        break
                    else:
                        correct_index_20_prosent = index_1
                        break
                    

    
            
            
            Median_Dose_per_slice.append(statistics.median(doses_in_slice))
            
            SeventyfivePercent_dose_per_slice.append(doses_in_slice[correct_index_80_prosent])
            TwentyfivePercent_dose_per_slice.append(doses_in_slice[correct_index_20_prosent])
            
            
            if statistics.median(doses_in_slice) > thr:
                length_of_full_dose += slice_thickness
        
    return Median_Dose_per_slice, TwentyPercent_dose_per_slice, EightyPercent_dose_per_slice, length_of_full_dose

def createCircumferenceMatrix(vol_matrix_1d, dose_matrix, indices, dim):
    circumference_indices = []
    for vol, ind in zip(vol_matrix_1d, indices):
        if vol < 0.0187:
            circumference_indices.append(ind)
            
    circumference_matrix_1d = np.zeros(dim[0]*dim[1]*dim[2])
    
    for index in circumference_indices:
        circumference_matrix_1d[index] = dose_matrix.flatten()[index]
        
    circumference_matrix = circumference_matrix_1d.reshape(dim)
    
    return circumference_matrix, np.array(circumference_indices)

