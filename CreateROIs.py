from connect import *

def createRoiFromDose(plan, case, threshold_level):
    ROI_name = str(threshold_level)+'Gy'
    ROI = case.PatientModel.CreateRoi(Name=ROI_name, Color='Yellow')
    ROI.CreateRoiGeometryFromDose(DoseDistribution= plan.TreatmentCourse.TotalDose, ThresholdLevel=threshold_level*100)

def createIntersecROI(ROI_name_1, ROI_name_2, case, examination):
    ROI_intersec = case.PatientModel.CreateRoi(Name= ROI_name_1 +'_intersec_' + ROI_name_2, Color="Red",
                                                     TissueName=None, RbeCellTypeName=None, RoiMaterial=None)

    ROI_intersec.SetAlgebraExpression(ExpressionA={'Operation': "Union", 'SourceRoiNames': [ROI_name_1],
                                                         'MarginSettings': {'Type': "Expand", 'Superior': 0,
                                                                            'Inferior': 0, 'Anterior': 0,
                                                                            'Posterior': 0, 'Right': 0, 'Left': 0}},
                                            ExpressionB={'Operation': "Union", 'SourceRoiNames': [ROI_name_2],
                                                         'MarginSettings': {'Type': "Expand", 'Superior': 0,
                                                                            'Inferior': 0, 'Anterior': 0,
                                                                            'Posterior': 0, 'Right': 0, 'Left': 0}},
                                            ResultOperation='Intersection',
                                            ResultMarginSettings={'Type': "Expand", 'Superior': 0, 'Inferior': 0,
                                                                  'Anterior': 0, 'Posterior': 0, 'Right': 0, 'Left': 0})
    ROI_intersec.UpdateDerivedGeometry(Examination=examination, Algorithm='Auto')
    
    
    