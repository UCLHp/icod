import os
import numpy as np
from vedo import *
from pydicom import *

# DICOM directory
dicomRoot = "O:\protons\Work in Progress\AlexG\icod\dcm"
# RT struct filename
structFilename = "RS.6skzOZqM4Icdg0Sv20WBrgG2F.OD CNS.dcm"
# RT struct name
structName = "BODY"

# RT struct import
dicomRoot = os.path.normpath(dicomRoot)
structFilename = os.path.join(dicomRoot,structFilename)
structFile = dcmread(structFilename)

# RT struct data
structSeq_name = structFile.StructureSetROISequence
structSeq_data = structFile.ROIContourSequence
numStruct = len(structSeq_name)
for s1, s2 in zip(structSeq_name,structSeq_data):
    if s1.ROIName == structName:
        roi_number = int(s1.ROINumber)
        contour_sequence = s2.ContourSequence
        contour_data = []
        for contour_slice in contour_sequence:
            contour_data.extend(contour_slice.ContourData)

X_data = np.array(contour_data[0::3])
Y_data = np.array(contour_data[1::3])
Z_data = np.array(contour_data[2::3])
struct_data = [X_data,Y_data,Z_data]
# load relevant Contour Data (3006,050)
pts1 = Points(struct_data)
# clean opints
pts1.clean(tol=0.005)
##surface from points
bounds=[X_data.min(), X_data.max(), Y_data.min(), Y_data.max(), Z_data.min(), Z_data.max()]
reco = recoSurface(pts1, dims=150, bounds=bounds, radius=10)
reco_cap = reco.cap(returnCap=True)
body_mesh = mesh.merge(reco, reco_cap)
plt = Plotter()
plt.show(body_mesh, axes=1, zoom=1.2, interactive=1).close()