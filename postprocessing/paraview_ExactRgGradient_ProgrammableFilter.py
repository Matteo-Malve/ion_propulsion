
# GradRg x,y components
input0 = inputs[0]
radiusDataArray = input0.PointData["radius"]
xDataArray = input0.PointData["x"]
yDataArray = input0.PointData["y"]
numPoints = input0.GetNumberOfPoints()
RgGrad_x_Array=numpy.zeros(numPoints)
RgGrad_y_Array=numpy.zeros(numPoints)
dRgdr_Array=numpy.zeros(numPoints)
R = 7.5e-4
Ve = 20000.0
for i in range(numPoints):
  if radiusDataArray[i] < 3*R:
    RgGrad_x_Array[i] = xDataArray[i]/radiusDataArray[i]
    dRgdr_Array[i] = Ve/R * ((radiusDataArray[i]-R)/(2*R)-1)
    RgGrad_y_Array[i] = Ve/R * ((radiusDataArray[i]-R)/(2*R)-1) * yDataArray[i]/radiusDataArray[i]
  else:
    RgGrad_x_Array[i] = 0.0
    RgGrad_y_Array[i] = 0.0
    dRgdr_Array[i] = 0.0
output.PointData.append(RgGrad_x_Array, "RgGrad_x")
output.PointData.append(RgGrad_y_Array, "RgGrad_y")
output.PointData.append(dRgdr_Array, "dRgdr")

# Then put components together in a vector with Calculator
