import os


import os
#Path for directories
pathresult4 = './results7/'
pathresult8 = './results8/'

#Initialize the list of files
lstresult4files  = []
lstresult8files = []


lstresult4  = []
lstresult8 = []

images = ["im01.ppm","im02.ppm","im03.ppm"]
 
#List all files in the directory pathImages
lstDir = os.walk(pathresult4)
for root, dirs, files in lstDir:
    for fichero in files:
        (nombreFichero, extension) = os.path.splitext(fichero)
	lstresult4files.append(nombreFichero+extension)

#List all the file in the directory pathkernel
lstDir = os.walk(pathresult8)
for root, dirs, files in lstDir:
    for fichero in files:
        (nombreFichero, extension) = os.path.splitext(fichero)
	lstresult8files.append(nombreFichero+extension)

for result4 in lstresult4files:
	if result4 != ".DS_Store":
		fresult4 = open(pathresult4+result4, "r")
		i=0
		dataSize = 0
		dataSizeFinish = False
		cores = 0
		for line in fresult4:
			if i==1:
				cores=line
			if i==3:
				lstresult4.append(result4.split("_")[0]+";"+line.replace(".",",").rstrip("\n")+"\n")
			i+=1
		fresult4.close()

for result8 in lstresult8files:
	if result8 != ".DS_Store":
		fresult8 = open(pathresult8+result8, "r")
		i=0
		cores = 0
		for line in fresult8:
			if i==1:
				cores=line
			if i==3:
				lstresult8.append(result8.split("_")[0]+";"+line.replace(".",",").rstrip("\n")+"\n")
			i+=1
		fresult8.close()
 
out_results = open("global_results.csv", "w")  

for result4 in lstresult4:
	out_results.write(result4)
out_results.write("\n")

for result8 in lstresult8:
	out_results.write(result8)

out_results.close()
