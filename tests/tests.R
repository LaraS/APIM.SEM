

res <- apim_sem(dname = "gender", xname = "neurotocism", yname="distress", 
            namex1 = "mneuro", namex2 = "fneuro", namey1 = "mdistr", 
            namey2 = "fdistr",labs1 = "husband", labs2 = "wive", 
            lab1 = "husbands", lab2 = "wives",data = campbell)

res
printTables(res)

# citation()
# figures(temploc= "H:/home/Doctoraat/Studie 2 - APIM/templates/")


