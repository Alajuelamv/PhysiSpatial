import pandas as pd

tabla = pd.read_csv("set_up_physiboss.csv")

tabla["col"] = tabla["col"] -min(tabla["col"])
tabla["row"] = tabla["row"] -min(tabla["row"])
tabla["row"] = (max(tabla["row"])+min(tabla["row"])-tabla["row"])
tabla["col"] = tabla["col"] +5
tabla["row"] = tabla["row"] +5

endo = tabla[tabla.kmeans == "Endothelial cells"]
no_endo = tabla[tabla.kmeans != "Endothelial cells"]

no_endo["col"] = no_endo["col"] *20 +10
no_endo["row"] = no_endo["row"] *20 +10

rango_col = (max(tabla["col"]) - min(tabla["col"]) + 10)*10
rango_row = (max(tabla["row"]) - min(tabla["row"]) +10 )*10

no_endo["col"] = no_endo["col"] - rango_col
no_endo["row"] = no_endo["row"] - rango_row

endo.drop(["kmeans"],axis=1,inplace=True)
endo.to_csv("setup_endo_norm.csv",index =False,header=None)

no_endo["kmeans"] = no_endo["kmeans"].astype(int) -1
no_endo.to_csv("setup_cells_norm.csv",index =False,header=None)


# here generating the setting up of the ecm

ecm = pd.DataFrame()

from itertools import product

def foo(rangecol,rangerow):
    for col, rows in product(rangecol,rangerow):
        yield [col,rows]

rcol = range(min(tabla["col"])-5,max(tabla["col"])+5)
rrow = range(min(tabla["row"])-5,max(tabla["row"])+5)

data = []
for row in foo(rcol,rrow):
    data.append(row)

df = pd.DataFrame(data)
df.columns = ['col','row']

tabla_coladnrow = pd.DataFrame()
tabla_coladnrow["col"] = tabla["col"]
tabla_coladnrow["row"] = tabla["row"]

ecm = pd.concat([df,tabla_coladnrow]).drop_duplicates(keep=False)
ecm["z"] = 0
ecm.reset_index(inplace=True)


false_ecm = list()
for i in range(1,len(ecm["index"]) -2):
    x = ecm.loc[i+1]["row"]
    y = ecm.loc[i-1]["row"]
    z = ecm.loc[i]["row"]
    if x-z != 1 and z-y != 1:
        false_ecm.append(i)
        
for j in false_ecm:
    ecm.drop(j,inplace=True)
   
   
ecm.drop("index",axis= 1,inplace=True)     
ecm.to_csv("setup_ecm.csv",index =False,header=None)
        
    


