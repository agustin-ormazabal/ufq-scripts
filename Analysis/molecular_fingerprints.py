import MDAnalysis as mda
from prolif import Molecule, Fingerprint
from MDAnalysis.analysis.align import alignto
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
import pandas as pd
import numpy as np

with open("agus_score.out", 'w') as file_out:
	file_out.write("Modelo"+"\t"+"Agus_score"+"\t"+"Agus_score2"+"\n")

interacciones_frecuentes = {}

for model in range(1, 101, 1):

	try:
		print("Model", model)

		lista_names_prot_rna = []
		lista_names_prot = []
		lista_tipo_interacciones = []
		referencia_colores = {}
		lista_colores = []

		model = str(model)

# Cargar la trayectoria
		u = mda.Universe("complejo_"+model+".prmtop", "md_complejo_"+model+".nc")

# Seleccionar proteína y ARN
		protein = u.select_atoms("protein")
		rna = u.select_atoms("resname A DT G C")

# Convertir cada residuo en una molécula de ProLIF
# Asegúrate de convertir correctamente cada residuo individualmente en un AtomGroup

		prot_mols = [Molecule.from_mda(protein.select_atoms(f"resid {res.resid}")) for res in protein.residues]
		rna_mols = [Molecule.from_mda(rna.select_atoms(f"resid {res.resid}")) for res in rna.residues]

# Alinear la trayectoria (opcional, si es necesario)
		alignto(u, select="protein", reference=protein)

# Calcular las interacciones
		fp = Fingerprint()
		interactions = fp.run(u.trajectory,protein,rna)

# Guardar el resultado
		df = interactions.to_dataframe()
		df.to_csv("interacciones_"+model+".csv")

# Cargar los datos de interacciones desde el CSV
		interaction_df = pd.read_csv("interacciones_"+model+".csv")
		nframes = int(len(interaction_df))-int(3)
		puntajes = int(nframes*10 + 10)

# Cargar los datos desde el CSV sin encabezados
		df = pd.read_csv("interacciones_"+model+".csv", header=None)

# Guarda los nombres de los residuos que interactúan
		n_interacciones = int(len(df.columns)) 
	
		for interaccion in range(1, n_interacciones, 1):

			names_prot = df.iloc[0,interaccion]
			names_rna = df.iloc[1,interaccion]
			names_nc  = names_rna[0]

			lista_names_prot.append(str(names_prot+"-"+names_nc))

			names_prot_rna = str(names_prot + "-" + names_rna)
			lista_names_prot_rna.append(names_prot_rna)

			tipo_interacciones = df.iloc[2,interaccion]
			lista_tipo_interacciones.append(tipo_interacciones)
			interacciones_unicas = list(set(lista_tipo_interacciones))
			
# Define los colores en función del tipo de interacción:

		for interaccion in interacciones_unicas:

			if interaccion == 'HBDonor': 
				referencia_colores[interaccion] = int(5)
			elif interaccion == 'HBAcceptor':
                        	referencia_colores[interaccion] = int(5)
			elif interaccion == 'PiStacking':
				referencia_colores[interaccion] = int(5)
			elif interaccion == 'CationPi':
				referencia_colores[interaccion] = int(10)
			elif interaccion == 'Cationic':
                        	referencia_colores[interaccion] = int(10)
			elif interaccion == 'Anionic':
				referencia_colores[interaccion] = int(10)
			elif interaccion == 'Hydrophobic':
				referencia_colores[interaccion] = int(1)
			elif interaccion == 'VdWContact':
				referencia_colores[interaccion] = int(1)

		for interaccion in lista_tipo_interacciones:
			lista_colores.append(referencia_colores[interaccion])

# Eliminar las primeras 3 filas y la primera columna
		df = df.iloc[3:, 1:]

# Convertir los valores a booleanos
		df = df.map(lambda x: True if str(x).lower() != 'false' else False)

# Convertir los valores booleanos a enteros (0 y 1)
		df = df.astype(int)

# Transponer el DataFrame
		df = df.T

		color_matrix = df.copy()

# Reemplazar valores en color_matrix
#  Si es True, se asigna el índice del color basado en lista_tipo_interacciones.
#  Si es False, se asigna -1 para que se muestre como blanco.

		positivos = 0
		positivos_totales = 0
		negativos = 0

		for i in range(len(df.index)):  # Recorrer filas
			for j in range(len(df.columns)):  # Recorrer columnas
				if df.iloc[i, j] == 1:
					color_matrix.iloc[i, j] = lista_colores[i]  # Asignar el color RGB
					if lista_colores[i] >= 5:
						positivos = positivos + 1
					if lista_colores[i] > 0:
						positivos_totales = positivos_totales + 1

				else:  
					color_matrix.iloc[i, j] = int(0)  # Asignar negro (o (1, 1, 1) para blanco)

		celdas_totales = int(len(df.index)) + int(len(df.columns))

		proporcion = int(positivos) / int(celdas_totales)
		proporcion_total = int(positivos_totales) / int(celdas_totales)

		print(positivos, positivos_totales,"LA PROPORCION ES",proporcion_total)

########################################################################################################################

		for i in range(len(df.index)):
			suma_fila = color_matrix.iloc[i, :].sum()
			proporcion = suma_fila / puntajes
			nombre_interaccion = lista_names_prot[i]

			if nombre_interaccion not in interacciones_frecuentes:			
				interacciones_frecuentes[nombre_interaccion] = proporcion
			else:
				interacciones_frecuentes[nombre_interaccion] = interacciones_frecuentes[nombre_interaccion] + proporcion
	
########################################################################################################################


# Graficar el heatmap con el colormap personalizado

		colors = ["#ffffff", "#ffcc66", "#ccff66", "#33cc33"]  # Blanco, anaranjado amarillento, verde amarillento, verde bonito
		bounds = [0, 1, 5, 10, 11]  # Límites para cada color

		cmap = mcolors.ListedColormap(colors)
		norm = mcolors.BoundaryNorm(bounds, cmap.N)

		plt.figure(figsize=(12, 16))
		sns.heatmap(color_matrix.astype(float),  # Asegurar que los valores sean numéricos
            	linewidths=.5,
            	linecolor='lightgrey',
	    	cmap=cmap,
	    	norm=norm,
            	cbar=False,  # Mostrar barra de colores para ver qué color es cada interacción
            	yticklabels=lista_names_prot_rna)

		plt.title('Interacciones Proteína-RNA: Modelo '+model, fontsize=16)
		plt.xlabel('#Frame', fontsize=12)  # Intercambiar las etiquetas de los ejes
		plt.ylabel('Interacción', fontsize=12)
		plt.xticks(rotation=45, ha='right')
		plt.yticks(rotation=0)

		plt.tight_layout()

# Guardar el heatmap como un archivo PNG
		plt.savefig('heatmap_'+model+'.png')
	
	except:

		print("El modelo",model,"está rancio")
		proporcion = 0
		proporcion_total = 0

	with open("agus_score.out", 'a') as file_out:
		file_out.write(model+"\t"+str(proporcion)+"\t"+str(proporcion_total)+"\n")
	

#  Crea un archivo de interacciones frecuentes en donde la primera columna
# muestra a la interacción en cuestión, y la segunda, el porcentaje de veces
# que fue observada sobre el total de snapshots de la trayectoria. 

with open("interacciones_frecuentes.txt", "w") as f:
    f.write("Interacción\tPrevalencia\n")  
    for key, value in interacciones_frecuentes.items():
        f.write(f"{key}\t{value}\n")
