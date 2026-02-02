import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
nnucleotidos= 40

interacciones_frecuentes = "interacciones_frecuentes_aptamero.txt"
interacciones_frecuentes_nucleotido = "interacciones_aptameros_nucleotido.txt"

# Leer archivo con las interacciones por posicion
interaccion_posicion = pd.read_csv(interacciones_frecuentes, delimiter='\t', header=None)
interaccion_posicion.columns = ['Posicion', 'Frecuencia']
interaccion_posicion = interaccion_posicion.sort_values(by='Posicion')

# Leer archivo con las interacciones por nucleotido
interaccion_nucleotido = pd.read_csv(interacciones_frecuentes_nucleotido, delimiter='\t', header=None)
interaccion_nucleotido.columns = ['Nucleotido', 'Frecuencia']


lista_A = []
lista_G = []
lista_T = []
lista_C = []

for posicion in interaccion_posicion['Posicion']:
	

	check_A = len(lista_A)
	check_G = len(lista_G)
	check_T = len(lista_T)
	check_C = len(lista_C)

	for nucleotido in interaccion_nucleotido['Nucleotido']:

		if int(posicion) < 10 and len(nucleotido) < 3 and str(posicion) in nucleotido:
			frecuencia = int(interaccion_nucleotido.loc[interaccion_nucleotido['Nucleotido'] == nucleotido, 'Frecuencia'].values[0])
			
			if nucleotido[0] == "A":
				lista_A.append(frecuencia)
			elif nucleotido[0] == "G":
				lista_G.append(frecuencia)
			elif nucleotido[0] == "T":
				lista_T.append(frecuencia)
			elif nucleotido[0] == "C":
				lista_C.append(frecuencia)	

		elif int(posicion) >= 10 and len(nucleotido) == 3 and str(posicion) in nucleotido:
			frecuencia = int(interaccion_nucleotido.loc[interaccion_nucleotido['Nucleotido'] == nucleotido, 'Frecuencia'].values[0])

			if nucleotido[0] == "A":
				lista_A.append(frecuencia)
			elif nucleotido[0] == "G":
				lista_G.append(frecuencia)
			elif nucleotido[0] == "T":
				lista_T.append(frecuencia)
			elif nucleotido[0] == "C":
				lista_C.append(frecuencia)
		
	if len(lista_A) == check_A:
		lista_A.append(int(0))
		check_A = len(lista_A)
	elif len(lista_G) == check_G:
		lista_G.append(int(0))
		check_G = len(lista_G)
	elif len(lista_T) == check_T:
		lista_T.append(int(0))
		check_T = len(lista_T)
	elif len(lista_C) == check_C:
		lista_C.append(int(0))
		check_C = len(lista_C) 

def rellenar(lista):
    return lista + [0] * (nnucleotidos - len(lista))

lista_A = rellenar(lista_A)
lista_G = rellenar(lista_G)
lista_T = rellenar(lista_T)
lista_C = rellenar(lista_C)

# Posiciones
x = list(range(1, nnucleotidos + 1))

# Setup del gráfico
fig, ax = plt.subplots(figsize=(12, 6))

# Graficar burbujas
ax.scatter(x, [3]*(nnucleotidos), s=[v*v*v for v in lista_A], label='A', alpha=0.5, color='blue', edgecolors='none')
ax.scatter(x, [2]*(nnucleotidos), s=[v*v*v for v in lista_G], label='G', alpha=0.5, color='green', edgecolors='none')
ax.scatter(x, [1]*(nnucleotidos), s=[v*v*v for v in lista_T], label='T', alpha=0.5, color='red', edgecolors='none')
ax.scatter(x, [0]*(nnucleotidos), s=[v*v*v for v in lista_C], label='C', alpha=0.5, color='orange', edgecolors='none')

# Etiquetas y estilo
ax.set_yticks([0,1,2,3])
ax.set_yticklabels(['C', 'T', 'G', 'A'])
ax.set_xlabel("Posición")
ax.set_title("Gráfico de burbujas por nucleótido y posición")
ax.grid(True, axis='x', linestyle='--', alpha=0.4)
ax.legend(loc='upper right')

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

plt.tight_layout()
plt.savefig("libreria.png", dpi=300, bbox_inches='tight')
