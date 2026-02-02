import numpy as np
from tqdm import tqdm
import MDAnalysis as mda
import matplotlib.pyplot as plt
from MDAnalysis.analysis.rms import rmsd
from MDAnalysis.analysis.align import AlignTraj

topology   = "complejo_1.prmtop" 
trajectory = "md_complejo_1.nc"
u = mda.Universe(topology, trajectory)

###################################
# Seleccionar los residuos de ADN #
###################################

dna = u.select_atoms('resname A G DT and name P')	##############################
n_frames = len(u.trajectory) 				# Guarda cuántos frames hay. #
n_atoms = len(dna)					# Guarda cuántos átomos hay. #
							##############################
# Guardar coordenadas de cada frame
coords_all = np.empty((n_frames, n_atoms, 3))

for i, ts in enumerate(u.trajectory):
    coords_all[i] = dna.positions.copy()

###########################
# Calcular matriz de RMSD #
###########################

#  El objetivo es determinar cuál es el frame que, en promedio,
# tiene el menor RMSD respecto a todos los demás.
								###################################################################################
rmsd_matrix = np.zeros((n_frames, n_frames))			# Arma una matriz de ceros, con dimensiones que dependen del número de snapshots. #
for i in tqdm(range(n_frames), desc="Calculando RMSD"):		# TQDM es una funcion para que te grafique barritas mientras hace el cálculo.     #
								###################################################################################

# Para cada frame, calcula su RMSD contra todos los demás

    for j in range(i, n_frames):						
        rms = rmsd(coords_all[i], coords_all[j], center=True, superposition=True)
        rmsd_matrix[i, j] = rms		#################################
        rmsd_matrix[j, i] = rms		# Porque la matriz es simétrica #
					#################################
   

plt.figure(figsize=(8, 6))
im = plt.imshow(rmsd_matrix, cmap='viridis', interpolation='nearest')
plt.title('RMSD Matrix Heatmap')
plt.xlabel('Frame index')
plt.ylabel('Frame index')
plt.colorbar(im, label='RMSD (Å)')
plt.tight_layout()
plt.savefig('rmsd_matrix.png', dpi=300)


# Calcular promedio por fila
mean_rmsd = rmsd_matrix.mean(axis=1)	################################################
ref_idx = np.argmin(mean_rmsd)		# Elige el snapshot con el menor RMSD promedio #
					################################################

########################
# ARRANCA EL PCA FRUTA #
########################

# 1. Alinear TODA la trayectoria al frame de referencia #

for i in tqdm(range(n_frames), desc="Alineando"):
    
    if i == ref_idx: continue 		# Saltar referencia (ya está alineada consigo misma)
        
    # Calcular matriz de rotación óptima
    R, rmsd_val = mda.analysis.align.rotation_matrix(coords_all[i], coords_all[ref_idx])

    # Aplicar rotación a las coordenadas
    coords_all[i] = np.dot(coords_all[i], R.T)

# 2. Centrar todas las estructuras en el origen

for i in tqdm(range(n_frames), desc="Centrando"):	#########################################################################################################
    centroid = np.mean(coords_all[i], axis=0)		# Acá tengo dudas, porque fija el origen con el promedio. ¿Debería ser con la estructura de referencia? #
    coords_all[i] -= centroid				#########################################################################################################

# 3. Calcular matriz de covarianza
							#################################################################	
n_coords = 3 * n_atoms					# Para fijar el número de coordenadas.				#
cov = np.zeros((n_coords, n_coords))			# Arma una matriz llena de ceros.				#
ref_structure = coords_all[ref_idx].reshape(-1)		# Fija la estructura de referencia para calcular la sumatoria.	#
							#################################################################

								#########################################
for i in tqdm(range(n_frames), desc="Calculando covarianza"):	#  Arranca un loop			#
    coords_flat = coords_all[i].reshape(-1)			# para que a cada frame.		#
    delta = coords_flat - ref_structure				# le reste la estructura de referencia.	#
    cov += np.outer(delta, delta)				# y calcule la covarianza.		#
								#########################################

		 ##################################################
cov /= n_frames  # Normalizar en función de la cantidad de frames #
		 ##################################################

# 4. Diagonaliza la matriz
eigenvalues, eigenvectors = np.linalg.eigh(cov)

					######################
idx = np.argsort(eigenvalues)[::-1]	#  Con esto retiene  #
eigenvalues = eigenvalues[idx]		# los autovalores    #
eigenvectors = eigenvectors[:, idx]	# y los autovectores #
					######################

# 5. Calcular el porcentaje de cada autovalores sobre el total	############################################
total_variance = np.sum(eigenvalues)				# Autovalores totales.			   #
explained_variance = 100 * eigenvalues / total_variance		# Porcentaje explicado por cada autovalor. #
								############################################

# 6. Proyección 						#############################################################################
k_components = 2						#  OJO, acá definís sobre qué subespacio proyectás.			    #
top_k_eigenvectors = eigenvectors[:, :k_components]		#  Esta función elige los autovectores con mayor varianza,		    #
ref_coords_flat = ref_structure					# porque así como está programado, no te quedan ordenados de menor a mayor  #
projections = np.zeros((n_frames, k_components))		#  Con esto armás una matriz de ceros para la proyección.		    #
								#############################################################################



for i in range(n_frames):
    current_coords = coords_all[i].reshape(-1)			##########################################
    deviation = current_coords - ref_coords_flat		# Este es el bloque con el que proyectás #
    projections[i] = np.dot(deviation, top_k_eigenvectors)	##########################################

# 7. Graficar la proyección 

print(projections[:, 0])
print(projections[:, 1])

plt.figure(figsize=(10, 8))

plt.plot(projections[:, 0], projections[:, 1],
         marker='o',        # Puntos
         linestyle='-',     # Línea
         color='blue',
         alpha=0.7,
         markersize=6)

plt.axhline(0, color='k', linestyle='--', alpha=0.3)
plt.axvline(0, color='k', linestyle='--', alpha=0.3)
plt.xlabel(f'PC1 ({explained_variance[0]:.2f}%)')
plt.ylabel(f'PC2 ({explained_variance[1]:.2f}%)')
plt.title('PCA fruta')
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('pca_projection.png', dpi=300)

# 8. Graficar el decaimiento de autovalores.

n_show = 20
explained_cumsum = np.cumsum(explained_variance[:n_show])

plt.figure(figsize=(10, 6))
bars = plt.bar(range(1, n_show + 1), explained_variance[:n_show],
               color='cornflowerblue', width=0.7)

for i, (bar, cum_val) in enumerate(zip(bars, explained_cumsum)):
    plt.text(bar.get_x() + bar.get_width()/2, 
             bar.get_height() + 0.5, 
             f"{cum_val:.1f}%", 
             ha='center', va='bottom', fontsize=9)

plt.xlabel('Autovalor', fontsize=12)
plt.ylabel('Porcentaje (%)', fontsize=12)
plt.title('Decaimiento de los primeros 20 autovalores + Porcentajes acumulados', fontsize=14)
plt.xticks(range(1, n_show + 1))
plt.grid(axis='y', linestyle='--', alpha=0.4)
plt.tight_layout()
plt.savefig('eigenvalue_decay.png', dpi=300)
