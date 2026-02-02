import barnaba as bb
from collections import Counter
import os

with open("pares_bases.csv", 'w') as f:
	f.write("%10s %10s %4s" % ("Residuo 1","Residuo 2", "Tipo"))
	f.write('\n')

with open("stacking.csv", 'w') as f2:
	f2.write("%10s %10s %4s" % ("Residuo 1","Residuo 2", "Tipo"))
	f2.write('\n')

with open("dot_plot.csv", 'w') as f3:
	f3.write(" ")

# annotate
for aptamero in range(1, 101, 1):
	pdb = str("md_complejo_"+str(aptamero)+".pdb")

	try:
		stackings, pairings, res = bb.annotate(pdb)

		dotbr, seq = bb.dot_bracket(pairings,res)

		with open("dot_plot.csv", 'a') as f3:
			for j in range(len(dotbr)):
				f3.write(dotbr[j])
				f3.write("\n")

		filename = str(str(aptamero)+"_dot_plot.csv")

		with open(filename,'w') as f4:
			for j in range(len(dotbr)):
				f4.write(dotbr[j])
				f4.write("\n")

# list base pairings
		for p in range(len(pairings[0][0])):
			res1 = res[pairings[0][0][p][0]]
			res2 = res[pairings[0][0][p][1]]
			interaction =  pairings[0][1][p]

			with open("pares_bases.csv", 'a') as f:
				f.write("%10s %10s %4s" % (res1,res2,interaction))
				f.write('\n')

# list base-stackings
		for p in range(len(stackings[0][0])):
			res1 = res[stackings[0][0][p][0]]
			res2 = res[stackings[0][0][p][1]]
			interaction =  stackings[0][1][p]

			with open("stacking.csv",'a') as f2:
				f2.write("%10s %10s %4s" % (res1,res2,interaction))
				f2.write('\n')

	except:

		print("Qué modelo pedorro el "+pdb)

os.system(f"grep -v '^[.]*$' dot_plot.csv > dot_plot_filtrado.csv")

def consenso_con_paréntesis(lines, umbral=0.4):
    columnas = zip(*lines)
    consenso = []
    total_lineas = len(lines)

    for col in columnas:
        c = Counter(col)
        parentesis_total = c["("] + c[")"]

        if parentesis_total / total_lineas >= umbral:
            # Hay suficientes paréntesis, elegimos el más frecuente entre ellos
            if c["("] >= c[")"]:
                consenso.append("(")
            else:
                consenso.append(")")
        else:
            # No hay suficientes paréntesis, usamos el carácter más común general
            consenso.append(c.most_common(1)[0][0])

    return ''.join(consenso)


def linea_mas_comun(lines):
    c = Counter(lines)
    mas_comun, cantidad = c.most_common(1)[0]
    return mas_comun, cantidad

archivo = "dot_plot_filtrado.csv"

with open(archivo) as f:
    lines = [line.strip() for line in f if len(line.strip()) == 40]

consenso = consenso_con_paréntesis(lines, umbral=0.20)
print(consenso)

linea, repeticiones = linea_mas_comun(lines)

print(f"Línea más común (repetida {repeticiones} veces):\n{linea}")

