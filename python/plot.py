import pandas as pd
import matplotlib.pyplot as plt

example = 2

# Die Datei einlesen
dateipfad = '../output/cG-cG/example-' + str(example) + '/final_condition=0/convergence_table.txt'  # Pfad zu deiner Datei angeben
df = pd.read_csv(dateipfad, sep='|')  # Trennzeichen ist '|' und überspringt die Kopfzeile

# Entfernen von führenden und nachfolgenden Leerzeichen aus den Spaltennamen
df.columns = df.columns.str.strip()
# print(df.columns)
# Daten für dofs(time) und L2 extrahieren
dofs_time = df['dofs(time)']
L2_1 = df['L2']

# Die Datei einlesen
dateipfad = '../output/cG-cG/example-' + str(example) + '/final_condition=1/convergence_table.txt'  # Pfad zu deiner Datei angeben
df = pd.read_csv(dateipfad, sep='|')  # Trennzeichen ist '|' und überspringt die Kopfzeile

# Entfernen von führenden und nachfolgenden Leerzeichen aus den Spaltennamen
df.columns = df.columns.str.strip()
# print(df.columns)
# Daten für dofs(time) und L2 extrahieren
dofs_time = df['dofs(time)']
L2_2 = df['L2']

dateipfad = '../output/dG/example-' + str(example) + '/final_condition=0/convergence_table.txt'  # Pfad zu deiner Datei angeben
df = pd.read_csv(dateipfad, sep='|')  # Trennzeichen ist '|' und überspringt die Kopfzeile

# Entfernen von führenden und nachfolgenden Leerzeichen aus den Spaltennamen
df.columns = df.columns.str.strip()
# print(df.columns)
# Daten für dofs(time) und L2 extrahieren
dofs_time_3 = df['dofs(time)']
L2_3 = df['L2']

dateipfad = '../output/cG/example-' + str(example) + '/final_condition=0/convergence_table.txt'  # Pfad zu deiner Datei angeben
df = pd.read_csv(dateipfad, sep='|')  # Trennzeichen ist '|' und überspringt die Kopfzeile

# Entfernen von führenden und nachfolgenden Leerzeichen aus den Spaltennamen
df.columns = df.columns.str.strip()
# print(df.columns)
# Daten für dofs(time) und L2 extrahieren
dofs_time_4 = df['dofs(time)']
L2_4 = df['L2']

# Diagramm erstellen
plt.figure(figsize=(10, 6))
plt.plot(dofs_time, L2_1, marker='s', linestyle='-', color='b', label='cG(1)/cG(1)-cG(1) (with normal conditions)')  # Kurve zeichnen
plt.plot(dofs_time, L2_2, marker='x', linestyle='-', color='r', label='cG(1)/cG(1)-cG(1) (with final condition)')  # Kurve zeichnen
plt.plot(dofs_time_3, L2_3, marker='o', linestyle='-', color='k', label='cG(1)/dG(0)')  # Kurve zeichnen
plt.plot(dofs_time_4, L2_4, marker='|', linestyle='-', color='g', label='cG(1)/cG(1)')  # Kurve zeichnen
plt.yscale('log')
plt.xscale('log')
# plt.grid(True)
plt.legend(loc='upper right')
plt.show()
