import os
import re
import csv
from Bio import Entrez, SeqIO

Entrez.email = "bvasquezp@utem.cl"

def obtener_secuencia_ncbi(accession_id, archivo_salida):
    """Descarga o lee la secuencia en formato FASTA desde el NCBI."""
    if not os.path.exists(archivo_salida):
        print(f"Iniciando solicitud al NCBI para el acceso: {accession_id}...")
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
            with open(archivo_salida, "w") as out_handle:
                out_handle.write(handle.read())
            handle.close()
        except Exception as e:
            print(f"Error de conexion o descarga: {e}")
            return None
    
    registro = SeqIO.read(archivo_salida, "fasta")
    return str(registro.seq).upper()

def calcular_contenido_gc(secuencia):
    """Calcula el porcentaje termodinamico de Guanina-Citosina."""
    g = secuencia.count('G')
    c = secuencia.count('C')
    return ((g + c) / len(secuencia)) * 100

def disenar_grnas(secuencia_adn):
    """Escanea el genoma buscando candidatos a gRNA utilizando el patron PAM."""
    print("\nIniciando escaneo enzimatico (Busqueda de patrones PAM)...")
    
    patron = re.compile(r'(?=([ACGT]{20})([ACGT]GG))')
    candidatos = []
    
    for match in patron.finditer(secuencia_adn):
        guia = match.group(1)
        pam = match.group(2)
        posicion = match.start()
        
        porcentaje_gc = calcular_contenido_gc(guia)
        
        if 40.0 <= porcentaje_gc <= 60.0:
            candidatos.append({
                'posicion': posicion,
                'guia': guia,
                'pam': pam,
                'gc_porcentaje': round(porcentaje_gc, 2)
            })
            
    return candidatos

def exportar_resultados_csv(candidatos, archivo_salida):
    """
    Cristaliza los resultados obtenidos en un archivo CSV estructurado.
    """
    if not candidatos:
        print("No hay candidatos viables para exportar.")
        return

    # Extraemos los nombres de las columnas directamente de las llaves de nuestro diccionario
    encabezados = ['posicion', 'guia', 'pam', 'gc_porcentaje']

    try:
        # Abrimos el archivo en modo escritura. 
        # newline='' previene saltos de linea dobles en sistemas Windows.
        with open(archivo_salida, mode='w', newline='', encoding='utf-8') as archivo:
            # DictWriter es una herramienta especializada para escribir diccionarios
            escritor = csv.DictWriter(archivo, fieldnames=encabezados)
            
            # Escribimos la fila superior con los nombres de las variables
            escritor.writeheader()
            
            # Iteramos sobre nuestros resultados y los escribimos fila por fila
            for candidato in candidatos:
                escritor.writerow(candidato)
                
        print(f"\nResultados cristalizados exitosamente en: {archivo_salida}")
    except Exception as e:
        print(f"Error durante la sintesis del archivo de salida: {e}")

def main():
    id_gen = "NM_000546.6" 
    ruta_archivo_datos = os.path.join("data", f"{id_gen}.fasta")
    
    # Definimos donde guardaremos el producto final
    ruta_archivo_salida = os.path.join("output", f"{id_gen}_gRNAs_optimos.csv")
    
    print("--- Fase 1: Adquisicion de Datos Genomicos ---")
    secuencia_adn = obtener_secuencia_ncbi(id_gen, ruta_archivo_datos)
    
    if secuencia_adn:
        print("--- Fase 2: Diseno de ARN guia ---")
        resultados = disenar_grnas(secuencia_adn)
        
        total_hallazgos = len(resultados)
        print(f"Analisis completado. Se encontraron {total_hallazgos} gRNAs optimos (40%-60% GC).")
        
        if total_hallazgos > 0:
            print("--- Fase 3: Exportacion de Resultados ---")
            exportar_resultados_csv(resultados, ruta_archivo_salida)

if __name__ == "__main__":
    main()