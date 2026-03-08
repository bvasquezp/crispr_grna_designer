import os
import re
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Entrez, SeqIO

Entrez.email = "tu_correo@ejemplo.com"

def obtener_secuencia_ncbi(accession_id, archivo_salida):
    if not os.path.exists(archivo_salida):
        print(f"[*] Conectando con NCBI para descargar {accession_id}...")
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
            with open(archivo_salida, "w") as out_handle:
                out_handle.write(handle.read())
            handle.close()
            print(f"[+] Descarga exitosa. Guardado en: {archivo_salida}")
        except Exception as e:
            print(f"[-] Error de conexion con NCBI: {e}")
            return None
    else:
        print(f"[*] Secuencia encontrada en cache local: {archivo_salida}")
    
    try:
        registro = SeqIO.read(archivo_salida, "fasta")
        return str(registro.seq).upper()
    except Exception as e:
         print(f"[-] Error parseando archivo FASTA: {e}")
         return None

def calcular_contenido_gc(secuencia):
    if not secuencia: return 0.0
    g = secuencia.count('G')
    c = secuencia.count('C')
    return ((g + c) / len(secuencia)) * 100

def disenar_grnas(secuencia_adn, pam_target="NGG"):
    print(f"[*] Escaneando genoma en busca de motivos PAM ({pam_target})...")
    regex_pam = pam_target.replace('N', '[ACGT]')
    patron = re.compile(rf'(?=([ACGT]{{20}})({regex_pam}))')
    
    candidatos = []
    for match in patron.finditer(secuencia_adn):
        guia = match.group(1)
        pam = match.group(2)
        posicion = match.start()
        porcentaje_gc = calcular_contenido_gc(guia)
        
        if 40.0 <= porcentaje_gc <= 60.0:
            candidatos.append({
                'Posicion': posicion,
                'Secuencia_Guia': guia,
                'PAM': pam,
                'GC_Porcentaje': round(porcentaje_gc, 2)
            })
    return candidatos

def generar_visualizacion(df_candidatos, id_gen, ruta_grafico):
    if df_candidatos.empty:
        print("[-] No hay datos para visualizar.")
        return

    print("[*] Generando visualizacion de datos termodinamicos...")
    plt.figure(figsize=(12, 6))
    
    scatter = plt.scatter(
        df_candidatos['Posicion'], 
        df_candidatos['GC_Porcentaje'], 
        c=df_candidatos['GC_Porcentaje'], 
        cmap='viridis',
        alpha=0.7,
        edgecolors='w',
        s=50
    )
    
    cbar = plt.colorbar(scatter)
    cbar.set_label('Contenido GC (%)')
    
    plt.axhline(y=40, color='red', linestyle='--', alpha=0.5, label='Limite Inferior (40%)')
    plt.axhline(y=60, color='green', linestyle='--', alpha=0.5, label='Limite Superior (60%)')
    
    plt.title(f'Distribucion de gRNAs Optimos en el gen {id_gen}', fontsize=14, pad=15)
    plt.xlabel('Posicion en la Secuencia (pb)', fontsize=12)
    plt.ylabel('Contenido GC (%)', fontsize=12)
    
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(loc='lower right')
    plt.tight_layout()
    
    plt.savefig(ruta_grafico, dpi=300)
    plt.close()
    print(f"[+] Grafico exportado exitosamente en: {ruta_grafico}")

def main():
    parser = argparse.ArgumentParser(description="Disenador CRISPR-Cas9.")
    parser.add_argument("-i", "--id", type=str, required=True, help="Accession ID del NCBI")
    args = parser.parse_args()
    id_gen = args.id

    os.makedirs("data", exist_ok=True)
    os.makedirs("output", exist_ok=True)
    
    ruta_fasta = os.path.join("data", f"{id_gen}.fasta")
    ruta_csv = os.path.join("output", f"{id_gen}_gRNAs.csv")
    ruta_grafico = os.path.join("output", f"{id_gen}_distribucion_GC.png")
    
    print("\n" + "="*60)
    print(f"  Analisis CRISPR-Cas9 Avanzado: {id_gen}")
    print("="*60 + "\n")
    
    secuencia_adn = obtener_secuencia_ncbi(id_gen, ruta_fasta)
    
    if secuencia_adn:
        print(f"[+] Longitud de secuencia analizada: {len(secuencia_adn)} pb.")
        resultados_brutos = disenar_grnas(secuencia_adn)
        df_resultados = pd.DataFrame(resultados_brutos)
        
        total_hallazgos = len(df_resultados)
        print(f"[*] Analisis completado. Se hallaron {total_hallazgos} candidatos optimos.")
        
        if total_hallazgos > 0:
            df_resultados.to_csv(ruta_csv, index=False)
            print(f"[+] Tabla de resultados exportada en: {ruta_csv}")
            generar_visualizacion(df_resultados, id_gen, ruta_grafico)
        else:
             print("[-] No se encontraron candidatos viables.")
    else:
         print("[-] Proceso abortado.")

# ESTE ES EL MOTOR DE ARRANQUE. DEBE ESTAR AL FINAL DEL ARCHIVO.
if __name__ == "__main__":
    main()