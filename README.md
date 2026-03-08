# CRISPR-Cas9 gRNA Designer


## Descripción del Proyecto
Este repositorio contiene una herramienta bioinformá
desarrollada en Python para el diseño automatizado d
de ARN guía (gRNA) orientadas al sistema de edición 
CRISPR-Cas9 (variante SpCas9). 


El algoritmo extrae secuencias genómicas directament
base de datos de nucleótidos del NCBI (National Cent
Biotechnology Information) y aplica expresiones regu
localizar secuencias adyacentes al motivo PAM (`NGG`
Posteriormente, evalúa la viabilidad termodinámica d
candidato calculando su porcentaje de Guanina-Citosi
filtrando y exportando únicamente aquellos con un ra
de hibridación (40% - 60%).


## Características Principales
* **Extracción Automatizada:** Integración con la AP
(Biopython) para descargas en formato FASTA.
* **Escaneo de Alta Precisión:** Uso de algoritmos d
con solapamiento (Lookahead Regex) para identificar 
sitios de corte posibles.
* **Filtro Termodinámico:** Cálculo estructural del 
para asegurar la estabilidad de la guía molecular.
* **Exportación de Datos:** Generación de reportes t
formato CSV listos para análisis estadístico.


## Estructura del Repositorio
```text
crispr_grna_designer/
├── data/                    # Archivos FASTA descar
el NCBI
├── src/                     # Código fuente princip
│   └── main.py              # Script de ejecución
├── output/                  # Resultados exportados
CSV)
├── .gitignore               # Reglas de exclusión p
control de versiones
├── requirements.txt         # Dependencias del proy
└── README.md                # Documentación técnica
