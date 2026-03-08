import os

def preparar_entorno_trabajo():
    """
    Sintetiza la estructura de directorios y archivos base para el proyecto CRISPR.
    """
    # 1. Definimos los compartimentos (carpetas)
    directorios = ["data", "src", "output"]
    
    # 2. Definimos los reactivos y bitacoras (archivos)
    archivos = [
        ".gitignore",
        "requirements.txt",
        "README.md",
        "src/__init__.py",
        "src/main.py"
    ]
    
    print("Creando entorno de trabajo")
    
    # 3. Ensamblamos los directorios
    for directorio in directorios:
        # exist_ok=True evita errores si la carpeta ya existe, 
        # similar a un buffer quimico que mantiene la estabilidad.
        os.makedirs(directorio, exist_ok=True)
        print(f"Directorio verificado/creado: {directorio}/")
        
    # 4. Sintetizamos los archivos vacios
    for archivo in archivos:
        # Abrimos el archivo en modo escritura ("w"). 
        # Al no escribir nada y cerrarlo, creamos un archivo en blanco.
        with open(archivo, "w") as f:
            pass 
        print(f"Archivo base generado: {archivo}")
        
    print("\nEstructura del proyecto generada con exito.")

if __name__ == "__main__":
    preparar_entorno_trabajo()