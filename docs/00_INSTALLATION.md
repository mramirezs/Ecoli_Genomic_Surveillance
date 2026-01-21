# 🛠️ Guía de Instalación y Configuración
### Bacterial Genomics Pipeline

---

## 📋 Tabla de Contenidos

1. [Requisitos Previos](#-requisitos-previos)
2. [Instalación de Miniforge/Mamba](#-instalación-de-miniforgemamba)
3. [Configuración de Canales Bioconda](#-configuración-de-canales-bioconda)
4. [Creación de Ambientes Conda](#-creación-de-ambientes-conda)
5. [Descarga de Bases de Datos](#-descarga-de-bases-de-datos)
6. [Verificación de Instalación](#-verificación-de-instalación)
7. [Configuración del Proyecto](#-configuración-del-proyecto)
8. [Solución de Problemas](#-solución-de-problemas)

---

## ⚙️ Requisitos Previos

### Sistema Operativo

✅ **Sistemas Soportados:**
- Ubuntu 20.04 LTS o superior
- Debian 10+
- CentOS 7+
- Rocky Linux 8+
- Cualquier distribución Linux moderna

❌ **No Soportado:**
- Windows (usar WSL2)
- macOS (algunas herramientas bioinformáticas no disponibles)

### Hardware Recomendado

| Componente | Mínimo | Recomendado | Óptimo |
|------------|--------|-------------|--------|
| **CPU** | 4 cores | 8 cores | 16+ cores |
| **RAM** | 16 GB | 32 GB | 64+ GB |
| **Almacenamiento** | 100 GB libres | 200 GB libres | SSD 500 GB |
| **Red** | 10 Mbps | 100 Mbps | 1 Gbps |

### Software Base Requerido

Antes de empezar, verifica que tengas instalado:

```bash
# Verificar bash
bash --version
# Requerido: bash 4.0+

# Verificar git
git --version
# Requerido: git 2.0+

# Verificar wget o curl
wget --version
curl --version
# Al menos uno de los dos

# Verificar permisos de escritura
cd ~
mkdir -p test_dir && rm -rf test_dir && echo "✓ Permisos OK"
```

---

## 📥 Instalación de Miniforge/Mamba

### ¿Por qué Miniforge y no Anaconda?

- ✅ **Más rápido:** Mamba resuelve dependencias 10-20x más rápido
- ✅ **Gratis y libre:** No requiere licencia comercial
- ✅ **Bioconda por defecto:** Canal principal para bioinformática
- ✅ **Menor tamaño:** Solo paquetes esenciales

### Paso 1: Descargar Miniforge

```bash
# Ir al directorio home
cd ~

# Descargar instalador para Linux x86_64
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh

# Verificar descarga
ls -lh Miniforge3-Linux-x86_64.sh
# Debe mostrar archivo de ~70-80 MB
```

### Paso 2: Instalar Miniforge

```bash
# Dar permisos de ejecución
chmod +x Miniforge3-Linux-x86_64.sh

# Instalar en modo batch (sin preguntas interactivas)
bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3

# Nota: 
# -b = batch mode (sin confirmaciones)
# -p = path de instalación
```

### Paso 3: Inicializar Conda

```bash
# Inicializar conda para bash
$HOME/miniforge3/bin/conda init bash

# Recargar configuración de bash
source ~/.bashrc

# Verificar instalación
conda --version
mamba --version

# Salida esperada:
# conda 24.x.x
# mamba 1.x.x
```

### Paso 4: Configuración Inicial de Conda

```bash
# Desactivar activación automática del ambiente base
conda config --set auto_activate_base false

# Configurar mamba como solver por defecto
conda config --set solver libmamba

# Verificar configuración
conda config --show-sources
```

**🎯 Verificación:**
```bash
# Después de reiniciar terminal, NO debería aparecer (base)
# antes del prompt

# Si aparece (base), ejecutar:
conda deactivate
```

---

## 🔧 Configuración de Canales Bioconda

### ¿Qué son los Canales?

Los canales son repositorios de paquetes. Para bioinformática necesitamos:
- **conda-forge:** Paquetes científicos generales
- **bioconda:** Herramientas bioinformáticas
- **defaults:** Paquetes base de conda

### Configurar Prioridad de Canales

```bash
# Agregar canales en orden de prioridad
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Establecer prioridad estricta (IMPORTANTE)
conda config --set channel_priority strict

# Verificar configuración
conda config --show channels
```

**Salida esperada:**
```yaml
channels:
  - conda-forge
  - bioconda
  - defaults
```

### ¿Por qué channel_priority strict?

`strict` asegura que conda/mamba use paquetes del canal de mayor prioridad primero, evitando conflictos de versiones entre canales.

---

## 🐍 Creación de Ambientes Conda

### ¿Por qué 3 Ambientes Separados?

Algunas herramientas bioinformáticas tienen **conflictos de dependencias** entre sí:
- **Prokka** requiere versiones específicas de Perl
- **RGI** necesita Python 3.11
- **SPAdes/Unicycler** funcionan mejor con Python 3.10

Por eso creamos 3 ambientes especializados:

| Ambiente | Propósito | Herramientas Principales |
|----------|-----------|-------------------------|
| `bact_main` | Pipeline principal | FastQC, SPAdes, BWA, Flye, Unicycler, AMRFinder |
| `bact_amr` | Anotación y AMR | Prokka, Abricate |
| `bact_rgi` | AMR avanzado | RGI (CARD database) |

---

### 🧬 Ambiente 1: `bact_main` (Principal)

Este ambiente contiene todas las herramientas para QC, ensamblaje, mapeo y detección básica de AMR.

#### Crear Ambiente Base

```bash
# Crear ambiente con Python 3.10
mamba create -n bact_main -c conda-forge -c bioconda \
  python=3.10 pip pigz openjdk=11 -y

# Tiempo estimado: 2-3 minutos
# Tamaño: ~500 MB
```

**Explicación de paquetes base:**
- `python=3.10`: Compatible con mayoría de herramientas
- `pip`: Para paquetes Python adicionales
- `pigz`: Compresión paralela (más rápido que gzip)
- `openjdk=11`: Requerido por algunas herramientas Java

#### Activar Ambiente

```bash
conda activate bact_main

# El prompt debe cambiar a:
# (bact_main) usuario@host:~$
```

#### Instalar Herramientas de Control de Calidad

```bash
# FastQC, MultiQC, fastp (Illumina)
mamba install -c bioconda fastqc multiqc fastp -y

# NanoPlot, Filtlong (Nanopore)
mamba install -c bioconda nanoplot filtlong -y

# Tiempo estimado: 3-5 minutos
```

**Verificar instalación:**
```bash
fastqc --version    # v0.12.1
multiqc --version   # v1.14
fastp --version     # 0.23.4
NanoPlot --version  # 1.41.0
filtlong --version  # v0.2.1
```

#### Instalar Herramientas de Mapeo y Análisis de Variantes

```bash
# BWA (Illumina), Minimap2 (Nanopore)
mamba install -c bioconda bwa minimap2 -y

# Samtools, BCFtools, BEDtools
mamba install -c bioconda samtools bcftools bedtools -y

# BLAST (para búsquedas de homología)
mamba install -c bioconda blast -y

# Tiempo estimado: 3-5 minutos
```

**Verificar instalación:**
```bash
bwa 2>&1 | head -3           # BWA para mapeo Illumina
minimap2 --version           # 2.24-r1122
samtools --version           # 1.17
bcftools --version           # 1.17
```

#### Instalar Ensambladores

```bash
# SPAdes (Illumina)
mamba install -c bioconda spades -y

# Flye (Nanopore)
mamba install -c bioconda flye -y

# Unicycler (Híbrido)
mamba install -c bioconda unicycler -y

# QUAST (Evaluación de calidad)
mamba install -c bioconda quast -y

# Bandage (Visualización de gráficos)
mamba install -c bioconda bandage -y

# Tiempo estimado: 5-8 minutos
```

**Verificar instalación:**
```bash
spades.py --version        # 3.15.5
flye --version             # 2.9.1
unicycler --version        # 0.5.0
quast.py --version         # 5.2.0
Bandage --version          # 0.8.1
```

#### Instalar Herramientas AMR y Typing

```bash
# AMRFinderPlus (NCBI)
mamba install -c bioconda ncbi-amrfinderplus -y

# Barrnap (rRNA prediction)
mamba install -c bioconda barrnap -y

# MLST (Multi-Locus Sequence Typing)
mamba install -c bioconda mlst -y

# Tiempo estimado: 2-3 minutos
```

**Verificar instalación:**
```bash
amrfinder --version        # 3.11.4
barrnap --version          # 0.9
mlst --version             # 2.23.0
```

#### Instalar Herramientas Adicionales

```bash
# seqtk (manipulación de secuencias)
mamba install -c bioconda seqtk -y

# Kraken2 (clasificación taxonómica - opcional)
mamba install -c bioconda kraken2 -y

# Tiempo estimado: 2-3 minutos
```

#### Actualizar Base de Datos AMRFinderPlus

```bash
# Crear directorio para base de datos
mkdir -p ~/bacterial_genomics/databases/amrfinder_db

# Descargar base de datos AMRFinderPlus
amrfinder_update --database ~/bacterial_genomics/databases/amrfinder_db

# Tiempo estimado: 5-10 minutos
# Tamaño: ~700 MB
```

**Verificar base de datos:**
```bash
amrfinder --database ~/bacterial_genomics/databases/amrfinder_db --version

# Salida esperada:
# Database version: 2024-01-31.1
```

#### Actualizar Base de Datos MLST

```bash
# Actualizar esquemas MLST
mlst --list

# Esto descarga esquemas para ~150 especies
# Tiempo estimado: 2-3 minutos

# Verificar que Klebsiella pneumoniae está disponible
mlst --list | grep pneumoniae
```

**✅ Ambiente `bact_main` completo**

```bash
# Verificación final
echo "=== VERIFICACIÓN BACT_MAIN ==="
which fastqc
which spades.py
which bwa
which amrfinder
echo "✓ Ambiente bact_main instalado correctamente"

# Desactivar ambiente
conda deactivate
```

---

### 🦠 Ambiente 2: `bact_amr` (Anotación y AMR)

Este ambiente está dedicado a **Prokka** y **Abricate**, que requieren versiones específicas de Perl.

#### Crear Ambiente

```bash
# Crear ambiente con Python 3.9 y herramientas Perl
mamba create -n bact_amr -c conda-forge -c bioconda \
  python=3.9 prokka abricate -y

# Tiempo estimado: 5-7 minutos
# Tamaño: ~800 MB (incluye dependencias Perl)
```

**¿Qué se instala?**
- **Prokka:** Anotación rápida de genomas bacterianos
- **Abricate:** Screening de genes de resistencia (múltiples DBs)
- Dependencias: Perl, BioPerl, BLAST+, etc.

#### Activar y Configurar

```bash
# Activar ambiente
conda activate bact_amr

# Configurar bases de datos de Abricate
abricate --setupdb

# Tiempo estimado: 3-5 minutos
# Descarga: ~150 MB
```

#### Verificar Bases de Datos Disponibles

```bash
# Listar bases de datos
abricate --list

# Salida esperada:
# DATABASE       SEQUENCES  DBTYPE  DATE
# argannot       2223       nucl    2023-Apr-17
# card           3094       nucl    2023-Aug-22
# ecoh           597        nucl    2023-Apr-17
# ecoli_vf       2701       nucl    2023-Apr-17
# megares        7635       nucl    2023-Apr-17
# ncbi           5386       nucl    2023-Jul-13
# plasmidfinder  460        nucl    2023-Apr-17
# resfinder      3077       nucl    2023-Apr-17
# vfdb           2597       nucl    2023-Apr-17
```

**Bases de datos importantes:**
- **card:** CARD (Comprehensive Antibiotic Resistance Database)
- **resfinder:** ResFinder (validado clínicamente)
- **ncbi:** NCBI AMR database
- **vfdb:** Factores de virulencia
- **plasmidfinder:** Detección de plásmidos

#### Verificar Prokka

```bash
# Verificar instalación
prokka --version

# Salida esperada:
# prokka 1.14.6

# Ver opciones disponibles
prokka --listdb

# Debe mostrar bases de datos de genes, proteínas, etc.
```

**✅ Ambiente `bact_amr` completo**

```bash
echo "=== VERIFICACIÓN BACT_AMR ==="
prokka --version
abricate --version
abricate --list | wc -l  # Debe mostrar ~9 bases de datos
echo "✓ Ambiente bact_amr instalado correctamente"

# Desactivar ambiente
conda deactivate
```

---

### 🧪 Ambiente 3: `bact_rgi` (AMR Avanzado)

Este ambiente está dedicado a **RGI** (Resistance Gene Identifier) con la base de datos **CARD**.

#### Crear Ambiente

```bash
# Crear ambiente con Python 3.11 (requerido por RGI)
mamba create -n bact_rgi -c conda-forge -c bioconda \
  python=3.11 rgi -y

# Tiempo estimado: 3-4 minutos
# Tamaño: ~400 MB
```

#### Activar Ambiente

```bash
conda activate bact_rgi

# Verificar instalación
rgi main --version

# Salida esperada:
# 6.0.2
```

#### Descargar y Cargar Base de Datos CARD

```bash
# Crear directorio
mkdir -p ~/bacterial_genomics/databases/card

# Ir al directorio
cd ~/bacterial_genomics/databases/card

# Descargar base de datos CARD
wget https://card.mcmaster.ca/latest/data

# Descomprimir
tar -xvf data

# Cargar base de datos en RGI (modo local)
rgi load --card_json card.json --local

# Tiempo estimado: 2-3 minutos
# Tamaño: ~50 MB

# Volver al directorio inicial
cd ~
```

#### Verificar Carga de Base de Datos

```bash
# Verificar versión de base de datos
rgi database --version --local

# Salida esperada:
# Database version: 3.2.7
# Database name: CARD
```

#### Actualizar Base de Datos CARD (Opcional)

```bash
# Para actualizar en el futuro
cd ~/bacterial_genomics/databases/card

# Descargar última versión
wget -O data_new https://card.mcmaster.ca/latest/data

# Reemplazar y recargar
mv data data_old
mv data_new data
tar -xvf data
rgi load --card_json card.json --local --local_database
```

**✅ Ambiente `bact_rgi` completo**

```bash
echo "=== VERIFICACIÓN BACT_RGI ==="
rgi main --version
rgi database --version --local
echo "✓ Ambiente bact_rgi instalado correctamente"

# Desactivar ambiente
conda deactivate
```

---

## 🚀 Script de Instalación Automatizada

Para facilitar la instalación, puedes usar este script que configura los 3 ambientes automáticamente.

### Crear Script de Instalación

```bash
# Crear directorio de scripts
mkdir -p ~/bacterial_genomics/scripts

# Crear script
cat > ~/bacterial_genomics/scripts/setup_environments.sh << 'EOF'
#!/bin/bash

set -e  # Salir si hay errores

echo "========================================"
echo "Instalación de Ambientes - Bacterial Genomics Pipeline"
echo "========================================"
echo ""

# Verificar que conda/mamba estén instalados
if ! command -v mamba &> /dev/null; then
    echo "❌ ERROR: mamba no está instalado"
    echo "Por favor instala Miniforge primero"
    exit 1
fi

echo "✓ mamba encontrado: $(mamba --version)"
echo ""

# Crear directorio de bases de datos
mkdir -p ~/bacterial_genomics/databases/{amrfinder_db,card}

#######################################
# AMBIENTE 1: bact_main
#######################################
echo "========================================" 
echo "[1/3] Creando ambiente: bact_main"
echo "========================================"

# Crear ambiente base
mamba create -n bact_main -c conda-forge -c bioconda \
  python=3.10 pip pigz openjdk=11 -y

# Activar ambiente
eval "$(conda shell.bash hook)"
conda activate bact_main

# Instalar herramientas QC
echo "  [1.1] Instalando herramientas QC..."
mamba install -c bioconda fastqc multiqc fastp nanoplot filtlong -y

# Instalar herramientas mapeo
echo "  [1.2] Instalando herramientas de mapeo..."
mamba install -c bioconda bwa minimap2 samtools bcftools bedtools blast -y

# Instalar ensambladores
echo "  [1.3] Instalando ensambladores..."
mamba install -c bioconda spades flye unicycler quast bandage -y

# Instalar AMR y typing
echo "  [1.4] Instalando AMRFinder y MLST..."
mamba install -c bioconda ncbi-amrfinderplus barrnap mlst -y

# Instalar adicionales
echo "  [1.5] Instalando herramientas adicionales..."
mamba install -c bioconda seqtk kraken2 -y

# Descargar base de datos AMRFinder
echo "  [1.6] Descargando base de datos AMRFinderPlus..."
amrfinder_update --database ~/bacterial_genomics/databases/amrfinder_db

# Actualizar MLST
echo "  [1.7] Actualizando esquemas MLST..."
mlst --list > /dev/null 2>&1

conda deactivate
echo "✓ Ambiente bact_main completado"
echo ""

#######################################
# AMBIENTE 2: bact_amr
#######################################
echo "========================================"
echo "[2/3] Creando ambiente: bact_amr"
echo "========================================"

mamba create -n bact_amr -c conda-forge -c bioconda \
  python=3.9 prokka abricate -y

conda activate bact_amr

echo "  [2.1] Configurando bases de datos Abricate..."
abricate --setupdb

conda deactivate
echo "✓ Ambiente bact_amr completado"
echo ""

#######################################
# AMBIENTE 3: bact_rgi
#######################################
echo "========================================"
echo "[3/3] Creando ambiente: bact_rgi"
echo "========================================"

mamba create -n bact_rgi -c conda-forge -c bioconda \
  python=3.11 rgi -y

conda activate bact_rgi

echo "  [3.1] Descargando base de datos CARD..."
cd ~/bacterial_genomics/databases/card
wget -q https://card.mcmaster.ca/latest/data
tar -xf data
rgi load --card_json card.json --local
cd ~

conda deactivate
echo "✓ Ambiente bact_rgi completado"
echo ""

#######################################
# VERIFICACIÓN FINAL
#######################################
echo "========================================"
echo "Verificación de Instalación"
echo "========================================"

conda activate bact_main
echo "✓ bact_main:"
echo "  - FastQC: $(fastqc --version 2>&1 | head -1)"
echo "  - SPAdes: $(spades.py --version 2>&1 | head -1)"
echo "  - AMRFinder: $(amrfinder --version 2>&1)"
conda deactivate

conda activate bact_amr
echo "✓ bact_amr:"
echo "  - Prokka: $(prokka --version 2>&1 | head -1)"
echo "  - Abricate: $(abricate --version 2>&1)"
conda deactivate

conda activate bact_rgi
echo "✓ bact_rgi:"
echo "  - RGI: $(rgi main --version 2>&1)"
conda deactivate

echo ""
echo "========================================"
echo "✓ INSTALACIÓN COMPLETADA"
echo "========================================"
echo ""
echo "Ambientes creados:"
echo "  1. bact_main  - Pipeline principal"
echo "  2. bact_amr   - Anotación y AMR"
echo "  3. bact_rgi   - AMR avanzado (CARD)"
echo ""
echo "Para activar un ambiente:"
echo "  conda activate bact_main"
echo ""
echo "Bases de datos en:"
echo "  ~/bacterial_genomics/databases/"
echo ""
echo "Siguiente paso:"
echo "  bash ~/bacterial_genomics/scripts/verify_installation.sh"
echo ""
EOF

# Dar permisos de ejecución
chmod +x ~/bacterial_genomics/scripts/setup_environments.sh
```

### Ejecutar Instalación Automatizada

```bash
# Ejecutar script
bash ~/bacterial_genomics/scripts/setup_environments.sh

# Tiempo total estimado: 45-60 minutos
# Depende de velocidad de internet y CPU
```

---

## ✅ Verificación de Instalación

### Script de Verificación

```bash
# Crear script de verificación
cat > ~/bacterial_genomics/scripts/verify_installation.sh << 'EOF'
#!/bin/bash

echo "========================================"
echo "Verificación de Instalación"
echo "Bacterial Genomics Pipeline"
echo "========================================"
echo ""

# Función para verificar comando
check_tool() {
    local env=$1
    local tool=$2
    local cmd=$3
    
    conda activate $env 2>/dev/null
    if command -v $tool &> /dev/null; then
        version=$($cmd 2>&1 | head -1)
        echo "  ✓ $tool: OK ($version)"
        status=0
    else
        echo "  ❌ $tool: NO ENCONTRADO"
        status=1
    fi
    conda deactivate 2>/dev/null
    return $status
}

errors=0

# Verificar ambiente bact_main
echo "[Ambiente: bact_main]"
check_tool bact_main fastqc "fastqc --version" || ((errors++))
check_tool bact_main fastp "fastp --version" || ((errors++))
check_tool bact_main bwa "bwa 2>&1 | head -1" || ((errors++))
check_tool bact_main samtools "samtools --version" || ((errors++))
check_tool bact_main spades.py "spades.py --version" || ((errors++))
check_tool bact_main flye "flye --version" || ((errors++))
check_tool bact_main unicycler "unicycler --version" || ((errors++))
check_tool bact_main quast.py "quast.py --version" || ((errors++))
check_tool bact_main amrfinder "amrfinder --version" || ((errors++))
check_tool bact_main mlst "mlst --version" || ((errors++))
echo ""

# Verificar ambiente bact_amr
echo "[Ambiente: bact_amr]"
check_tool bact_amr prokka "prokka --version" || ((errors++))
check_tool bact_amr abricate "abricate --version" || ((errors++))
echo ""

# Verificar ambiente bact_rgi
echo "[Ambiente: bact_rgi]"
check_tool bact_rgi rgi "rgi main --version" || ((errors++))
echo ""

# Verificar bases de datos
echo "[Bases de Datos]"
if [ -d ~/bacterial_genomics/databases/amrfinder_db ]; then
    echo "  ✓ AMRFinderPlus DB: Instalada"
else
    echo "  ❌ AMRFinderPlus DB: NO ENCONTRADA"
    ((errors++))
fi

conda activate bact_amr 2>/dev/null
db_count=$(abricate --list 2>/dev/null | wc -l)
if [ $db_count -gt 5 ]; then
    echo "  ✓ Abricate DBs: $db_count bases disponibles"
else
    echo "  ❌ Abricate DBs: Incompletas"
    ((errors++))
fi
conda deactivate 2>/dev/null

if [ -f ~/bacterial_genomics/databases/card/card.json ]; then
    echo "  ✓ CARD DB: Instalada"
else
    echo "  ❌ CARD DB: NO ENCONTRADA"
    ((errors++))
fi
echo ""

# Resumen final
echo "========================================"
if [ $errors -eq 0 ]; then
    echo "✓ TODAS LAS VERIFICACIONES PASARON"
    echo "El sistema está listo para usar"
else
    echo "❌ SE ENCONTRARON $errors ERRORES"
    echo "Revisa los mensajes arriba"
fi
echo "========================================"
echo ""

exit $errors
EOF

chmod +x ~/bacterial_genomics/scripts/verify_installation.sh
```

### Ejecutar Verificación

```bash
bash ~/bacterial_genomics/scripts/verify_installation.sh
```

**Salida esperada si todo está bien:**
```
========================================
Verificación de Instalación
Bacterial Genomics Pipeline
========================================

[Ambiente: bact_main]
  ✓ fastqc: OK (FastQC v0.12.1)
  ✓ fastp: OK (fastp 0.23.4)
  ✓ bwa: OK (Version: 0.7.17-r1188)
  ✓ samtools: OK (samtools 1.17)
  ✓ spades.py: OK (SPAdes v3.15.5)
  ✓ flye: OK (2.9.1-b1780)
  ✓ unicycler: OK (Unicycler v0.5.0)
  ✓ quast.py: OK (QUAST v5.2.0)
  ✓ amrfinder: OK (3.11.4)
  ✓ mlst: OK (mlst 2.23.0)

[Ambiente: bact_amr]
  ✓ prokka: OK (prokka 1.14.6)
  ✓ abricate: OK (abricate 1.0.1)

[Ambiente: bact_rgi]
  ✓ rgi: OK (6.0.2)

[Bases de Datos]
  ✓ AMRFinderPlus DB: Instalada
  ✓ Abricate DBs: 9 bases disponibles
  ✓ CARD DB: Instalada

========================================
✓ TODAS LAS VERIFICACIONES PASARON
El sistema está listo para usar
========================================
```

---

## 📦 Exportar Ambientes (Reproducibilidad)

### ¿Por qué Exportar?

Exportar los ambientes te permite:
- ✅ Recrear la instalación en otro servidor
- ✅ Compartir configuración con colaboradores
- ✅ Documentar versiones exactas de software
- ✅ Reproducir resultados

### Exportar los 3 Ambientes

```bash
# Crear directorio para archivos YAML
mkdir -p ~/bacterial_genomics/envs

# Exportar bact_main
conda activate bact_main
conda env export --no-builds > ~/bacterial_genomics/envs/bact_main.yml
conda deactivate

# Exportar bact_amr
conda activate bact_amr
conda env export --no-builds > ~/bacterial_genomics/envs/bact_amr.yml
conda deactivate

# Exportar bact_rgi
conda activate bact_rgi
conda env export --no-builds > ~/bacterial_genomics/envs/bact_rgi.yml
conda deactivate

echo "✓ Ambientes exportados en: ~/bacterial_genomics/envs/"
ls -lh ~/bacterial_genomics/envs/
```

### Recrear Ambientes desde YAML (En otro servidor
