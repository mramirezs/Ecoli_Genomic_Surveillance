# 🧬 Pipeline de Vigilancia Genómica y Análisis de Resistencia Antimicrobiana en Bacterias

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Bioinformatics](https://img.shields.io/badge/Bioinformatics-Pipeline-blue.svg)]()
[![Status](https://img.shields.io/badge/Status-Production-green.svg)]()

Este repositorio documenta un flujo de trabajo bioinformático completo para el análisis de genomas bacterianos clínicos utilizando datos de secuenciación de nueva generación (NGS). El pipeline integra tres estrategias de ensamblaje complementarias: **Ensamblaje con Illumina**, **Ensamblaje con Nanopore** y **Ensamblaje Híbrido (Illumina + Nanopore)**, junto con detección exhaustiva de genes de resistencia a antimicrobianos (AMR) y análisis de variantes genómicas.

---

## 📋 Tabla de Contenidos

- [⚠️ Antes de Comenzar](#️-antes-de-comenzar)
- [Características del Pipeline](#-características-del-pipeline)
- [Estructura del Proyecto](#-estructura-del-proyecto)
- [Requisitos del Sistema](#-requisitos-del-sistema)
- [Instalación y Configuración](#️-instalación-y-configuración)
- [Dataset de Prueba](#-dataset-de-prueba)
- [Flujo de Trabajo](#-flujo-de-trabajo)
- [Resultados Esperados](#-resultados-esperados)
- [Interpretación de Resultados](#-interpretación-de-resultados)
- [Checklist de Validación](#-checklist-de-validación)
- [Solución de Problemas](#-solución-de-problemas)
- [Casos de Uso](#-casos-de-uso)
- [Limitaciones Conocidas](#️-limitaciones-conocidas)
- [Referencias](#-referencias)

---

## ⚠️ Antes de Comenzar

### Requisitos Previos

- [ ] **Datos de secuenciación** en formato FASTQ (Illumina y/o Nanopore)
- [ ] **~100-200 GB** de espacio libre en disco por muestra
- [ ] **Sistema Linux/Unix** (Ubuntu 20.04+, CentOS 7+, o similar)
- [ ] **Acceso a internet** para descargar herramientas y bases de datos
- [ ] **Tiempo estimado**: 4-8 horas por muestra (dependiendo de hardware)

### 🚀 Inicio Rápido (Usuarios Experimentados)

```bash
# 1. Clonar repositorio
git clone https://github.com/tu-usuario/Bacterial_Genomics_Project.git
cd Bacterial_Genomics_Project

# 2. Configurar ambientes (primera vez)
bash scripts/setup_environments.sh

# 3. Preparar datos
ln -s /ruta/a/tus/datos/*.fastq.gz 00_raw_data/illumina/
ln -s /ruta/a/tus/datos/*.fastq.gz 00_raw_data/nanopore/

# 4. Ejecutar pipeline completo
bash scripts/run_full_pipeline.sh sample_name

# 5. Ver resultados
firefox 08_results/FINAL_REPORT.html
```

### 📊 ¿Qué Puedo Hacer con Este Pipeline?

✅ **Ensamblar genomas bacterianos** de alta calidad  
✅ **Identificar genes de resistencia** a antibióticos (AMR)  
✅ **Detectar variantes genómicas** (SNPs, INDELs)  
✅ **Anotar genes y funciones** biológicas  
✅ **Comparar diferentes estrategias** de ensamblaje  
✅ **Tipificar cepas** (MLST, detección de plásmidos)  
✅ **Generar reportes automatizados** para vigilancia epidemiológica  

---

## 🎯 Características del Pipeline

### Tecnologías Soportadas
- **Illumina** (lecturas cortas, paired-end): Alta precisión, ideal para SNPs/INDELs
- **Oxford Nanopore** (lecturas largas): Ensamblajes contiguos, detección de variantes estructurales
- **Híbrido** (Illumina + Nanopore): Combina precisión y continuidad

### Análisis Incluidos
- ✅ Control de calidad exhaustivo (raw y trimmed reads)
- ✅ Tres estrategias de ensamblaje independientes
- ✅ Mapeo contra genoma de referencia y llamado de variantes
- ✅ Detección de genes AMR con múltiples bases de datos
- ✅ Anotación funcional de genomas
- ✅ Evaluación de calidad de ensamblajes
- ✅ MLST typing (tipificación molecular)
- ✅ Detección de plásmidos y factores de virulencia
- ✅ Visualización y reportes integrados

---

## 📂 Estructura del Proyecto

```text
Bacterial_Genomics_Project/
├── 00_raw_data/                    # Datos crudos de secuenciación
│   ├── illumina/                   # Lecturas paired-end (R1, R2)
│   │   ├── sample_R1.fastq.gz
│   │   └── sample_R2.fastq.gz
│   └── nanopore/                   # Lecturas largas ONT
│       └── sample_ont.fastq.gz
│
├── 01_reference/                   # Genomas de referencia (opcional)
│   ├── reference.fasta
│   └── reference.gff
│
├── 02_qc/                          # Control de calidad
│   ├── 01_illumina_raw/            # FastQC de datos crudos Illumina
│   ├── 02_illumina_trimmed/        # FastQC post-trimming + reportes fastp
│   ├── 03_nanopore_raw/            # NanoPlot de datos crudos ONT
│   ├── 04_nanopore_filtered/       # NanoPlot post-filtrado
│   └── 05_multiqc/                 # Reporte consolidado MultiQC
│
├── 03_assembly/                    # Ensamblajes de novo
│   ├── 01_illumina_only/           # SPAdes (solo Illumina)
│   ├── 02_nanopore_only/           # Flye (solo Nanopore)
│   ├── 03_hybrid/                  # Unicycler (Illumina + Nanopore)
│   └── 04_quast_evaluation/        # Evaluación comparativa QUAST
│
├── 04_mapping/                     # Mapeo y análisis de variantes
│   ├── 01_illumina/                # BWA + Samtools
│   ├── 02_nanopore/                # Minimap2 + Samtools
│   └── 03_variants/                # BCFtools variant calling
│
├── 05_annotation/                  # Anotación funcional
│   ├── 01_prokka/                  # Anotación Prokka
│   └── 02_bakta/                   # Anotación Bakta (alternativa)
│
├── 06_amr_screening/               # Detección de genes AMR
│   ├── amrfinder_db/               # Base de datos local AMRFinderPlus
│   ├── 01_amrfinder/               # Resultados AMRFinderPlus (NCBI)
│   ├── 02_abricate/                # Resultados Abricate (múltiples DBs)
│   └── 03_rgi/                     # Resultados RGI/CARD
│
├── 07_typing/                      # Tipificación molecular
│   ├── mlst/                       # MLST typing
│   ├── plasmids/                   # Detección de plásmidos
│   └── virulence/                  # Factores de virulencia
│
├── 08_results/                     # Resultados consolidados y figuras
│   ├── assembly_comparison.png
│   ├── amr_summary.xlsx
│   ├── quality_dashboard.html
│   └── FINAL_REPORT.html
│
├── envs/                           # Archivos YAML de ambientes Conda
│   ├── bact_main.yml
│   ├── bact_amr.yml
│   └── bact_rgi.yml
│
├── scripts/                        # Scripts de automatización
│   ├── setup_environments.sh       # Instalación de ambientes
│   ├── verify_installation.sh      # Verificación de instalación
│   ├── run_full_pipeline.sh        # Pipeline completo
│   └── utils/                      # Scripts auxiliares
│       ├── calculate_metrics.sh
│       ├── compare_amr_tools.py
│       ├── generate_plots.py
│       └── extract_plasmids.sh
│
├── test_data/                      # Datos de prueba
│   └── ecoli_test/
│
├── logs/                           # Logs de ejecución
│   └── [timestamp]_pipeline.log
│
├── README.md                       # Este archivo
└── LICENSE                         # Licencia MIT
```

---

## 💻 Requisitos del Sistema

### Hardware Recomendado

| Componente | Mínimo | Recomendado | Óptimo |
|------------|--------|-------------|--------|
| **CPU** | 4 cores | 8 cores | 16+ cores |
| **RAM** | 16 GB | 32 GB | 64+ GB |
| **Almacenamiento** | 50 GB/muestra | 100 GB/muestra | SSD 200 GB/muestra |
| **Red** | 10 Mbps | 100 Mbps | 1 Gbps |

### Software Base
- **Sistema Operativo**: Linux/Unix (Ubuntu 20.04+, CentOS 7+, Debian 10+)
- **Shell**: Bash 4.0+
- **Git**: 2.0+
- **Conexión a internet**: Requerida para instalación inicial

### Tiempo de Ejecución Estimado

| Análisis | Hardware Mínimo | Hardware Recomendado |
|----------|-----------------|---------------------|
| QC Completo | 30-60 min | 15-30 min |
| Ensamblaje Illumina | 2-4 horas | 1-2 horas |
| Ensamblaje Nanopore | 1-2 horas | 30-60 min |
| Ensamblaje Híbrido | 4-8 horas | 2-4 horas |
| Detección AMR | 30-60 min | 15-30 min |
| **Pipeline Completo** | **8-15 horas** | **4-8 horas** |

---

## 🛠️ Instalación y Configuración

### Paso 1: Instalar Miniforge (Gestor de Paquetes)

Si aún no tienes un gestor de ambientes Conda instalado:

```bash
# Descargar Miniforge para Linux x86_64
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"

# Instalar
bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3

# Inicializar
$HOME/miniforge3/bin/conda init bash
source ~/.bashrc

# Verificar instalación
mamba --version
conda --version
```

### Paso 2: Configurar Canales de Bioconda

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

### Paso 3: Crear los Tres Ambientes Especializados

Debido a conflictos de dependencias entre herramientas bioinformáticas, el pipeline utiliza **tres ambientes Conda separados** para garantizar compatibilidad y reproducibilidad.

#### 🧬 Ambiente 1: `bact_main` (Pipeline Principal)

Contiene herramientas para QC, mapeo, ensamblaje y detección básica de AMR.

```bash
# Crear ambiente base
mamba create -n bact_main -c conda-forge -c bioconda -c defaults \
  python=3.10 pip pigz openjdk=11 -y

# Activar
conda activate bact_main

# Instalar herramientas de control de calidad
mamba install fastqc multiqc fastp nanoplot filtlong -y

# Instalar herramientas de mapeo y análisis de variantes
mamba install bwa minimap2 samtools bcftools bedtools blast -y

# Instalar ensambladores
mamba install unicycler flye spades quast bandage -y

# Instalar herramientas AMR y typing
mamba install ncbi-amrfinderplus barrnap mlst -y

# Instalar herramientas adicionales
mamba install seqtk kraken2 -y

# Configurar base de datos AMRFinderPlus (primera vez)
mkdir -p 06_amr_screening/amrfinder_db
amrfinder_update --database 06_amr_screening/amrfinder_db

# Actualizar base de datos MLST
mlst --list
```

**⏱️ Tiempo de instalación**: ~20 minutos  
**📦 Descarga de base de datos**: ~700 MB adicionales

#### 🦠 Ambiente 2: `bact_amr` (Anotación y AMR)

Dedicado a Prokka y Abricate, que requieren versiones específicas de Perl.

```bash
# Crear ambiente
mamba create -n bact_amr -c conda-forge -c bioconda -c defaults \
  python=3.9 prokka abricate -y

# Activar y configurar bases de datos
conda activate bact_amr
abricate --setupdb

# Verificar bases de datos disponibles
abricate --list
```

**⏱️ Tiempo de instalación**: ~15 minutos  
**📦 Descarga de bases de datos**: ~150 MB adicionales

#### 🧪 Ambiente 3: `bact_rgi` (AMR Avanzado)

Para RGI (Resistance Gene Identifier) con base de datos CARD.

```bash
# Crear ambiente
mamba create -n bact_rgi -c conda-forge -c bioconda -c defaults \
  python=3.11 rgi -y

# Activar
conda activate bact_rgi

# Descargar y cargar base de datos CARD
mkdir -p 06_amr_screening/rgi
cd 06_amr_screening/rgi
wget https://card.mcmaster.ca/latest/data
tar -xvf data ./card.json
rgi load --card_json card.json --local
cd ../..

# Verificar carga
rgi database --version --local
```

**⏱️ Tiempo de instalación**: ~10 minutos  
**📦 Descarga de base de datos CARD**: ~50 MB

### Paso 4: Script de Instalación Automatizada

Crea el archivo `scripts/setup_environments.sh`:

```bash
#!/bin/bash
# setup_environments.sh - Instalación automatizada de ambientes

set -euo pipefail

echo "========================================"
echo "Configuración de Ambientes Bioinformáticos"
echo "Inicio: $(date)"
echo "========================================"

# Verificar Conda/Mamba instalado
if ! command -v mamba &> /dev/null; then
    echo "ERROR: Mamba no está instalado. Instalar Miniforge primero."
    exit 1
fi

# Configurar canales
echo "[1/4] Configurando canales Conda..."
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# Crear ambiente bact_main
echo "[2/4] Creando ambiente bact_main..."
mamba create -n bact_main -c conda-forge -c bioconda \
  python=3.10 fastqc multiqc fastp nanoplot filtlong \
  bwa minimap2 samtools bcftools bedtools blast \
  unicycler flye spades quast bandage \
  ncbi-amrfinderplus barrnap mlst seqtk kraken2 -y

# Crear ambiente bact_amr
echo "[3/4] Creando ambiente bact_amr..."
mamba create -n bact_amr -c conda-forge -c bioconda \
  python=3.9 prokka abricate -y

# Crear ambiente bact_rgi
echo "[4/4] Creando ambiente bact_rgi..."
mamba create -n bact_rgi -c conda-forge -c bioconda \
  python=3.11 rgi -y

echo ""
echo "✓ Ambientes creados exitosamente"
echo ""
echo "Descargando bases de datos..."

# Configurar bases de datos
conda activate bact_main
mkdir -p 06_amr_screening/amrfinder_db
amrfinder_update --database 06_amr_screening/amrfinder_db
mlst --list > /dev/null 2>&1

conda activate bact_amr
abricate --setupdb

conda activate bact_rgi
mkdir -p 06_amr_screening/rgi
cd 06_amr_screening/rgi
wget -q https://card.mcmaster.ca/latest/data -O card_data.tar.bz2
tar -xjf card_data.tar.bz2 ./card.json
rgi load --card_json card.json --local
cd ../..

echo ""
echo "========================================"
echo "Instalación Completada Exitosamente"
echo "Fin: $(date)"
echo "========================================"
echo ""
echo "Para verificar la instalación, ejecuta:"
echo "  bash scripts/verify_installation.sh"
```

### Paso 5: Script de Verificación

Crea el archivo `scripts/verify_installation.sh`:

```bash
#!/bin/bash
# verify_installation.sh - Verificar instalación completa

set -euo pipefail

echo "========================================"
echo "Verificación de Instalación"
echo "========================================"
echo ""

ERRORS=0

# Función para verificar comando
check_command() {
    local env=$1
    local cmd=$2
    local name=$3
    
    conda activate $env 2>/dev/null
    if command -v $cmd &> /dev/null; then
        version=$($cmd --version 2>&1 | head -1 || echo "instalado")
        echo "✓ $name: OK"
    else
        echo "✗ $name: NO ENCONTRADO"
        ((ERRORS++))
    fi
}

# Verificar bact_main
echo "[Ambiente: bact_main]"
check_command bact_main fastqc "FastQC"
check_command bact_main fastp "fastp"
check_command bact_main bwa "BWA"
check_command bact_main samtools "Samtools"
check_command bact_main spades.py "SPAdes"
check_command bact_main flye "Flye"
check_command bact_main unicycler "Unicycler"
check_command bact_main quast.py "QUAST"
check_command bact_main amrfinder "AMRFinderPlus"
check_command bact_main mlst "MLST"
echo ""

# Verificar bact_amr
echo "[Ambiente: bact_amr]"
check_command bact_amr prokka "Prokka"
check_command bact_amr abricate "Abricate"
echo ""

# Verificar bact_rgi
echo "[Ambiente: bact_rgi]"
check_command bact_rgi rgi "RGI"
echo ""

# Verificar bases de datos
echo "[Bases de Datos]"
conda activate bact_main
if [ -d "06_amr_screening/amrfinder_db/latest" ]; then
    echo "✓ AMRFinderPlus DB: Instalada"
else
    echo "✗ AMRFinderPlus DB: NO ENCONTRADA"
    ((ERRORS++))
fi

conda activate bact_amr
DB_COUNT=$(abricate --list 2>/dev/null | wc -l)
if [ $DB_COUNT -gt 1 ]; then
    echo "✓ Abricate DBs: $DB_COUNT bases disponibles"
else
    echo "✗ Abricate DBs: NO ENCONTRADAS"
    ((ERRORS++))
fi

conda activate bact_rgi
if rgi database --version --local &>/dev/null; then
    echo "✓ CARD DB: Instalada"
else
    echo "✗ CARD DB: NO ENCONTRADA"
    ((ERRORS++))
fi

echo ""
echo "========================================"
if [ $ERRORS -eq 0 ]; then
    echo "✓ TODAS LAS VERIFICACIONES PASARON"
    echo "El sistema está listo para usar"
else
    echo "✗ SE ENCONTRARON $ERRORS ERRORES"
    echo "Revisa los componentes faltantes arriba"
fi
echo "========================================"

exit $ERRORS
```

### Paso 6: Ejecutar Instalación y Verificación

```bash
# Dar permisos de ejecución
chmod +x scripts/setup_environments.sh
chmod +x scripts/verify_installation.sh

# Ejecutar instalación
bash scripts/setup_environments.sh

# Verificar instalación
bash scripts/verify_installation.sh
```

### Paso 7: Exportar Ambientes (Reproducibilidad)

```bash
# Crear directorio
mkdir -p envs

# Exportar ambientes
conda activate bact_main
conda env export --no-builds > envs/bact_main.yml

conda activate bact_amr
conda env export --no-builds > envs/bact_amr.yml

conda activate bact_rgi
conda env export --no-builds > envs/bact_rgi.yml

echo "Ambientes exportados en envs/"
```

---

## 🧪 Dataset de Prueba

Para probar el pipeline sin usar tus propios datos, descarga este dataset público de *E. coli* O157:H7:

```bash
# Crear directorio de prueba
mkdir -p test_data/ecoli_o157h7
cd test_data/ecoli_o157h7

# Descargar lecturas Illumina (SRA: SRR1616829)
conda activate bact_main
fastq-dump --split-files --gzip SRR1616829

# Renombrar archivos
mv SRR1616829_1.fastq.gz ecoli_R1.fastq.gz
mv SRR1616829_2.fastq.gz ecoli_R2.fastq.gz

# Volver al directorio principal
cd ../..

# Crear enlaces simbólicos
ln -s $(pwd)/test_data/ecoli_o157h7/ecoli_R1.fastq.gz 00_raw_data/illumina/
ln -s $(pwd)/test_data/ecoli_o157h7/ecoli_R2.fastq.gz 00_raw_data/illumina/

# Descargar genoma de referencia E. coli O157:H7 Sakai
cd 01_reference
wget -O reference.fasta.gz \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/865/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz"
gunzip reference.fasta.gz
cd ..
```

### Resultados Esperados del Dataset de Prueba

| Métrica | Valor Esperado |
|---------|----------------|
| Tamaño del Ensamblaje | ~5.5 Mb |
| Número de Contigs | 10-30 |
| N50 | >200 kb |
| Genes AMR | blaCTX-M, tetA, tetB, sul2, aadA |
| MLST | ST11 (E. coli) |

---

**🔔 Esta es la Parte 1 del README. Continúa en la Parte 2 con el Flujo de Trabajo completo.**
