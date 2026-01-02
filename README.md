# 🧬 Pipeline de Vigilancia Genómica y Análisis de Resistencia Antimicrobiana en Bacterias

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Bioinformatics](https://img.shields.io/badge/Bioinformatics-Pipeline-blue.svg)]()
[![Status](https://img.shields.io/badge/Status-Production-green.svg)]()

Este repositorio documenta un flujo de trabajo bioinformático completo para el análisis de genomas bacterianos clínicos utilizando datos de secuenciación de nueva generación (NGS). El pipeline integra tres estrategias de ensamblaje complementarias: **Ensamblaje con Illumina**, **Ensamblaje con Nanopore** y **Ensamblaje Híbrido (Illumina + Nanopore)**, junto con detección exhaustiva de genes de resistencia a antimicrobianos (AMR) y análisis de variantes genómicas.

**🎯 Caso de Estudio**: *Klebsiella pneumoniae* URO5550422 con genoma multi-secuencia (1 cromosoma + 6 plásmidos)

---

## 📋 Tabla de Contenidos

- [⚠️ Antes de Comenzar](#️-antes-de-comenzar)
- [Características del Pipeline](#-características-del-pipeline)
- [Estructura del Proyecto](#-estructura-del-proyecto)
- [Requisitos del Sistema](#-requisitos-del-sistema)
- [Instalación y Configuración](#️-instalación-y-configuración)
- [Configuración del Proyecto](#-configuración-del-proyecto)
- [Dataset de Ejemplo](#-dataset-de-ejemplo)
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

### 🚀 Inicio Rápido

```bash
# 1. Clonar repositorio
git clone https://github.com/tu-usuario/Bacterial_Genomics_Project.git
cd Bacterial_Genomics_Project

# 2. Crear estructura y descargar referencia
bash setup_project_structure.sh

# 3. Configurar ambientes Conda (primera vez - ~45 minutos)
bash scripts/setup_environments.sh

# 4. Verificar instalación
bash scripts/verify_installation.sh

# 5. Enlazar datos de secuenciación
bash scripts/link_raw_data.sh /ruta/illumina /ruta/nanopore

# 6. Ejecutar pipeline completo
bash scripts/run_full_pipeline.sh URO5550422

# 7. Ver resultados
firefox 08_results/FINAL_REPORT.html
```

### 📊 ¿Qué Puedo Hacer con Este Pipeline?

✅ **Ensamblar genomas bacterianos** de alta calidad  
✅ **Identificar genes de resistencia** a antibióticos (AMR)  
✅ **Detectar variantes genómicas** (SNPs, INDELs)  
✅ **Anotar genes y funciones** biológicas  
✅ **Comparar diferentes estrategias** de ensamblaje  
✅ **Analizar cromosomas y plásmidos** por separado  
✅ **Tipificar cepas** (MLST, detección de plásmidos)  
✅ **Generar reportes automatizados** para vigilancia epidemiológica  

---

## 🎯 Características del Pipeline

### Tecnologías Soportadas
- **Illumina** (lecturas cortas, paired-end): Alta precisión, ideal para SNPs/INDELs
- **Oxford Nanopore** (lecturas largas): Ensamblajes contiguos, cierre de plásmidos
- **Híbrido** (Illumina + Nanopore): Combina precisión y continuidad

### Análisis Incluidos
- ✅ Control de calidad exhaustivo (raw y trimmed reads)
- ✅ Tres estrategias de ensamblaje independientes
- ✅ Mapeo contra genoma multi-secuencia (cromosoma + plásmidos)
- ✅ Análisis de cobertura por secuencia individual
- ✅ Detección de genes AMR con múltiples bases de datos
- ✅ Anotación funcional de genomas
- ✅ Evaluación de calidad de ensamblajes
- ✅ MLST typing y detección de plásmidos
- ✅ Identificación de factores de virulencia
- ✅ Visualización y reportes integrados

### Características Especiales para *K. pneumoniae*
- 🔬 Análisis separado de cromosoma y 6 plásmidos
- 🧬 Detección de genes AMR en elementos móviles
- 📊 Perfiles de resistencia específicos de la especie
- 🗺️ Mapeo optimizado para genomas multi-secuencia

---

## 📂 Estructura del Proyecto

```text
Bacterial_Genomics_Project/
├── 00_raw_data/                    # Datos crudos de secuenciación
│   ├── illumina/                   # Lecturas paired-end
│   │   ├── URO5550422_1.fastq.gz  # Forward reads
│   │   └── URO5550422_2.fastq.gz  # Reverse reads
│   ├── nanopore/                   # Lecturas largas ONT
│   │   └── URO5550422_1.fastq.gz  # Long reads (nota: mismo nombre, diferente tecnología)
│   └── sample_metadata.txt         # Metadata de la muestra
│
├── 01_reference/                   # Genoma de referencia K. pneumoniae
│   ├── GCF_000240185.1_ASM24018v2_genomic.fna  # Referencia completa
│   ├── reference.fasta             # Enlace simbólico
│   └── reference_sequences.txt     # Índice: 1 cromosoma + 6 plásmidos
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
│   │   ├── contigs.fasta
│   │   ├── scaffolds.fasta
│   │   └── assembly_graph.fastg
│   ├── 02_nanopore_only/           # Flye (solo Nanopore)
│   │   ├── assembly.fasta
│   │   ├── assembly_info.txt
│   │   └── assembly_graph.gfa
│   ├── 03_hybrid/                  # Unicycler (Illumina + Nanopore)
│   │   ├── assembly.fasta
│   │   └── assembly.gfa
│   └── 04_quast_evaluation/        # Evaluación comparativa QUAST
│       └── report.html
│
├── 04_mapping/                     # Mapeo y análisis de variantes
│   ├── 01_illumina/                # BWA + Samtools
│   │   ├── aligned_sorted.bam
│   │   ├── flagstat.txt
│   │   └── coverage.txt
│   ├── 02_nanopore/                # Minimap2 + Samtools
│   │   ├── aligned_sorted.bam
│   │   └── coverage.txt
│   ├── 03_variants/                # BCFtools variant calling
│   │   ├── illumina_variants.vcf
│   │   ├── nanopore_variants.vcf
│   │   └── consensus.fasta
│   └── 04_coverage_analysis/       # Cobertura por cromosoma/plásmidos
│       ├── Chromosome.bam          # Cobertura solo cromosoma
│       ├── Plasmid_pKPHS1.bam      # Cobertura plásmido 1
│       ├── [...más plásmidos...]
│       └── coverage_summary.txt    # Resumen por secuencia
│
├── 05_annotation/                  # Anotación funcional
│   ├── 01_prokka/                  # Anotación Prokka
│   │   ├── URO5550422.gff
│   │   ├── URO5550422.gbk
│   │   ├── URO5550422.faa
│   │   └── URO5550422.ffn
│   ├── 02_bakta/                   # Anotación Bakta (alternativa)
│   └── prokka_config.txt           # Configuración específica K. pneumoniae
│
├── 06_amr_screening/               # Detección de genes AMR
│   ├── amrfinder_db/               # Base de datos local AMRFinderPlus
│   │   └── latest/
│   ├── 01_amrfinder/               # Resultados AMRFinderPlus (NCBI)
│   │   ├── amrfinder_results.tsv
│   │   └── amrfinder_summary.txt
│   ├── 02_abricate/                # Resultados Abricate (múltiples DBs)
│   │   ├── card_results.tsv
│   │   ├── resfinder_results.tsv
│   │   ├── ncbi_results.tsv
│   │   └── abricate_summary.tsv
│   └── 03_rgi/                     # Resultados RGI/CARD
│       ├── rgi_results.txt
│       └── rgi_heatmap.png
│
├── 07_typing/                      # Tipificación molecular
│   ├── mlst/                       # MLST typing
│   │   └── mlst_results.txt
│   ├── plasmids/                   # Detección de plásmidos
│   │   ├── plasmidfinder_results.txt
│   │   └── plasmid_reconstruction/
│   └── virulence/                  # Factores de virulencia
│       └── vfdb_results.txt
│
├── 08_results/                     # Resultados consolidados y figuras
│   ├── figures/
│   │   ├── assembly_comparison.png
│   │   ├── coverage_plot.png
│   │   └── amr_heatmap.png
│   ├── tables/
│   │   ├── amr_summary.xlsx
│   │   └── variant_summary.xlsx
│   └── reports/
│       ├── quality_dashboard.html
│       └── FINAL_REPORT.html
│
├── envs/                           # Archivos YAML de ambientes Conda
│   ├── bact_main.yml
│   ├── bact_amr.yml
│   └── bact_rgi.yml
│
├── scripts/                        # Scripts de automatización
│   ├── setup_environments.sh       # Instalación de ambientes
│   ├── verify_installation.sh      # Verificación de instalación
│   ├── link_raw_data.sh            # Enlazar datos crudos
│   ├── run_full_pipeline.sh        # Pipeline completo
│   ├── 01_qc_illumina.sh
│   ├── 02_qc_nanopore.sh
│   ├── 03_assembly_illumina.sh
│   ├── 04_assembly_nanopore.sh
│   ├── 05_assembly_hybrid.sh
│   ├── 06_mapping.sh
│   ├── 07_annotation.sh
│   ├── 08_amr_screening.sh
│   ├── 09_typing.sh
│   └── utils/
│       ├── analyze_coverage_per_sequence.sh  # Análisis cromosoma/plásmidos
│       ├── calculate_metrics.sh
│       ├── compare_amr_tools.py
│       ├── generate_plots.py
│       └── extract_plasmids.sh
│
├── test_data/                      # Datos de prueba
│
├── logs/                           # Logs de ejecución
│   └── [timestamp]_pipeline.log
│
├── setup_project_structure.sh      # Script de configuración inicial
├── PROJECT_CONFIG.md               # Configuración del proyecto
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
- **Wget/Curl**: Para descargas
- **Conexión a internet**: Requerida para instalación inicial

### Tiempo de Ejecución Estimado

| Análisis | Hardware Mínimo | Hardware Recomendado |
|----------|-----------------|---------------------|
| QC Completo | 30-60 min | 15-30 min |
| Ensamblaje Illumina | 2-4 horas | 1-2 horas |
| Ensamblaje Nanopore | 1-2 horas | 30-60 min |
| Ensamblaje Híbrido | 4-8 horas | 2-4 horas |
| Mapeo + Variantes | 1-2 horas | 30-60 min |
| Detección AMR | 30-60 min | 15-30 min |
| Anotación | 30-60 min | 15-30 min |
| **Pipeline Completo** | **10-18 horas** | **5-9 horas** |

---

## 🛠️ Instalación y Configuración

### Paso 1: Clonar el Repositorio

```bash
# Clonar repositorio
git clone https://github.com/tu-usuario/Bacterial_Genomics_Project.git
cd Bacterial_Genomics_Project

# Verificar contenido
ls -lh
```

### Paso 2: Instalar Miniforge (Gestor de Paquetes)

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

### Paso 3: Configurar Canales de Bioconda

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

### Paso 4: Crear los Tres Ambientes Especializados

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

### Paso 5: Script de Instalación Automatizada (Recomendado)

En lugar de instalar manualmente cada ambiente, usa el script automatizado:

```bash
# Dar permisos de ejecución
chmod +x scripts/setup_environments.sh

# Ejecutar instalación automatizada
bash scripts/setup_environments.sh

# Tiempo total estimado: ~45 minutos
```

Este script:
- ✅ Configura los 3 ambientes automáticamente
- ✅ Descarga todas las bases de datos necesarias
- ✅ Verifica que todo esté correctamente instalado
- ✅ Muestra un resumen al finalizar

### Paso 6: Verificar Instalación

```bash
# Ejecutar script de verificación
bash scripts/verify_installation.sh

# Salida esperada:
# ========================================
# Verificación de Instalación
# ========================================
# 
# [Ambiente: bact_main]
# ✓ FastQC: OK
# ✓ fastp: OK
# ✓ BWA: OK
# ✓ Samtools: OK
# ✓ SPAdes: OK
# ✓ Flye: OK
# ✓ Unicycler: OK
# ✓ QUAST: OK
# ✓ AMRFinderPlus: OK
# ✓ MLST: OK
# 
# [Ambiente: bact_amr]
# ✓ Prokka: OK
# ✓ Abricate: OK
# 
# [Ambiente: bact_rgi]
# ✓ RGI: OK
# 
# [Bases de Datos]
# ✓ AMRFinderPlus DB: Instalada
# ✓ Abricate DBs: 8 bases disponibles
# ✓ CARD DB: Instalada
# 
# ========================================
# ✓ TODAS LAS VERIFICACIONES PASARON
# El sistema está listo para usar
# ========================================
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

**💡 Uso de ambientes exportados**:

```bash
# En otro servidor, recrear ambientes desde archivos YAML
mamba env create -f envs/bact_main.yml
mamba env create -f envs/bact_amr.yml
mamba env create -f envs/bact_rgi.yml

# Luego configurar bases de datos
bash scripts/setup_databases.sh
```

---

## 🔧 Configuración del Proyecto

### Paso 1: Crear Estructura y Descargar Referencia

```bash
# Ejecutar script de configuración inicial
bash setup_project_structure.sh

# Este script automáticamente:
# 1. Crea todos los directorios necesarios
# 2. Descarga el genoma de referencia K. pneumoniae (GCF_000240185.1)
# 3. Crea archivo de metadata
# 4. Genera documentación del proyecto
# 5. Crea scripts auxiliares
```

**Salida esperada**:

```
========================================
Configuración del Proyecto de Genómica Bacteriana
========================================

Muestra: URO5550422
Organismo: Klebsiella pneumoniae
Cepa de referencia: HS11286

[Paso 1/9] Creando estructura de directorios
✓ Creado: 00_raw_data/illumina
✓ Creado: 00_raw_data/nanopore
✓ Creado: 01_reference
[...más directorios...]

[Paso 2/9] Creando archivo de metadata
✓ Archivo de metadata creado

[Paso 3/9] Descargando genoma de referencia
ℹ Descargando GCF_000240185.1_ASM24018v2_genomic.fna.gz...
✓ Descarga completada
✓ Genoma de referencia listo

[Paso 4/9] Creando índice de secuencias de referencia
✓ Índice creado: 01_reference/reference_sequences.txt

[...más pasos...]

========================================
✓ Configuración Completada
========================================

Próximos pasos:
1. Enlazar datos de secuenciación
2. Ejecutar pipeline completo
```

### Paso 2: Revisar Información del Genoma de Referencia

```bash
# Ver información de las secuencias
cat 01_reference/reference_sequences.txt
```

**Contenido esperado**:

```
# Secuencias del Genoma de Referencia
# Klebsiella pneumoniae HS11286
# Accession: GCF_000240185.1

SeqID           Length      Type            Description
NC_016845.1     5333942     Chromosome      Cromosoma principal
NC_016838.1     122799      Plasmid         Plásmido pKPHS1
NC_016846.1     111195      Plasmid         Plásmido pKPHS2
NC_016839.1     105974      Plasmid         Plásmido pKPHS3
NC_016840.1     3751        Plasmid         Plásmido pKPHS4
NC_016847.1     3353        Plasmid         Plásmido pKPHS5
NC_016841.1     1308        Plasmid         Plásmido pKPHS6

# Total Genome Size: 5,682,322 bp
# Chromosome: 5,333,942 bp (93.9%)
# Plasmids: 348,380 bp (6.1%)
```

### Paso 3: Leer Configuración Completa del Proyecto

```bash
# Ver documentación completa
cat PROJECT_CONFIG.md

# O abrirlo con editor
nano PROJECT_CONFIG.md
```

Este archivo contiene:
- ✅ Información detallada de la muestra
- ✅ Descripción de las 7 secuencias (cromosoma + plásmidos)
- ✅ Consideraciones importantes para el análisis
- ✅ Comandos específicos para K. pneumoniae
- ✅ Referencias y próximos pasos

---

## 📊 Dataset de Ejemplo: URO5550422

### Información de la Muestra

- **ID**: URO5550422
- **Organismo**: *Klebsiella pneumoniae*
- **Origen**: Aislado clínico (urinario)
- **Referencia**: K. pneumoniae subsp. pneumoniae HS11286 (GCF_000240185.1)

### Archivos de Secuenciación

#### Illumina (Paired-end)
```
00_raw_data/illumina/
├── URO5550422_1.fastq.gz    # Forward reads (R1)
└── URO5550422_2.fastq.gz    # Reverse reads (R2)
```

**Especificaciones**:
- Plataforma: Illumina (MiSeq/NextSeq/NovaSeq)
- Química: Paired-end
- Longitud esperada: 150-300 bp
- Cobertura esperada: >50x

#### Nanopore (Long reads)
```
00_raw_data/nanopore/
└── URO5550422_1.fastq.gz    # Long reads
```

**⚠️ NOTA IMPORTANTE**: Este archivo tiene el mismo nombre que el R1 de Illumina, pero corresponde a **tecnología Nanopore**. Los archivos deben estar en directorios separados.

**Especificaciones**:
- Plataforma: Oxford Nanopore (MinION/GridION)
- Longitud esperada: 1-50 kb
- Cobertura esperada: >30x
- Calidad esperada: Q10-Q15

### Genoma de Referencia

**Archivo**: `GCF_000240185.1_ASM24018v2_genomic.fna`

**Composición genómica**:

| Secuencia | Accesión | Longitud (bp) | Tipo | % del Genoma |
|-----------|----------|---------------|------|--------------|
| Chromosome | NC_016845.1 | 5,333,942 | Cromosoma | 93.9% |
| pKPHS1 | NC_016838.1 | 122,799 | Plásmido | 2.2% |
| pKPHS2 | NC_016846.1 | 111,195 | Plásmido | 2.0% |
| pKPHS3 | NC_016839.1 | 105,974 | Plásmido | 1.9% |
| pKPHS4 | NC_016840.1 | 3,751 | Plásmido | 0.07% |
| pKPHS5 | NC_016847.1 | 3,353 | Plásmido | 0.06% |
| pKPHS6 | NC_016841.1 | 1,308 | Plásmido | 0.02% |
| **TOTAL** | - | **5,682,322** | - | **100%** |

### Enlazar Datos de Secuenciación

Una vez que tengas tus archivos de secuenciación, enlázalos al proyecto:

```bash
# Opción 1: Archivos en directorios separados (RECOMENDADO)
bash scripts/link_raw_data.sh /ruta/illumina /ruta/nanopore

# Opción 2: Archivos en el mismo directorio
# (El script los diferenciará por subdirectorio de destino)
bash scripts/link_raw_data.sh /ruta/datos /ruta/datos

# Verificar enlaces
ls -lh 00_raw_data/illumina/
ls -lh 00_raw_data/nanopore/

# Salida esperada:
# 00_raw_data/illumina/URO5550422_1.fastq.gz -> /ruta/real/URO5550422_1.fastq.gz
# 00_raw_data/illumina/URO5
