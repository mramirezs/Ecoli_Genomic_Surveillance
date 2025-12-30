# 🧬 Pipeline de Vigilancia Genómica y Análisis de Resistencia Antimicrobiana en Bacterias

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Bioinformatics](https://img.shields.io/badge/Bioinformatics-Pipeline-blue.svg)]()

Este repositorio documenta un flujo de trabajo bioinformático completo para el análisis de genomas bacterianos clínicos utilizando datos de secuenciación de nueva generación (NGS). El pipeline integra tres estrategias de análisis complementarias: **Ensamblaje con Illumina**, **Ensamblaje con Nanopore** y **Ensamblaje Híbrido (Illumina + Nanopore)**, junto con detección exhaustiva de genes de resistencia a antimicrobianos (AMR) y análisis de variantes genómicas.

---

## 📋 Tabla de Contenidos

- [Características del Pipeline](#-características-del-pipeline)
- [Estructura del Proyecto](#-estructura-del-proyecto)
- [Requisitos del Sistema](#-requisitos-del-sistema)
- [Instalación y Configuración](#-instalación-y-configuración)
- [Flujo de Trabajo](#-flujo-de-trabajo)
- [Resultados Esperados](#-resultados-esperados)
- [Interpretación de Resultados](#-interpretación-de-resultados)
- [Solución de Problemas](#-solución-de-problemas)
- [Referencias](#-referencias)

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
│   └── 03_variants/                # BCFtools variant calling
│       ├── illumina_variants.vcf
│       ├── nanopore_variants.vcf
│       └── consensus.fasta
│
├── 05_annotation/                  # Anotación funcional
│   ├── 01_prokka/                  # Anotación Prokka
│   │   ├── genome.gff
│   │   ├── genome.gbk
│   │   ├── genome.faa
│   │   └── genome.ffn
│   └── 02_bakta/                   # Anotación Bakta (alternativa)
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
├── 07_results/                     # Resultados consolidados y figuras
│   ├── assembly_comparison.png
│   ├── amr_summary.xlsx
│   └── final_report.html
│
├── envs/                           # Archivos YAML de ambientes Conda
│   ├── bact_main.yml
│   ├── bact_amr.yml
│   └── bact_rgi.yml
│
├── scripts/                        # Scripts de automatización
│   ├── 01_qc_illumina.sh
│   ├── 02_qc_nanopore.sh
│   ├── 03_assembly_illumina.sh
│   ├── 04_assembly_nanopore.sh
│   ├── 05_assembly_hybrid.sh
│   ├── 06_mapping.sh
│   ├── 07_annotation.sh
│   ├── 08_amr_screening.sh
│   └── run_full_pipeline.sh
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
- **CPU**: Mínimo 8 cores (16+ cores recomendado para ensamblaje híbrido)
- **RAM**: Mínimo 16 GB (32+ GB recomendado)
- **Almacenamiento**: 50-100 GB libres por muestra (dependiendo de la cobertura)

### Software Base
- Linux/Unix (Ubuntu 20.04+, CentOS 7+, o similar)
- Bash shell
- Git
- Conexión a internet (para instalación de herramientas)

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
conda create -n bact_main -c conda-forge -c bioconda -c defaults \
  python=3.10 pip pigz openjdk=11 -y

# Activar
conda activate bact_main

# Instalar herramientas de control de calidad
conda install fastqc multiqc fastp nanoplot filtlong -y

# Instalar herramientas de mapeo y análisis de variantes
conda install bwa minimap2 samtools bcftools bedtools blast -y

# Instalar ensambladores
conda install unicycler flye spades quast bandage -y

# Instalar herramientas AMR
conda install ncbi-amrfinderplus barrnap -y

# Configurar base de datos AMRFinderPlus (primera vez)
mkdir -p 06_amr_screening/amrfinder_db
amrfinder_update --database 06_amr_screening/amrfinder_db
```

**⏱️ Tiempo de instalación**: ~15 minutos  
**📦 Descarga de base de datos**: ~500 MB adicionales

#### 🦠 Ambiente 2: `bact_amr` (Anotación y AMR)

Dedicado a Prokka y Abricate, que requieren versiones específicas de Perl.

```bash
# Crear ambiente
mamba create -n bact_amr -c conda-forge -c bioconda -c defaults \
  python=3.9 prokka abricate -y

# Activar y configurar bases de datos
mamba activate bact_amr
abricate --setupdb
```

**⏱️ Tiempo de instalación**: ~10 minutos  
**📦 Descarga de bases de datos**: ~100 MB adicionales

#### 🧪 Ambiente 3: `bact_rgi` (AMR Avanzado)

Para RGI (Resistance Gene Identifier) con base de datos CARD.

```bash
# Crear ambiente
mamba create -n bact_rgi -c conda-forge -c bioconda -c defaults \
  python=3.11 rgi -y

# Activar
mamba activate bact_rgi

# Descargar y cargar base de datos CARD
mkdir -p 06_amr_screening/rgi
cd 06_amr_screening/rgi
wget https://card.mcmaster.ca/latest/data
tar -xvf data
rgi load --card_json card.json --local
cd ../..
```

**⏱️ Tiempo de instalación**: ~10 minutos  
**📦 Descarga de base de datos CARD**: ~50 MB

### Paso 4: Verificar Instalaciones

```bash
# Verificar bact_main
conda activate bact_main
fastqc --version
bwa 2>&1 | head -3
samtools --version
unicycler --version
spades.py --version
flye --version
quast --version
amrfinder --version

# Verificar bact_amr
conda activate bact_amr
prokka --version
abricate --version
abricate --list

# Verificar bact_rgi
conda activate bact_rgi
rgi main --version
rgi database --version --local
```

### Paso 5: Exportar Ambientes (Reproducibilidad)

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
```

### Paso 6: Clonar o Replicar en Otro Servidor

```bash
# Opción A: Clonar repositorio
git clone https://github.com/tu-usuario/Bacterial_Genomics_Project.git
cd Bacterial_Genomics_Project

# Opción B: Copiar archivos YML
scp envs/*.yml usuario@servidor:/ruta/proyecto/envs/

# Crear ambientes desde YML
mamba env create -f envs/bact_main.yml
mamba env create -f envs/bact_amr.yml
mamba env create -f envs/bact_rgi.yml

# Configurar bases de datos
conda activate bact_main
amrfinder_update --database 06_amr_screening/amrfinder_db

conda activate bact_amr
abricate --setupdb

conda activate bact_rgi
# Descargar CARD y ejecutar: rgi load --card_json card.json --local
```

---

## 🔬 Flujo de Trabajo

### Fase 1: Preparación de Datos

#### 1.1 Crear Enlaces Simbólicos a Datos Crudos

```bash
# Crear directorio de datos crudos
mkdir -p 00_raw_data/illumina 00_raw_data/nanopore

# Crear enlaces simbólicos (evita duplicar datos)
ln -s /ruta/absoluta/datos/sample_R1.fastq.gz 00_raw_data/illumina/
ln -s /ruta/absoluta/datos/sample_R2.fastq.gz 00_raw_data/illumina/
ln -s /ruta/absoluta/datos/sample_ont.fastq.gz 00_raw_data/nanopore/
```

#### 1.2 Descargar Genoma de Referencia (Opcional)

Para análisis de mapeo y detección de variantes:

```bash
mkdir -p 01_reference

# Ejemplo: Descargar E. coli K-12 MG1655 desde NCBI
# Para otras bacterias, buscar en NCBI Genome: https://www.ncbi.nlm.nih.gov/genome/
wget -O 01_reference/reference.fasta.gz \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"

gunzip 01_reference/reference.fasta.gz
```

---

### Fase 2: Control de Calidad (QC)

#### 2.1 QC de Lecturas Illumina

```bash
conda activate bact_main

# Crear directorios
mkdir -p 02_qc/01_illumina_raw 02_qc/02_illumina_trimmed

# FastQC en datos crudos
fastqc 00_raw_data/illumina/*.fastq.gz \
  -o 02_qc/01_illumina_raw/ \
  -t 8

# Limpieza y recorte con fastp
fastp \
  -i 00_raw_data/illumina/sample_R1.fastq.gz \
  -I 00_raw_data/illumina/sample_R2.fastq.gz \
  -o 02_qc/02_illumina_trimmed/sample_R1_trimmed.fastq.gz \
  -O 02_qc/02_illumina_trimmed/sample_R2_trimmed.fastq.gz \
  --detect_adapter_for_pe \
  --cut_front --cut_tail \
  --trim_poly_g \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 40 \
  --n_base_limit 5 \
  --length_required 50 \
  --thread 8 \
  --html 02_qc/02_illumina_trimmed/fastp_report.html \
  --json 02_qc/02_illumina_trimmed/fastp_report.json

# FastQC en datos limpios
fastqc 02_qc/02_illumina_trimmed/*_trimmed.fastq.gz \
  -o 02_qc/02_illumina_trimmed/ \
  -t 8
```

**📊 Resultados QC Illumina**

_[Incluir aquí capturas de pantalla o estadísticas clave]_

| Métrica | Raw Reads | Trimmed Reads |
|---------|-----------|---------------|
| Total Reads | | |
| % Bases ≥Q30 | | |
| GC Content (%) | | |
| Duplicación (%) | | |
| Adaptadores Detectados | | |

---

#### 2.2 QC de Lecturas Nanopore

```bash
conda activate bact_main

# Crear directorios
mkdir -p 02_qc/03_nanopore_raw 02_qc/04_nanopore_filtered

# NanoPlot en datos crudos
NanoPlot \
  --fastq 00_raw_data/nanopore/sample_ont.fastq.gz \
  -o 02_qc/03_nanopore_raw/ \
  -t 8 \
  --plots kde

# Filtrado con Filtlong
filtlong \
  --min_length 1000 \
  --keep_percent 90 \
  --target_bases 500000000 \
  00_raw_data/nanopore/sample_ont.fastq.gz | \
  pigz > 02_qc/04_nanopore_filtered/sample_ont_filtered.fastq.gz

# NanoPlot en datos filtrados
NanoPlot \
  --fastq 02_qc/04_nanopore_filtered/sample_ont_filtered.fastq.gz \
  -o 02_qc/04_nanopore_filtered/ \
  -t 8 \
  --plots kde
```

**📊 Resultados QC Nanopore**

_[Incluir aquí gráficos de distribución de longitud y calidad]_

| Métrica | Raw Reads | Filtered Reads |
|---------|-----------|----------------|
| Total Reads | | |
| Mean Read Length (bp) | | |
| Median Read Length (bp) | | |
| Mean Quality Score | | |
| N50 (bp) | | |
| Total Bases (Gb) | | |

---

#### 2.3 Reporte Consolidado con MultiQC

```bash
conda activate bact_main

mkdir -p 02_qc/05_multiqc

# Generar reporte integrado
multiqc 02_qc/ \
  -o 02_qc/05_multiqc/ \
  --filename multiqc_report_complete
```

**📊 Reporte MultiQC**

_[Enlace a reporte HTML o capturas de pantalla clave]_

---

### Fase 3: Estrategias de Ensamblaje

#### 3.1 Ensamblaje Solo Illumina (SPAdes)

```bash
conda activate bact_main

mkdir -p 03_assembly/01_illumina_only

# Ensamblaje con SPAdes
spades.py \
  -1 02_qc/02_illumina_trimmed/sample_R1_trimmed.fastq.gz \
  -2 02_qc/02_illumina_trimmed/sample_R2_trimmed.fastq.gz \
  -o 03_assembly/01_illumina_only/ \
  --careful \
  -t 8 -m 16

# Copiar contigs finales
cp 03_assembly/01_illumina_only/contigs.fasta \
   03_assembly/01_illumina_only/assembly_illumina.fasta
```

**📊 Estadísticas Ensamblaje Illumina**

| Métrica | Valor |
|---------|-------|
| Número de Contigs | |
| Tamaño Total del Ensamblaje (bp) | |
| Contig Más Largo (bp) | |
| N50 (bp) | |
| L50 | |
| GC Content (%) | |

---

#### 3.2 Ensamblaje Solo Nanopore (Flye)

```bash
conda activate bact_main

mkdir -p 03_assembly/02_nanopore_only

# Ensamblaje con Flye
flye \
  --nano-raw 02_qc/04_nanopore_filtered/sample_ont_filtered.fastq.gz \
  --out-dir 03_assembly/02_nanopore_only/ \
  --threads 8 \
  --genome-size 5m

# Copiar ensamblaje final
cp 03_assembly/02_nanopore_only/assembly.fasta \
   03_assembly/02_nanopore_only/assembly_nanopore.fasta
```

**📊 Estadísticas Ensamblaje Nanopore**

| Métrica | Valor |
|---------|-------|
| Número de Contigs | |
| Tamaño Total del Ensamblaje (bp) | |
| Contig Más Largo (bp) | |
| N50 (bp) | |
| L50 | |
| GC Content (%) | |
| Circularidad Detectada | |

---

#### 3.3 Ensamblaje Híbrido (Unicycler)

```bash
conda activate bact_main

mkdir -p 03_assembly/03_hybrid

# Ensamblaje híbrido con Unicycler
unicycler \
  -1 02_qc/02_illumina_trimmed/sample_R1_trimmed.fastq.gz \
  -2 02_qc/02_illumina_trimmed/sample_R2_trimmed.fastq.gz \
  -l 02_qc/04_nanopore_filtered/sample_ont_filtered.fastq.gz \
  -o 03_assembly/03_hybrid/ \
  -t 8

# Copiar ensamblaje final
cp 03_assembly/03_hybrid/assembly.fasta \
   03_assembly/03_hybrid/assembly_hybrid.fasta
```

**📊 Estadísticas Ensamblaje Híbrido**

| Métrica | Valor |
|---------|-------|
| Número de Contigs | |
| Tamaño Total del Ensamblaje (bp) | |
| Contig Más Largo (bp) | |
| N50 (bp) | |
| L50 | |
| GC Content (%) | |
| Circularidad Detectada | |

---

#### 3.4 Evaluación Comparativa de Ensamblajes (QUAST)

```bash
conda activate bact_main

mkdir -p 03_assembly/04_quast_evaluation

# Evaluación con QUAST (con referencia)
quast.py \
  03_assembly/01_illumina_only/assembly_illumina.fasta \
  03_assembly/02_nanopore_only/assembly_nanopore.fasta \
  03_assembly/03_hybrid/assembly_hybrid.fasta \
  -r 01_reference/reference.fasta \
  -o 03_assembly/04_quast_evaluation/ \
  --threads 8 \
  --labels "Illumina,Nanopore,Hybrid"

# Si no tienes referencia, omite el parámetro -r
```

**📊 Comparación de Ensamblajes (QUAST)**

_[Incluir tabla comparativa generada por QUAST]_

| Métrica | Illumina | Nanopore | Híbrido |
|---------|----------|----------|---------|
| Contigs (≥500 bp) | | | |
| Tamaño Total (bp) | | | |
| Contig Más Largo (bp) | | | |
| N50 (bp) | | | |
| L50 | | | |
| GC (%) | | | |
| Genes Predichos | | | |
| % Genoma Cubierto | | | |
| Mismatches por 100 kb | | | |

**🎯 Recomendación de Ensamblaje:**

_[Seleccionar el mejor ensamblaje basado en métricas QUAST]_

---

### Fase 4: Mapeo y Análisis de Variantes

#### 4.1 Mapeo de Lecturas Illumina

```bash
conda activate bact_main

mkdir -p 04_mapping/01_illumina

# Indexar referencia (solo primera vez)
bwa index 01_reference/reference.fasta

# Mapeo con BWA-MEM
bwa mem -t 8 \
  01_reference/reference.fasta \
  02_qc/02_illumina_trimmed/sample_R1_trimmed.fastq.gz \
  02_qc/02_illumina_trimmed/sample_R2_trimmed.fastq.gz | \
  samtools view -Sb - | \
  samtools sort -@ 8 -o 04_mapping/01_illumina/aligned_sorted.bam

# Indexar BAM
samtools index 04_mapping/01_illumina/aligned_sorted.bam

# Estadísticas de mapeo
samtools flagstat 04_mapping/01_illumina/aligned_sorted.bam > \
  04_mapping/01_illumina/flagstat.txt

samtools coverage 04_mapping/01_illumina/aligned_sorted.bam > \
  04_mapping/01_illumina/coverage.txt

samtools depth 04_mapping/01_illumina/aligned_sorted.bam | \
  awk '{sum+=$3} END {print "Mean Depth:", sum/NR}' > \
  04_mapping/01_illumina/mean_depth.txt
```

**📊 Estadísticas de Mapeo Illumina**

| Métrica | Valor |
|---------|-------|
| Total Reads | |
| Reads Mapeadas (%) | |
| Reads Paired (%) | |
| Cobertura Media | |
| Duplicados (%) | |

---

#### 4.2 Mapeo de Lecturas Nanopore

```bash
conda activate bact_main

mkdir -p 04_mapping/02_nanopore

# Mapeo con Minimap2
minimap2 -ax map-ont -t 8 \
  01_reference/reference.fasta \
  02_qc/04_nanopore_filtered/sample_ont_filtered.fastq.gz | \
  samtools view -Sb - | \
  samtools sort -@ 8 -o 04_mapping/02_nanopore/aligned_sorted.bam

# Indexar BAM
samtools index 04_mapping/02_nanopore/aligned_sorted.bam

# Estadísticas
samtools flagstat 04_mapping/02_nanopore/aligned_sorted.bam > \
  04_mapping/02_nanopore/flagstat.txt

samtools coverage 04_mapping/02_nanopore/aligned_sorted.bam > \
  04_mapping/02_nanopore/coverage.txt
```

**📊 Estadísticas de Mapeo Nanopore**

| Métrica | Valor |
|---------|-------|
| Total Reads | |
| Reads Mapeadas (%) | |
| Cobertura Media | |

---

#### 4.3 Llamado de Variantes y Consenso

```bash
conda activate bact_main

mkdir -p 04_mapping/03_variants

# Llamado de variantes Illumina
bcftools mpileup -Ou -f 01_reference/reference.fasta \
  04_mapping/01_illumina/aligned_sorted.bam | \
  bcftools call -mv -Oz -o 04_mapping/03_variants/illumina_variants.vcf.gz

bcftools index 04_mapping/03_variants/illumina_variants.vcf.gz

# Llamado de variantes Nanopore
bcftools mpileup -Ou -f 01_reference/reference.fasta \
  04_mapping/02_nanopore/aligned_sorted.bam | \
  bcftools call -mv -Oz -o 04_mapping/03_variants/nanopore_variants.vcf.gz

bcftools index 04_mapping/03_variants/nanopore_variants.vcf.gz

# Generar secuencia consenso (Illumina)
bcftools consensus -f 01_reference/reference.fasta \
  04_mapping/03_variants/illumina_variants.vcf.gz > \
  04_mapping/03_variants/consensus_illumina.fasta

# Estadísticas de variantes
bcftools stats 04_mapping/03_variants/illumina_variants.vcf.gz > \
  04_mapping/03_variants/illumina_variants_stats.txt

bcftools stats 04_mapping/03_variants/nanopore_variants.vcf.gz > \
  04_mapping/03_variants/nanopore_variants_stats.txt
```

**📊 Variantes Detectadas**

| Tipo de Variante | Illumina | Nanopore |
|------------------|----------|----------|
| SNPs | | |
| INDELs | | |
| Variantes en Genes | | |

---

### Fase 5: Anotación Funcional

#### 5.1 Anotación con Prokka

```bash
conda activate bact_amr

mkdir -p 05_annotation/01_prokka

# Anotar el mejor ensamblaje (en este ejemplo, el híbrido)
prokka \
  --outdir 05_annotation/01_prokka/ \
  --prefix sample_genome \
  --kingdom Bacteria \
  --genus [Género] \
  --species [especie] \
  --strain [cepa] \
  --gram [neg/pos] \
  --usegenus \
  --addgenes \
  --addmrna \
  --rfam \
  --cpus 8 \
  03_assembly/03_hybrid/assembly_hybrid.fasta
```

**💡 Nota**: Ajusta los parámetros `--genus`, `--species`, `--strain` y `--gram` según tu bacteria de estudio.

**📊 Estadísticas de Anotación**

| Característica | Cantidad |
|----------------|----------|
| Secuencias Anotadas | |
| Genes (CDS) | |
| rRNA | |
| tRNA | |
| tmRNA | |
| CRISPR arrays | |
| Tamaño Total (bp) | |
| GC Content (%) | |

**🗂️ Archivos Generados**:
- `sample_genome.gff`: Anotaciones en formato GFF3
- `sample_genome.gbk`: GenBank format
- `sample_genome.faa`: Secuencias proteicas
- `sample_genome.ffn`: Secuencias de genes
- `sample_genome.txt`: Resumen de anotación

---

#### 5.2 Anotación con Bakta (Alternativa Moderna)

```bash
conda activate bact_main

# Instalar Bakta (si no está instalado)
# conda install bakta -y

mkdir -p 05_annotation/02_bakta

# Descargar base de datos Bakta (primera vez, ~30 GB)
# bakta_db download --output 05_annotation/bakta_db

# Anotar genoma
bakta \
  --db 05_annotation/bakta_db \
  --output 05_annotation/02_bakta/ \
  --prefix sample_genome \
  --locus-tag SAMPLE \
  --threads 8 \
  03_assembly/03_hybrid/assembly_hybrid.fasta
```

**💡 Ventaja de Bakta**: Anotación más actualizada y rápida que Prokka, con mejor integración de bases de datos modernas.

---

### Fase 6: Detección de Genes de Resistencia Antimicrobiana (AMR)

#### 6.1 AMRFinderPlus (NCBI - Recomendado)

```bash
conda activate bact_main

mkdir -p 06_amr_screening/01_amrfinder

# Verificar base de datos actualizada
amrfinder --database 06_amr_screening/amrfinder_db --list_organisms

# Ejecutar AMRFinderPlus en el ensamblaje híbrido
amrfinder \
  --nucleotide 03_assembly/03_hybrid/assembly_hybrid.fasta \
  --database 06_amr_screening/amrfinder_db \
  --organism [Género] \
  --output 06_amr_screening/01_amrfinder/amrfinder_results.tsv \
  --plus \
  --name sample_hybrid \
  --threads 8

# Si tienes archivo de proteínas de Prokka
amrfinder \
  --protein 05_annotation/01_prokka/sample_genome.faa \
  --database 06_amr_screening/amrfinder_db \
  --organism [Género] \
  --output 06_amr_screening/01_amrfinder/amrfinder_protein_results.tsv \
  --plus \
  --threads 8

# Generar resumen
grep -v "^#" 06_amr_screening/01_amrfinder/amrfinder_results.tsv | \
  cut -f5,6,7,9,11,12 | \
  sort -u > 06_amr_screening/01_amrfinder/amrfinder_summary.txt
```

**📊 Genes AMR Detectados (AMRFinderPlus)**

_[Tabla resumen de genes encontrados]_

| Gen | Clase de Antibiótico | % Identity | % Coverage | Método |
|-----|---------------------|------------|------------|--------|
| | | | | |
| | | | | |
| | | | | |

**🦠 Perfil de Resistencia:**

_[Describir fenotipos de resistencia esperados basados en genes detectados]_

---

#### 6.2 Abricate (Múltiples Bases de Datos)

```bash
conda activate bact_amr

mkdir -p 06_amr_screening/02_abricate

# Ejecutar contra múltiples bases de datos
# CARD database
abricate --db card \
  03_assembly/03_hybrid/assembly_hybrid.fasta > \
  06_amr_screening/02_abricate/card_results.tsv

# ResFinder database
abricate --db resfinder \
  03_assembly/03_hybrid/assembly_hybrid.fasta > \
  06_amr_screening/02_abricate/resfinder_results.tsv

# NCBI database
abricate --db ncbi \
  03_assembly/03_hybrid/assembly_hybrid.fasta > \
  06_amr_screening/02_abricate/ncbi_results.tsv

# ARG-ANNOT database
abricate --db argannot \
  03_assembly/03_hybrid/assembly_hybrid.fasta > \
  06_amr_screening/02_abricate/argannot_results.tsv

# MEGARes database
abricate --db megares \
  03_assembly/03_hybrid/assembly_hybrid.fasta > \
  06_amr_screening/02_abricate/megares_results.tsv

# Resumen consolidado
abricate --summary 06_amr_screening/02_abricate/*.tsv > \
  06_amr_screening/02_abricate/abricate_summary.tsv
```

**📊 Comparación entre Bases de Datos**

| Base de Datos | Genes Detectados | Cobertura Promedio (%) | Identidad Promedio (%) |
|---------------|------------------|------------------------|------------------------|
| CARD | | | |
| ResFinder | | | |
| NCBI | | | |
| ARG-ANNOT | | | |
| MEGARes | | | |

---

#### 6.3 RGI (CARD - Análisis Avanzado)

```bash
conda activate bact_rgi

mkdir -p 06_amr_screening/03_rgi

# Verificar base de datos cargada
rgi database --version --local

# Ejecutar análisis RGI
rgi main \
  --input_sequence 03_assembly/03_hybrid/assembly_hybrid.fasta \
  --output_file 06_amr_screening/03_rgi/rgi_results \
  --input_type contig \
  --local \
  --clean \
  --num_threads 8

# Generar heatmap
rgi heatmap \
  --input 06_amr_screening/03_rgi/rgi_results.txt \
  --output 06_amr_screening/03_rgi/rgi_heatmap

# Análisis de BWT (lectura de variantes en tiempo real, opcional)
# rgi bwt --help
```

**📊 Análisis RGI/CARD**

_[Incluir heatmap generado y tabla de genes]_

| Gen | ARO Accession | Mecanismo de Resistencia | Drug Class | % Identity |
|-----|---------------|--------------------------|------------|------------|
| | | | | |
| | | | | |

---

### Fase 7: Análisis Comparativo y Consolidación

#### 7.1 Comparación de Resultados AMR

```bash
# Crear directorio de resultados consolidados
mkdir -p 07_results

# Script Python para consolidar resultados (ejemplo básico)
cat > 07_results/consolidate_amr.py << 'EOF'
import pandas as pd
import sys

# Leer resultados AMRFinderPlus
amrf = pd.read_csv('06_amr_screening/01_amrfinder/amrfinder_results.tsv', sep='\t', comment='#')
amrf_genes = set(amrf['Gene symbol'].dropna())

# Leer resultados Abricate CARD
abr_card = pd.read_csv('06_amr_screening/02_abricate/card_results.tsv', sep='\t')
abr_genes = set(abr_card['GENE'].dropna())

# Leer resultados RGI
rgi = pd.read_csv('06_amr_screening/03_rgi/rgi_results.txt', sep='\t')
rgi_genes = set(rgi['Best_Hit_ARO'].dropna())

# Genes comunes
common_genes = amrf_genes & abr_genes & rgi_genes
print(f"Genes AMR detectados por las 3 herramientas: {len(common_genes)}")
print(common_genes)

# Genes únicos por herramienta
print(f"\nGenes únicos AMRFinderPlus: {amrf_genes - abr_genes - rgi_genes}")
print(f"Genes únicos Abricate: {abr_genes - amrf_genes - rgi_genes}")
print(f"Genes únicos RGI: {rgi_genes - amrf_genes - abr_genes}")
EOF

python 07_results/consolidate_amr.py > 07_results/amr_comparison.txt
```

**📊 Consolidación de Resultados AMR**

_[Diagrama de Venn o tabla comparativa]_

| Herramienta | Genes Detectados | Genes Únicos | Consenso con Otras |
|-------------|------------------|--------------|-------------------|
| AMRFinderPlus | | | |
| Abricate (CARD) | | | |
| RGI | | | |
| **Consenso (3 herramientas)** | | | |

**🎯 Genes AMR de Alta Confianza** (detectados por ≥2 herramientas):

_[Listar genes con descripción]_

---

#### 7.2 Visualización de Ensamblajes

```bash
conda activate bact_main

# Visualizar gráficos de ensamblaje con Bandage
# Para cada ensamblaje GFA disponible

# Ensamblaje Nanopore (Flye)
Bandage image 03_assembly/02_nanopore_only/assembly_graph.gfa \
  07_results/assembly_nanopore_graph.png \
  --height 2000 --width 2000

# Ensamblaje Híbrido (Unicycler)
Bandage image 03_assembly/03_hybrid/assembly.gfa \
  07_results/assembly_hybrid_graph.png \
  --height 2000 --width 2000
```

**📊 Gráficos de Ensamblaje**

_[Incluir imágenes generadas por Bandage]_

---

#### 7.3 Reporte Final Integrado

```bash
# Crear reporte Markdown consolidado
cat > 07_results/FINAL_REPORT.md << 'EOF'
# Reporte de Análisis Genómico Bacteriano

## Información de la Muestra
- **ID Muestra**: [sample_name]
- **Bacteria**: [Género especie]
- **Fecha de Análisis**: [fecha]
- **Tecnologías de Secuenciación**: Illumina + Nanopore

## Resumen Ejecutivo

### Control de Calidad
- **Illumina**: [X] millones de lecturas paired-end, Q30 > [Y]%
- **Nanopore**: [Z] Mb de lecturas largas, N50 = [N] kb

### Ensamblaje Seleccionado
- **Estrategia**: [Illumina/Nanopore/Híbrido]
- **Tamaño**: [X.X] Mb
- **Número de Contigs**: [N]
- **N50**: [X] kb
- **Completitud**: [X]% del genoma de referencia

### Genes AMR Detectados
- **Total de Genes AMR**: [N]
- **Clases de Antibióticos Afectadas**: [N]
- **Fenotipos de Resistencia Predichos**: [lista]

### Hallazgos Clave
- [Hallazgo 1]
- [Hallazgo 2]
- [Hallazgo 3]

## Detalles por Sección
[Completar con secciones previas]

## Conclusiones y Recomendaciones
[Escribir conclusiones]

EOF
```

---

## 📊 Resultados Esperados

### Archivos Principales Generados

| Archivo | Descripción | Ubicación |
|---------|-------------|-----------|
| `multiqc_report_complete.html` | Reporte consolidado de QC | `02_qc/05_multiqc/` |
| `assembly_[method].fasta` | Ensamblajes finales | `03_assembly/01-03_*/` |
| `quast_report.html` | Comparación de ensamblajes | `03_assembly/04_quast_evaluation/` |
| `aligned_sorted.bam` | Archivos de mapeo | `04_mapping/01-02_*/` |
| `*_variants.vcf.gz` | Variantes llamadas | `04_mapping/03_variants/` |
| `sample_genome.gff` | Anotación funcional | `05_annotation/01_prokka/` |
| `amrfinder_results.tsv` | Genes AMR (NCBI) | `06_amr_screening/01_amrfinder/` |
| `card_results.tsv` | Genes AMR (CARD) | `06_amr_screening/02_abricate/` |
| `rgi_results.txt` | Análisis AMR avanzado | `06_amr_screening/03_rgi/` |
| `FINAL_REPORT.md` | Reporte consolidado | `07_results/` |

---

## 🔍 Interpretación de Resultados

### 1. Evaluación de Calidad del Ensamblaje

**Métricas Clave**:
- **N50**: Cuanto mayor, mejor. Valores >50 kb son excelentes para genomas bacterianos.
- **Número de Contigs**: Menos es mejor. Idealmente <50 para un genoma bacteriano de ~5 Mb.
- **Cobertura**: Mínimo 30x para Illumina, 20x para Nanopore.
- **Completitud**: Comparar con genoma de referencia (QUAST).

**Interpretación**:
- **Ensamblaje Excelente**: N50 >100 kb, <20 contigs, cobertura >50x
- **Ensamblaje Bueno**: N50 50-100 kb, 20-50 contigs, cobertura 30-50x
- **Ensamblaje Aceptable**: N50 20-50 kb, 50-100 contigs, cobertura >20x
- **Requiere Mejora**: N50 <20 kb, >100 contigs, cobertura <20x

---

### 2. Interpretación de Genes AMR

#### Categorías de Resistencia

Los genes AMR detectados confieren resistencia a diferentes clases de antibióticos:

| Clase de Antibiótico | Genes Comunes | Impacto Clínico |
|---------------------|---------------|-----------------|
| **Beta-lactámicos** | blaCTX-M, blaTEM, blaOXA | Alto (tratamiento de primera línea) |
| **Aminoglicósidos** | aac, aph, ant | Moderado-Alto |
| **Quinolonas** | qnr, aac(6')-Ib-cr | Alto |
| **Tetraciclinas** | tet(A), tet(B) | Moderado |
| **Sulfonamidas** | sul1, sul2 | Moderado |
| **Trimetoprim** | dfrA | Moderado |
| **Fenicoles** | catA, cmlA | Moderado |
| **Macrólidos** | erm, mef | Moderado |

#### Niveles de Confianza

- **Alta Confianza**: Gen detectado por ≥2 herramientas, identidad >95%, cobertura >95%
- **Confianza Media**: Gen detectado por 1-2 herramientas, identidad 90-95%
- **Confianza Baja**: Gen detectado por 1 herramienta, identidad <90% o fragmentario

---

### 3. Variantes Clínicamente Relevantes

**Tipos de Variantes**:
- **SNPs en genes AMR**: Pueden conferir resistencia (ej. mutaciones en gyrA/parC para quinolonas)
- **INDELs**: Pueden inactivar genes (resistencia o virulencia)
- **Variantes estructurales**: Duplicaciones, deleciones grandes (detectadas mejor con Nanopore)

**Validación**:
- Comparar variantes entre Illumina y Nanopore (mayor confianza si concuerdan)
- Verificar cobertura de lecturas en posición de variante (mínimo 10x)
- Confirmar variantes AMR con literatura científica

---

### 4. Selección del Mejor Ensamblaje

**Criterios de Decisión**:

1. **Para SNPs/pequeñas variantes**: Preferir ensamblaje Illumina o Híbrido (mayor precisión)
2. **Para genes completos y plásmidos**: Preferir Híbrido (mejor continuidad)
3. **Para detección de elementos móviles**: Preferir Nanopore o Híbrido (lecturas largas)
4. **Para análisis filogenético**: Preferir Híbrido (equilibrio precisión/continuidad)

**Recomendación General**: Use el **ensamblaje híbrido** como principal para reportes, y valide hallazgos críticos con los otros ensamblajes.

---

## 🔧 Solución de Problemas

### Problema 1: Error "Could not solve for environment specs"

**Causa**: Conflictos de dependencias entre herramientas.

**Solución**:
```bash
# NO mezclar las herramientas de diferentes ambientes
# Usar los 3 ambientes separados como se describe
# Verificar canales configurados correctamente:
conda config --show channels
```

---

### Problema 2: Ensamblaje muy fragmentado

**Posibles Causas y Soluciones**:

| Causa | Solución |
|-------|----------|
| Baja cobertura | Incrementar profundidad de secuenciación (>30x) |
| Baja calidad de reads | Mejorar filtrado en QC (fastp/filtlong) |
| Genoma complejo (repeticiones) | Usar ensamblaje híbrido o solo Nanopore |
| Parámetros incorrectos | Ajustar `--careful` en SPAdes, `--genome-size` en Flye |

---

### Problema 3: Genes AMR no detectados

**Verificar**:
```bash
# 1. Base de datos actualizada
conda activate bact_main
amrfinder --database 06_amr_screening/amrfinder_db --database_version

# 2. Organismo correcto
amrfinder --list_organisms

# 3. Parámetros de identidad/cobertura
# Reducir umbral si sospechas genes divergentes:
amrfinder --nucleotide [input] --ident_min 80 --coverage_min 70
```

---

### Problema 4: Memoria insuficiente

**Soluciones**:
```bash
# 1. Limitar memoria en SPAdes
spades.py -m 16  # Usar máximo 16 GB

# 2. Reducir threads
--threads 4  # En lugar de 8

# 3. Subsamplear reads si cobertura >100x
seqtk sample -s100 input.fastq.gz 0.5 > subsampled.fastq.gz
```

---

### Problema 5: Base de datos CARD/RGI desactualizada

```bash
conda activate bact_rgi

# Descargar última versión
cd 06_amr_screening/rgi
wget -O card_data.tar.bz2 https://card.mcmaster.ca/latest/data
tar -xvf card_data.tar.bz2

# Recargar base de datos
rgi load --card_json card.json --local

# Verificar versión
rgi database --version --local
```

---

## 🚀 Automatización Completa

### Script Maestro

Crea un script para ejecutar todo el pipeline automáticamente:

```bash
#!/bin/bash
# run_full_pipeline.sh

set -euo pipefail  # Salir si hay errores

# Variables
SAMPLE="sample_name"
THREADS=8
MEMORY=16
ILLUMINA_R1="00_raw_data/illumina/${SAMPLE}_R1.fastq.gz"
ILLUMINA_R2="00_raw_data/illumina/${SAMPLE}_R2.fastq.gz"
NANOPORE="00_raw_data/nanopore/${SAMPLE}_ont.fastq.gz"
REFERENCE="01_reference/reference.fasta"

echo "========================================="
echo "Pipeline de Análisis Genómico Bacteriano"
echo "Muestra: ${SAMPLE}"
echo "Inicio: $(date)"
echo "========================================="

# Fase 1: QC Illumina
echo "[$(date)] Fase 1: Control de Calidad Illumina"
conda activate bact_main
bash scripts/01_qc_illumina.sh

# Fase 2: QC Nanopore
echo "[$(date)] Fase 2: Control de Calidad Nanopore"
bash scripts/02_qc_nanopore.sh

# Fase 3: Ensamblajes
echo "[$(date)] Fase 3: Ensamblaje Illumina"
bash scripts/03_assembly_illumina.sh

echo "[$(date)] Fase 4: Ensamblaje Nanopore"
bash scripts/04_assembly_nanopore.sh

echo "[$(date)] Fase 5: Ensamblaje Híbrido"
bash scripts/05_assembly_hybrid.sh

# Fase 4: Evaluación
echo "[$(date)] Fase 6: Evaluación de Ensamblajes"
bash scripts/06_quast_evaluation.sh

# Fase 5: Mapeo
echo "[$(date)] Fase 7: Mapeo y Variantes"
bash scripts/07_mapping.sh

# Fase 6: Anotación
echo "[$(date)] Fase 8: Anotación con Prokka"
conda activate bact_amr
bash scripts/08_annotation.sh

# Fase 7: AMR
echo "[$(date)] Fase 9: Detección de Genes AMR"
conda activate bact_main
bash scripts/09_amrfinder.sh

conda activate bact_amr
bash scripts/10_abricate.sh

conda activate bact_rgi
bash scripts/11_rgi.sh

# Fase 8: Consolidación
echo "[$(date)] Fase 10: Generación de Reportes"
conda activate bact_main
bash scripts/12_generate_reports.sh

echo "========================================="
echo "Pipeline Completado Exitosamente"
echo "Fin: $(date)"
echo "========================================="
```

---

## 📚 Referencias

### Herramientas Bioinformáticas

- **FastQC**: Andrews S. (2010). FastQC. [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- **fastp**: Chen et al. (2018). Bioinformatics. DOI: 10.1093/bioinformatics/bty560
- **NanoPlot**: De Coster et al. (2018). Bioinformatics. DOI: 10.1093/bioinformatics/bty149
- **BWA**: Li & Durbin (2009). Bioinformatics. DOI: 10.1093/bioinformatics/btp324
- **Minimap2**: Li (2018). Bioinformatics. DOI: 10.1093/bioinformatics/bty191
- **SPAdes**: Bankevich et al. (2012). J Comput Biol. DOI: 10.1089/cmb.2012.0021
- **Flye**: Kolmogorov et al. (2019). Nat Biotechnol. DOI: 10.1038/s41587-019-0072-8
- **Unicycler**: Wick et al. (2017). PLoS Comput Biol. DOI: 10.1371/journal.pcbi.1005595
- **QUAST**: Gurevich et al. (2013). Bioinformatics. DOI: 10.1093/bioinformatics/btt086
- **Prokka**: Seemann (2014). Bioinformatics. DOI: 10.1093/bioinformatics/btu153
- **AMRFinderPlus**: Feldgarden et al. (2021). Sci Rep. DOI: 10.1038/s41598-021-91456-0
- **Abricate**: Seemann T. [https://github.com/tseemann/abricate](https://github.com/tseemann/abricate)
- **RGI**: Alcock et al. (2020). Nucleic Acids Res. DOI: 10.1093/nar/gkz935

### Bases de Datos

- **NCBI Genome**: [https://www.ncbi.nlm.nih.gov/genome/](https://www.ncbi.nlm.nih.gov/genome/)
- **CARD**: [https://card.mcmaster.ca/](https://card.mcmaster.ca/)
- **ResFinder**: [https://cge.food.dtu.dk/services/ResFinder/](https://cge.food.dtu.dk/services/ResFinder/)
- **NCBI AMRFinderPlus Database**: [https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/)

### Tutoriales y Documentación

- **Bioconda**: [https://bioconda.github.io/](https://bioconda.github.io/)
- **Conda Documentation**: [https://docs.conda.io/](https://docs.conda.io/)
- **Mamba Documentation**: [https://mamba.readthedocs.io/](https://mamba.readthedocs.io/)

---

## 📄 Licencia

Este proyecto está bajo la licencia MIT. Ver archivo `LICENSE` para más detalles.

---

## 👥 Contribuciones

Las contribuciones son bienvenidas. Para contribuir:

1. Fork el repositorio
2. Crea una rama para tu feature (`git checkout -b feature/nueva-funcionalidad`)
3. Commit tus cambios (`git commit -am 'Añadir nueva funcionalidad'`)
4. Push a la rama (`git push origin feature/nueva-funcionalidad`)
5. Abre un Pull Request

---

## ✉️ Contacto y Soporte

**Autor**: [Tu Nombre]  
**Email**: [tu-email@institucion.edu]  
**Institución**: [Tu Institución]  
**GitHub**: [https://github.com/tu-usuario](https://github.com/tu-usuario)

Para reportar problemas o solicitar nuevas funcionalidades, abre un [Issue en GitHub](https://github.com/tu-usuario/Bacterial_Genomics_Project/issues).

---

## 🎓 Citación

Si utilizas este pipeline en tu investigación, por favor cita:

```
[Tu Nombre] (2025). Pipeline de Vigilancia Genómica y Análisis de Resistencia 
Antimicrobiana en Bacterias. GitHub repository: 
https://github.com/tu-usuario/Bacterial_Genomics_Project
```

---
