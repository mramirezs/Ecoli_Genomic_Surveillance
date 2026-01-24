# 📕 Pipeline Híbrido (Illumina + Nanopore)
### Análisis de Genomas Bacterianos con Ensamblaje Híbrido

---

## 📋 Tabla de Contenidos

1. [Introducción](#-introducción)
2. [Prerrequisitos](#-prerrequisitos)
3. [Visión General del Pipeline](#-visión-general-del-pipeline)
4. [Fase 1: Control de Calidad](#-fase-1-control-de-calidad)
5. [Fase 2: Ensamblaje Híbrido](#-fase-2-ensamblaje-híbrido)
6. [Fase 3: Evaluación Comparativa](#-fase-3-evaluación-comparativa)
7. [Fase 4: Mapeo y Validación](#-fase-4-mapeo-y-validación)
8. [Fase 5: Consenso de Alta Calidad](#-fase-5-consenso-de-alta-calidad)
9. [Fase 6: Análisis de Cobertura](#-fase-6-análisis-de-cobertura)
10. [Comparación de Estrategias](#-comparación-de-estrategias)
11. [Interpretación de Resultados](#-interpretación-de-resultados)
12. [Solución de Problemas](#-solución-de-problemas)

---

## 🎯 Introducción

### ¿Por Qué Usar Ensamblaje Híbrido?

El pipeline híbrido combina **lecturas cortas de Illumina** (alta precisión) con **lecturas largas de Nanopore** (alta continuidad) para producir ensamblajes de **calidad excepcional**.

### ⭐ Lo Mejor de Ambos Mundos

| Característica | Illumina Solo | Nanopore Solo | **Híbrido** |
|----------------|---------------|---------------|-------------|
| **Continuidad** | ⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
| **Precisión** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
| **Cromosoma cerrado** | ❌ | ✅ | ✅ |
| **Plásmidos cerrados** | ❌ | ✅ | ✅ |
| **SNP calling confiable** | ✅ | ⚠️ | ✅ |
| **Regiones repetitivas** | ❌ | ✅ | ✅ |
| **# Contigs** | 50-200 | 2-10 | 1-5 |
| **Tasa de errores** | <0.01% | ~2% | <0.01% |

### ✅ Cuándo Usar Este Pipeline

**Ideal para:**
- ✅ Genomas de referencia de alta calidad
- ✅ Publicaciones científicas
- ✅ Caracterización completa de estructura genómica
- ✅ Tipificación precisa de plásmidos
- ✅ Estudios de epidemiología molecular avanzada
- ✅ Cuando necesitas **lo mejor posible**

**Requisitos:**
- Datos Illumina paired-end (≥50x cobertura)
- Datos Nanopore (≥50x cobertura)
- Tiempo de cómputo mayor (5-8 horas)
- Mayor espacio en disco (~150-200 GB)

### 🎯 Resultados Esperados

Con un ensamblaje híbrido bien ejecutado obtendrás:
- **1-5 contigs** (cromosoma + plásmidos principales)
- **N50 >5 Mb** (cromosoma completo)
- **Precisión >99.99%** (corregido con Illumina)
- **Elementos circulares cerrados** (cromosoma + plásmidos)
- **SNPs confiables** (validados por ambas tecnologías)

---

## ✅ Prerrequisitos

### Antes de Empezar

- [ ] Instalación completa según [00_INSTALLATION.md](00_INSTALLATION.md)
- [ ] Ambiente `bact_main` activado
- [ ] **Datos Illumina** paired-end (≥50x cobertura)
- [ ] **Datos Nanopore** (≥50x cobertura)
- [ ] ~150-200 GB de espacio libre en disco
- [ ] Opcional pero recomendado: haber ejecutado pipelines individuales primero

### Verificar Instalación

```bash
# Activar ambiente
conda activate bact_main

# Verificar herramientas críticas para híbrido
unicycler --version
fastqc --version
NanoPlot --version
bwa
minimap2 --version
samtools --version

# Si todo está bien, continuar
```

### Estructura de Datos Esperada

```
00_raw_data/
├── illumina/
│   ├── SAMPLE_1.fastq.gz    # R1 (forward)
│   └── SAMPLE_2.fastq.gz    # R2 (reverse)
└── nanopore/
    └── SAMPLE_1.fastq.gz    # Long reads
```

**⚠️ IMPORTANTE**: Aunque los archivos pueden tener nombres similares, deben estar en directorios separados.

---

## 🔄 Visión General del Pipeline

```
┌──────────────────────────────────────────────────────────────┐
│                    PIPELINE HÍBRIDO                          │
└──────────────────────────────────────────────────────────────┘

1. DATOS CRUDOS
   ├─ Illumina R1/R2 (paired-end)
   └─ Nanopore (long reads)
   │
   ▼
2. CONTROL DE CALIDAD (PARALELO)
   ├─ Illumina: FastQC → fastp → FastQC
   └─ Nanopore: NanoPlot → Filtlong → NanoPlot
   │
   ▼
3. ENSAMBLAJE HÍBRIDO CON UNICYCLER
   ├─ Entrada: Illumina trimmed + Nanopore filtered
   ├─ Proceso: SPAdes + Miniasm + Racon + Pilon
   └─ Salida: Assembly híbrido de alta calidad
   │
   ▼
4. EVALUACIÓN COMPARATIVA (3-WAY)
   ├─ QUAST: Illumina vs Nanopore vs Híbrido
   └─ Identificar mejor ensamblaje
   │
   ▼
5. MAPEO CRUZADO
   ├─ Illumina → Ensamblaje híbrido
   ├─ Nanopore → Ensamblaje híbrido
   └─ Validación de ambas tecnologías
   │
   ▼
6. CONSENSO DE ALTA CALIDAD
   ├─ Variant calling con ambas tecnologías
   ├─ Validación cruzada de variantes
   └─ Secuencia consenso final
   │
   ▼
7. ANÁLISIS DE COBERTURA
   ├─ Cobertura Illumina por secuencia
   ├─ Cobertura Nanopore por secuencia
   └─ Validación de estructura
   │
   ▼
8. RESULTADOS FINALES
   ├─ Ensamblaje híbrido (1-5 contigs)
   ├─ Elementos circulares identificados
   ├─ Variantes validadas
   └─ Reportes integrados
```

**⏱️ Tiempo estimado total:** 5-8 horas  
**💾 Espacio requerido:** ~150-200 GB por muestra

---

## 🔬 Fase 1: Control de Calidad

### Objetivo

Realizar QC completo de ambas tecnologías antes del ensamblaje híbrido.

### Paso 1.1: QC de Datos Illumina

Si ya ejecutaste el pipeline Illumina ([01_ILLUMINA_PIPELINE.md](01_ILLUMINA_PIPELINE.md)), puedes reutilizar los datos limpios. Si no:

```bash
# Activar ambiente
conda activate bact_main

# Variables
SAMPLE="URO5550422"
R1="00_raw_data/illumina/${SAMPLE}_1.fastq.gz"
R2="00_raw_data/illumina/${SAMPLE}_2.fastq.gz"
THREADS=8

echo "========================================"
echo "QC Illumina - Pipeline Híbrido"
echo "Muestra: ${SAMPLE}"
echo "========================================"

# FastQC raw
mkdir -p 02_qc/01_illumina_raw
fastqc ${R1} ${R2} -o 02_qc/01_illumina_raw/ -t ${THREADS} -q

# Limpieza con fastp
mkdir -p 02_qc/02_illumina_trimmed
fastp \
  -i ${R1} -I ${R2} \
  -o 02_qc/02_illumina_trimmed/${SAMPLE}_R1_trimmed.fastq.gz \
  -O 02_qc/02_illumina_trimmed/${SAMPLE}_R2_trimmed.fastq.gz \
  --detect_adapter_for_pe --cut_front --cut_tail \
  --cut_window_size 4 --cut_mean_quality 20 --trim_poly_g \
  --qualified_quality_phred 20 --unqualified_percent_limit 40 \
  --n_base_limit 5 --length_required 50 --thread ${THREADS} \
  --html 02_qc/02_illumina_trimmed/${SAMPLE}_fastp_report.html \
  --json 02_qc/02_illumina_trimmed/${SAMPLE}_fastp_report.json

# FastQC trimmed
fastqc 02_qc/02_illumina_trimmed/${SAMPLE}_R*_trimmed.fastq.gz \
  -o 02_qc/02_illumina_trimmed/ -t ${THREADS} -q

echo "✓ QC Illumina completado"
```

### Paso 1.2: QC de Datos Nanopore

Si ya ejecutaste el pipeline Nanopore ([02_NANOPORE_PIPELINE.md](02_NANOPORE_PIPELINE.md)), puedes reutilizar los datos filtrados. Si no:

```bash
# Variables
NANOPORE="00_raw_data/nanopore/${SAMPLE}_1.fastq.gz"

echo "========================================"
echo "QC Nanopore - Pipeline Híbrido"
echo "========================================"

# NanoPlot raw
mkdir -p 02_qc/03_nanopore_raw
NanoPlot --fastq ${NANOPORE} \
  -o 02_qc/03_nanopore_raw/ -t ${THREADS} \
  --plots kde dot --N50 \
  --title "${SAMPLE} - Raw Nanopore"

# Filtlong
mkdir -p 02_qc/04_nanopore_filtered
filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 \
  ${NANOPORE} | \
  pigz -p ${THREADS} > 02_qc/04_nanopore_filtered/${SAMPLE}_ont_filtered.fastq.gz

# NanoPlot filtered
NanoPlot --fastq 02_qc/04_nanopore_filtered/${SAMPLE}_ont_filtered.fastq.gz \
  -o 02_qc/04_nanopore_filtered/ -t ${THREADS} \
  --plots kde dot --N50 \
  --title "${SAMPLE} - Filtered Nanopore"

echo "✓ QC Nanopore completado"
```

### Paso 1.3: Verificar Cobertura de Ambas Tecnologías

```bash
echo "========================================"
echo "Verificación de Cobertura"
echo "========================================"

GENOME_SIZE=5700000  # K. pneumoniae

# Cobertura Illumina
ILLUMINA_BASES=$(zcat 02_qc/02_illumina_trimmed/${SAMPLE}_R1_trimmed.fastq.gz | \
  paste - - - - | cut -f2 | tr -d '\n' | wc -c)
ILLUMINA_BASES=$((ILLUMINA_BASES * 2))  # x2 porque es paired-end
ILLUMINA_COV=$(echo "$ILLUMINA_BASES / $GENOME_SIZE" | bc)

echo "Illumina:"
echo "  Total bases: $(printf "%'d" $ILLUMINA_BASES)"
echo "  Cobertura estimada: ${ILLUMINA_COV}x"

# Cobertura Nanopore
NANOPORE_BASES=$(grep "Total bases:" 02_qc/04_nanopore_filtered/NanoStats.txt | \
  awk '{print $NF}' | tr -d ',')
NANOPORE_COV=$(echo "$NANOPORE_BASES / $GENOME_SIZE" | bc)

echo ""
echo "Nanopore:"
echo "  Total bases: $(printf "%'d" $NANOPORE_BASES)"
echo "  Cobertura estimada: ${NANOPORE_COV}x"

echo ""
echo "=== EVALUACIÓN ==="
if [ $ILLUMINA_COV -ge 50 ] && [ $NANOPORE_COV -ge 50 ]; then
    echo "✓ Coberturas adecuadas para ensamblaje híbrido"
else
    echo "⚠️  Advertencia: cobertura baja detectada"
    [ $ILLUMINA_COV -lt 50 ] && echo "  - Illumina: ${ILLUMINA_COV}x (recomendado ≥50x)"
    [ $NANOPORE_COV -lt 50 ] && echo "  - Nanopore: ${NANOPORE_COV}x (recomendado ≥50x)"
fi
```

### Paso 1.4: Reporte MultiQC Integrado

```bash
echo "========================================"
echo "Reporte MultiQC Integrado"
echo "========================================"

mkdir -p 02_qc/05_multiqc

multiqc 02_qc/ \
  -o 02_qc/05_multiqc/ \
  --filename ${SAMPLE}_hybrid_multiqc_report \
  --title "Hybrid QC Report - ${SAMPLE}" \
  --comment "Illumina + Nanopore for hybrid assembly" \
  --force

echo "✓ Reporte MultiQC generado"
firefox 02_qc/05_multiqc/${SAMPLE}_hybrid_multiqc_report.html &
```

---

## 🧬 Fase 2: Ensamblaje Híbrido

### Objetivo

Ensamblar el genoma usando Unicycler, que integra inteligentemente lecturas cortas y largas.

### ¿Cómo Funciona Unicycler?

Unicycler ejecuta varios pasos automáticamente:

1. **Ensamblaje inicial con SPAdes** (usando Illumina)
2. **Bridging con lecturas largas** (Nanopore cierra gaps)
3. **Polishing con Racon** (corrige errores Nanopore)
4. **Polishing final con Pilon** (usa Illumina para máxima precisión)

### Paso 2.1: Ensamblaje Híbrido con Unicycler

```bash
echo "========================================"
echo "Ensamblaje Híbrido con Unicycler"
echo "Muestra: ${SAMPLE}"
echo "Inicio: $(date)"
echo "========================================"

# Variables
R1_TRIM="02_qc/02_illumina_trimmed/${SAMPLE}_R1_trimmed.fastq.gz"
R2_TRIM="02_qc/02_illumina_trimmed/${SAMPLE}_R2_trimmed.fastq.gz"
NANOPORE_FILT="02_qc/04_nanopore_filtered/${SAMPLE}_ont_filtered.fastq.gz"
THREADS=8

# Crear directorio
mkdir -p 03_assembly/03_hybrid

# Ejecutar Unicycler
unicycler \
  -1 ${R1_TRIM} \
  -2 ${R2_TRIM} \
  -l ${NANOPORE_FILT} \
  -o 03_assembly/03_hybrid/ \
  --threads ${THREADS} \
  --mode normal \
  --min_fasta_length 200 \
  --keep 2

echo "✓ Ensamblaje híbrido completado"
echo "  Fin: $(date)"
```

**⚙️ Parámetros de Unicycler:**

| Parámetro | Función |
|-----------|---------|
| `-1 / -2` | Lecturas Illumina paired-end |
| `-l` | Lecturas largas (Nanopore) |
| `--mode normal` | Balance entre velocidad y calidad |
| `--min_fasta_length 200` | Descartar contigs <200 bp |
| `--keep 2` | Guardar archivos intermedios (nivel medio) |

**🎯 Modos de Unicycler:**

- `--mode conservative`: Más lento, máxima calidad (usar para publicaciones)
- `--mode normal`: Balance (recomendado para mayoría)
- `--mode bold`: Más rápido, puede ser menos preciso

### Paso 2.2: Archivos Generados por Unicycler

```bash
echo "========================================"
echo "Archivos Generados"
echo "========================================"

ls -lh 03_assembly/03_hybrid/

# Archivos principales:
# assembly.fasta - Ensamblaje final (USAR ESTE)
# assembly.gfa - Grafo de ensamblaje
# unicycler.log - Log detallado
```

**📁 Archivos importantes:**

```
03_assembly/03_hybrid/
├── assembly.fasta           # ⭐ ENSAMBLAJE FINAL
├── assembly.gfa             # Grafo (visualizar con Bandage)
├── unicycler.log            # Log del proceso
└── [varios archivos SAM/BAM de polishing]
```

### Paso 2.3: Renombrar y Analizar Ensamblaje

```bash
# Copiar con nombre estándar
cp 03_assembly/03_hybrid/assembly.fasta \
   03_assembly/03_hybrid/${SAMPLE}_hybrid_assembly.fasta

echo "========================================"
echo "Estadísticas del Ensamblaje Híbrido"
echo "========================================"

ASSEMBLY="03_assembly/03_hybrid/${SAMPLE}_hybrid_assembly.fasta"

# Número de contigs
echo -n "Número de contigs: "
grep -c ">" ${ASSEMBLY}

# Tamaño total
echo -n "Tamaño total: "
grep -v ">" ${ASSEMBLY} | tr -d '\n' | wc -c | awk '{printf "%'"'"'d bp\n", $1}'

# Contigs y longitudes
echo ""
echo "=== CONTIGS ENSAMBLADOS ==="
grep ">" ${ASSEMBLY} | while read header; do
    contig_name=$(echo $header | sed 's/>//' | awk '{print $1}')
    length=$(echo $header | grep -oP 'length=\K[0-9]+')
    depth=$(echo $header | grep -oP 'depth=\K[0-9.]+')
    circular=$(echo $header | grep -o 'circular=true' || echo "linear")
    
    printf "%-15s %12s bp  Depth: %6sx  %s\n" \
      "$contig_name" "$(printf "%'d" $length)" "$depth" "$circular"
done
```

**📊 Salida Esperada:**

```
Número de contigs: 4

=== CONTIGS ENSAMBLADOS ===
1              5,334,567 bp  Depth:   65.2x  circular=true
2                122,799 bp  Depth:   54.1x  circular=true
3                111,195 bp  Depth:   48.7x  circular=true
4                105,974 bp  Depth:   51.3x  circular=true
```

---

## 📊 Fase 3: Evaluación Comparativa

### Objetivo

Comparar el ensamblaje híbrido contra los ensamblajes individuales (Illumina y Nanopore) para validar la mejora.

### Paso 3.1: Preparar Ensamblajes para Comparación

```bash
echo "========================================"
echo "Preparando Comparación 3-Way"
echo "========================================"

# Si NO tienes los ensamblajes individuales, crearlos
# (omitir si ya los ejecutaste)

# Ensamblaje Illumina (si no existe)
if [ ! -f "03_assembly/01_illumina_only/${SAMPLE}_illumina_assembly.fasta" ]; then
    echo "Ejecutando ensamblaje Illumina..."
    bash scripts/run_illumina_assembly_only.sh ${SAMPLE}
fi

# Ensamblaje Nanopore (si no existe)
if [ ! -f "03_assembly/02_nanopore_only/${SAMPLE}_nanopore_polished.fasta" ]; then
    echo "Ejecutando ensamblaje Nanopore..."
    bash scripts/run_nanopore_assembly_only.sh ${SAMPLE}
fi
```

### Paso 3.2: Evaluación con QUAST (3-Way)

```bash
echo "========================================"
echo "Evaluación QUAST - Comparación 3-Way"
echo "========================================"

ILLUMINA_ASM="03_assembly/01_illumina_only/${SAMPLE}_illumina_assembly.fasta"
NANOPORE_ASM="03_assembly/02_nanopore_only/${SAMPLE}_nanopore_polished.fasta"
HYBRID_ASM="03_assembly/03_hybrid/${SAMPLE}_hybrid_assembly.fasta"
REFERENCE="01_reference/reference.fasta"

mkdir -p 03_assembly/04_quast_evaluation/hybrid_comparison

quast.py \
  ${ILLUMINA_ASM} \
  ${NANOPORE_ASM} \
  ${HYBRID_ASM} \
  -r ${REFERENCE} \
  -o 03_assembly/04_quast_evaluation/hybrid_comparison/ \
  --threads ${THREADS} \
  --labels "Illumina,Nanopore,Hybrid" \
  --glimmer \
  --min-contig 200 \
  --plots-format png \
  --circos

echo "✓ Evaluación QUAST completada"
firefox 03_assembly/04_quast_evaluation/hybrid_comparison/report.html &
```

### Paso 3.3: Tabla Comparativa

```bash
echo "========================================"
echo "Tabla Comparativa - 3 Estrategias"
echo "========================================"

# Mostrar reporte en terminal
cat 03_assembly/04_quast_evaluation/hybrid_comparison/report.txt

# Generar tabla resumida
cat > 03_assembly/04_quast_evaluation/hybrid_comparison/summary_table.txt << EOF
# Comparación de Estrategias de Ensamblaje
# Muestra: ${SAMPLE}
# Fecha: $(date)

EOF

echo "Métrica                         | Illumina  | Nanopore | Híbrido  | Mejor" >> \
  03_assembly/04_quast_evaluation/hybrid_comparison/summary_table.txt
echo "--------------------------------|-----------|----------|----------|-------" >> \
  03_assembly/04_quast_evaluation/hybrid_comparison/summary_table.txt

# Extraer métricas clave
grep "# contigs (>= 0 bp)" 03_assembly/04_quast_evaluation/hybrid_comparison/report.txt | \
  awk '{printf "%-31s | %9s | %8s | %8s |\n", "# contigs", $4, $5, $6}' >> \
  03_assembly/04_quast_evaluation/hybrid_comparison/summary_table.txt

grep "Largest contig" 03_assembly/04_quast_evaluation/hybrid_comparison/report.txt | \
  awk '{printf "%-31s | %9s | %8s | %8s |\n", "Largest contig", $3, $4, $5}' >> \
  03_assembly/04_quast_evaluation/hybrid_comparison/summary_table.txt

grep "Total length" 03_assembly/04_quast_evaluation/hybrid_comparison/report.txt | head -1 | \
  awk '{printf "%-31s | %9s | %8s | %8s |\n", "Total length", $3, $4, $5}' >> \
  03_assembly/04_quast_evaluation/hybrid_comparison/summary_table.txt

grep "N50" 03_assembly/04_quast_evaluation/hybrid_comparison/report.txt | head -1 | \
  awk '{printf "%-31s | %9s | %8s | %8s |\n", "N50", $2, $3, $4}' >> \
  03_assembly/04_quast_evaluation/hybrid_comparison/summary_table.txt

grep "L50" 03_assembly/04_quast_evaluation/hybrid_comparison/report.txt | head -1 | \
  awk '{printf "%-31s | %9s | %8s | %8s |\n", "L50", $2, $3, $4}' >> \
  03_assembly/04_quast_evaluation/hybrid_comparison/summary_table.txt

# Mostrar tabla
cat 03_assembly/04_quast_evaluation/hybrid_comparison/summary_table.txt
```

**📊 Ejemplo de Tabla Comparativa:**

```
Métrica                         | Illumina  | Nanopore | Híbrido  | Mejor
--------------------------------|-----------|----------|----------|-------
# contigs                       |        98 |        7 |        4 | Híbrido
Largest contig                  |   387,234 | 5,334,567| 5,334,567| Híbrido
Total length                    | 5,612,345 | 5,723,892| 5,689,234| Híbrido
N50                             |   145,678 | 5,334,567| 5,334,567| Híbrido
L50                             |        12 |        1 |        1 | Híbrido
```

---

## 🗺️ Fase 4: Mapeo y Validación

### Objetivo

Mapear AMBAS tecnologías contra el ensamblaje híbrido para validación cruzada.

### Paso 4.1: Preparar Ensamblaje Híbrido como Referencia

```bash
echo "========================================"
echo "Preparando Ensamblaje Híbrido"
echo "========================================"

HYBRID_ASM="03_assembly/03_hybrid/${SAMPLE}_hybrid_assembly.fasta"

# Índices para mapeo
bwa index ${HYBRID_ASM}
samtools faidx ${HYBRID_ASM}

echo "✓ Índices creados"
```

### Paso 4.2: Mapeo de Lecturas Illumina

```bash
echo "========================================"
echo "Mapeo Illumina → Ensamblaje Híbrido"
echo "========================================"

R1_TRIM="02_qc/02_illumina_trimmed/${SAMPLE}_R1_trimmed.fastq.gz"
R2_TRIM="02_qc/02_illumina_trimmed/${SAMPLE}_R2_trimmed.fastq.gz"

mkdir -p 04_mapping/01_illumina

# Mapeo
bwa mem -t ${THREADS} \
  -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
  ${HYBRID_ASM} ${R1_TRIM} ${R2_TRIM} | \
  samtools view -Sb - | \
  samtools sort -@ ${THREADS} -o 04_mapping/01_illumina/${SAMPLE}_hybrid_sorted.bam

# Indexar
samtools index 04_mapping/01_illumina/${SAMPLE}_hybrid_sorted.bam

# Estadísticas
samtools flagstat 04_mapping/01_illumina/${SAMPLE}_hybrid_sorted.bam > \
  04_mapping/01_illumina/${SAMPLE}_hybrid_flagstat.txt

samtools coverage 04_mapping/01_illumina/${SAMPLE}_hybrid_sorted.bam > \
  04_mapping/01_illumina/${SAMPLE}_hybrid_coverage.txt

echo "✓ Mapeo Illumina completado"
cat 04_mapping/01_illumina/${SAMPLE}_hybrid_flagstat.txt
```

### Paso 4.3: Mapeo de Lecturas Nanopore

```bash
echo "========================================"
echo "Mapeo Nanopore → Ensamblaje Híbrido"
echo "========================================"

NANOPORE_FILT="02_qc/04_nanopore_filtered/${SAMPLE}_ont_filtered.fastq.gz"

mkdir -p 04_mapping/02_nanopore

# Mapeo
minimap2 -ax map-ont -t ${THREADS} \
  ${HYBRID_ASM} ${NANOPORE_FILT} | \
  samtools view -Sb - | \
  samtools sort -@ ${THREADS} -o 04_mapping/02_nanopore/${SAMPLE}_hybrid_sorted.bam

# Indexar
samtools index 04_mapping/02_nanopore/${SAMPLE}_hybrid_sorted.bam

# Estadísticas
samtools flagstat 04_mapping/02_nanopore/${SAMPLE}_hybrid_sorted.bam > \
  04_mapping/02_nanopore/${SAMPLE}_hybrid_flagstat.txt

samtools coverage 04_mapping/02_nanopore/${SAMPLE}_hybrid_sorted.bam > \
  04_mapping/02_nanopore/${SAMPLE}_hybrid_coverage.txt

echo "✓ Mapeo Nanopore completado"
cat 04_mapping/02_nanopore/${SAMPLE}_hybrid_flagstat.txt
```

---

## 🎯 Fase 5: Consenso de Alta Calidad

### Objetivo

Generar secuencia consenso validada por ambas tecnologías.

### Paso 5.1: Variant Calling con Illumina

```bash
echo "========================================"
echo "Variant Calling - Illumina"
echo "========================================"

BAM_ILLUMINA="04_mapping/01_illumina/${SAMPLE}_hybrid_sorted.bam"

mkdir -p 04_mapping/03_variants

# Call variants
bcftools
