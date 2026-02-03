# 🧬 Bacterial Genomics Pipeline
### Análisis Completo de Genomas Bacterianos con NGS

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Bioinformatics](https://img.shields.io/badge/Bioinformatics-Pipeline-blue.svg)]()
[![Status](https://img.shields.io/badge/Status-Production-green.svg)]()

---

## 📋 Descripción

Pipeline modular para análisis de genomas bacterianos utilizando datos de secuenciación de nueva generación (NGS). Soporta **tres estrategias independientes** de ensamblaje según los datos disponibles:

- 🔵 **Solo Illumina** - Lecturas cortas de alta precisión
- 🟢 **Solo Nanopore** - Lecturas largas para mayor continuidad  
- 🟣 **Híbrido** - Combina ambas tecnologías (recomendado)

Además incluye análisis exhaustivo de **resistencia antimicrobiana (AMR)**, anotación funcional y tipificación molecular.

---

## 🚀 Inicio Rápido

### ¿Qué tipo de datos tienes?

| Tus Datos | Pipeline Recomendado | Tiempo Estimado | Ir a Documentación |
|-----------|---------------------|-----------------|-------------------|
| 📘 Solo Illumina | Pipeline Illumina | 3-5 horas | [Ver guía →](docs/01_ILLUMINA_PIPELINE.md) |
| 📗 Solo Nanopore | Pipeline Nanopore | 2-4 horas | [Ver guía →](docs/02_NANOPORE_PIPELINE.md) |
| 📕 Illumina + Nanopore | **Pipeline Híbrido** ⭐ | 5-8 horas | [Ver guía →](docs/03_HYBRID_PIPELINE.md) |

> **💡 Recomendación**: Si tienes ambos tipos de datos, usa el pipeline híbrido para obtener la **mejor calidad** (continuidad de Nanopore + precisión de Illumina).

---

## ⚡ Instalación Rápida (3 Pasos)

### Paso 1: Configurar Estructura del Proyecto

```bash
# Clonar repositorio
git clone https://github.com/tu-usuario/Bacterial_Genomics_Pipeline.git
cd Bacterial_Genomics_Pipeline

# Ejecutar configuración automática
bash setup_project_structure.sh

# O personalizar nombre del proyecto
bash setup_project_structure.sh mi_proyecto URO5550422
```

**Esto crea:**
- ✅ Estructura completa de 14 directorios principales
- ✅ 40+ subdirectorios organizados
- ✅ Descarga genoma de referencia K. pneumoniae
- ✅ Archivos de configuración y metadata
- ✅ Scripts auxiliares

📚 **Guía completa:** [SETUP_PROJECT_GUIDE.md](docs/SETUP_PROJECT_GUIDE.md)

### Paso 2: Instalar Ambientes Conda

```bash
# Configurar ambientes especializados (~45 min)
bash scripts/setup_environments.sh

# Verificar instalación
bash scripts/verify_installation.sh
```

📚 **Guía completa:** [00_INSTALLATION.md](docs/00_INSTALLATION.md)

### Paso 3: Agregar tus Datos y Ejecutar

```bash
# Enlazar datos de secuenciación
bash scripts/link_raw_data.sh /ruta/illumina /ruta/nanopore

# Ejecutar pipeline según tus datos
bash scripts/run_hybrid_pipeline.sh URO5550422

# Ver resultados
firefox 08_results/FINAL_REPORT.html
```

---

## 📚 Documentación Completa

### 🛠️ Configuración e Instalación

| Documento | Descripción | Tiempo |
|-----------|-------------|--------|
| **[SETUP_PROJECT_GUIDE.md](docs/SETUP_PROJECT_GUIDE.md)** | Configuración automática de estructura | ~5 min |
| **[00_INSTALLATION.md](docs/00_INSTALLATION.md)** | Instalación completa de ambientes y bases de datos | ~45 min |

### 🔬 Pipelines de Análisis

| Pipeline | Descripción | Cuándo Usar | Documentación |
|----------|-------------|-------------|---------------|
| **📘 Illumina** | Ensamblaje con lecturas cortas | Solo tienes datos Illumina | [01_ILLUMINA_PIPELINE.md](docs/01_ILLUMINA_PIPELINE.md) |
| **📗 Nanopore** | Ensamblaje con lecturas largas | Solo tienes datos Nanopore | [02_NANOPORE_PIPELINE.md](docs/02_NANOPORE_PIPELINE.md) |
| **📕 Híbrido** ⭐ | Combina Illumina + Nanopore | Tienes ambos tipos (mejor calidad) | [03_HYBRID_PIPELINE.md](docs/03_HYBRID_PIPELINE.md) |

### 🛡️ Análisis Downstream

| Documento | Descripción |
|-----------|-------------|
| **[04_AMR_TYPING.md](docs/04_AMR_TYPING.md)** | Detección AMR, anotación, MLST, plásmidos |
| **[05_TROUBLESHOOTING.md](docs/05_TROUBLESHOOTING.md)** | Solución de problemas comunes |

---

## 📊 ¿Qué Puedo Hacer con Este Pipeline?

✅ **Ensamblar genomas bacterianos** de alta calidad  
✅ **Identificar genes de resistencia** a antibióticos (AMR)  
✅ **Detectar variantes genómicas** (SNPs, INDELs)  
✅ **Anotar genes y funciones** biológicas  
✅ **Comparar diferentes estrategias** de ensamblaje  
✅ **Analizar cromosomas y plásmidos** por separado  
✅ **Tipificar cepas** (MLST, detección de plásmidos)  
✅ **Generar reportes automatizados** para vigilancia epidemiológica  

---

## 🎯 Caso de Estudio: *Klebsiella pneumoniae* URO5550422

Todos los pipelines están documentados usando un caso real:

- **Organismo:** *Klebsiella pneumoniae*
- **Muestra:** URO5550422 (aislado clínico urinario)
- **Referencia:** K. pneumoniae HS11286 (GCF_000240185.1)
- **Genoma:** 5.7 Mb (1 cromosoma + 6 plásmidos)
- **Datos disponibles:** Illumina paired-end + Nanopore long-reads

---

## 📂 Estructura del Proyecto

Después de ejecutar `setup_project_structure.sh`:

```
bacterial_genomics/
├── 00_raw_data/          # Datos de secuenciación (FASTQ)
│   ├── illumina/         # Lecturas paired-end
│   └── nanopore/         # Lecturas largas
├── 01_reference/         # Genoma de referencia
├── 02_qc/                # Control de calidad
├── 03_assembly/          # Ensamblajes (Illumina/Nanopore/Híbrido)
├── 04_mapping/           # Mapeos y variantes
├── 05_annotation/        # Anotación funcional
├── 06_amr_screening/     # Genes de resistencia AMR
├── 07_typing/            # Tipificación molecular (MLST, plásmidos)
├── 08_results/           # Resultados finales y reportes
├── databases/            # Bases de datos locales
├── envs/                 # Ambientes conda (YAML)
├── scripts/              # Scripts de automatización
└── logs/                 # Logs de ejecución
```

---

## 💻 Requisitos del Sistema

### Hardware Recomendado

| Componente | Mínimo | Recomendado | Óptimo |
|------------|--------|-------------|--------|
| **CPU** | 4 cores | 8 cores | 16+ cores |
| **RAM** | 16 GB | 32 GB | 64+ GB |
| **Almacenamiento** | 100 GB/muestra | 200 GB/muestra | SSD 500 GB |
| **Red** | 10 Mbps | 100 Mbps | 1 Gbps |

### Software

- **Sistema Operativo**: Linux/Unix (Ubuntu 20.04+, CentOS 7+, Debian 10+)
- **Conda/Mamba**: Para gestión de ambientes
- **Git**: Para clonar repositorio
- **Conexión a internet**: Para instalación inicial y descarga de bases de datos

---

## 📊 Comparación de Estrategias

| Característica | Illumina | Nanopore | Híbrido |
|---------------|----------|----------|---------|
| **Número de contigs** | 50-150 | 2-10 | 1-10 |
| **N50** | 100-300 kb | 5+ Mb | 5+ Mb |
| **Precisión** | >99.9% | ~95-98% | >99.99% |
| **Continuidad** | Baja | Alta | Alta |
| **Plásmidos cerrados** | No | Sí | Sí |
| **Costo computacional** | Bajo | Medio | Alto |
| **Tiempo ejecución** | 3-5h | 2-4h | 5-8h |
| **SNP calling** | Excelente | Regular | Excelente |
| **Mejor para** | Variantes | Estructura | Todo |

---

## 🔄 Flujo de Trabajo General

```
┌─────────────────────────────────────────────┐
│  1. Configurar Estructura del Proyecto     │
│     bash setup_project_structure.sh        │
└─────────────────┬───────────────────────────┘
                  │
┌─────────────────▼───────────────────────────┐
│  2. Instalar Ambientes Conda               │
│     bash scripts/setup_environments.sh     │
└─────────────────┬───────────────────────────┘
                  │
┌─────────────────▼───────────────────────────┐
│  3. Agregar Datos de Secuenciación         │
│     bash scripts/link_raw_data.sh          │
└─────────────────┬───────────────────────────┘
                  │
┌─────────────────▼───────────────────────────┐
│  4. Elegir Pipeline                        │
│     • Illumina                             │
│     • Nanopore                             │
│     • Híbrido ⭐                           │
└─────────────────┬───────────────────────────┘
                  │
┌─────────────────▼───────────────────────────┐
│  5. Control de Calidad                     │
│     FastQC, fastp, NanoPlot                │
└─────────────────┬───────────────────────────┘
                  │
┌─────────────────▼───────────────────────────┐
│  6. Ensamblaje                             │
│     SPAdes / Flye / Unicycler              │
└─────────────────┬───────────────────────────┘
                  │
┌─────────────────▼───────────────────────────┐
│  7. Análisis Downstream                    │
│     • Anotación (Prokka)                   │
│     • Detección AMR                        │
│     • MLST typing                          │
│     • Plásmidos                            │
└─────────────────┬───────────────────────────┘
                  │
┌─────────────────▼───────────────────────────┐
│  8. Resultados y Reportes                  │
│     Visualización, tablas, gráficos        │
└─────────────────────────────────────────────┘
```

---

## ✅ Checklist de Decisión

### ¿Qué pipeline debo usar?

- [ ] **¿Tengo datos Illumina paired-end?**
  - Sí → Puedes usar pipeline Illumina o Híbrido
  - No → Usa pipeline Nanopore

- [ ] **¿Tengo datos Nanopore long-reads?**
  - Sí → Puedes usar pipeline Nanopore o Híbrido
  - No → Usa pipeline Illumina

- [ ] **¿Tengo AMBOS tipos de datos?**
  - ✅ Sí → **USA PIPELINE HÍBRIDO** (mejor opción)

- [ ] **¿Necesito plásmidos cerrados?**
  - Sí → Requiere Nanopore o Híbrido
  - No → Illumina es suficiente

- [ ] **¿Priorizo precisión en SNPs?**
  - Sí → Illumina o Híbrido
  - No → Nanopore puede ser suficiente

---

## 🎓 Para Empezar

### Usuarios Nuevos

1. **Configurar estructura:** `bash setup_project_structure.sh`
2. **Instalar ambientes:** Ver [00_INSTALLATION.md](docs/00_INSTALLATION.md)
3. **Agregar datos:** `bash scripts/link_raw_data.sh`
4. **Elegir pipeline:** Según datos disponibles
5. **Ejecutar análisis:** Seguir guía del pipeline elegido
6. **Analizar resultados:** Detección AMR y tipificación

### Usuarios Avanzados

- Revisar documentación específica de tu pipeline
- Modificar scripts según necesidades
- Integrar con tus propios workflows
- Contribuir con mejoras (pull requests bienvenidos)

---

## 📖 Referencias y Recursos

### Herramientas Principales

- **FastQC/fastp:** Control de calidad
- **SPAdes:** Ensamblaje Illumina
- **Flye:** Ensamblaje Nanopore
- **Unicycler:** Ensamblaje híbrido
- **BWA/Minimap2:** Mapeo de lecturas
- **Prokka:** Anotación funcional
- **AMRFinderPlus/CARD:** Detección AMR

### Bases de Datos

- NCBI RefSeq
- CARD (Comprehensive Antibiotic Resistance Database)
- ResFinder
- VFDB (Virulence Factors)
- PubMLST

### Publicaciones Clave

- Wick et al. (2017) - Unicycler: https://doi.org/10.1371/journal.pcbi.1005595
- Kolmogorov et al. (2019) - Flye: https://doi.org/10.1038/s41587-019-0072-8
- Bankevich et al. (2012) - SPAdes: https://doi.org/10.1089/cmb.2012.0021

---

## 🤝 Contribuir

¿Encontraste un bug? ¿Tienes una sugerencia?

1. Abre un **Issue** describiendo el problema
2. Envía un **Pull Request** con mejoras
3. Comparte tus casos de uso
4. Ayuda a mejorar la documentación

---

## 📄 Licencia

Este proyecto está licenciado bajo MIT License - ver archivo [LICENSE](LICENSE)

---

## 📧 Contacto y Soporte

- **Issues:** [GitHub Issues](https://github.com/tu-usuario/Bacterial_Genomics_Pipeline/issues)
- **Discusiones:** [GitHub Discussions](https://github.com/tu-usuario/Bacterial_Genomics_Pipeline/discussions)
- **Email:** tu-email@ejemplo.com

---

## 🌟 Agradecimientos

Este pipeline integra herramientas desarrolladas por la comunidad científica y bioinformática. Agradecemos a todos los desarrolladores de:

- Bioconda project
- Galaxy project  
- NCBI
- CARD
- PubMLST
- Y todos los creadores de herramientas open-source

---

<div align="center">

**¿Listo para empezar?**

[🛠️ Configurar Proyecto](docs/SETUP_PROJECT_GUIDE.md) | [📚 Instalación](docs/00_INSTALLATION.md) | [📘 Pipeline Illumina](docs/01_ILLUMINA_PIPELINE.md) | [📗 Pipeline Nanopore](docs/02_NANOPORE_PIPELINE.md) | [📕 Pipeline Híbrido](docs/03_HYBRID_PIPELINE.md)

---

⭐ **Si este proyecto te fue útil, considera darle una estrella en GitHub** ⭐

</div>
