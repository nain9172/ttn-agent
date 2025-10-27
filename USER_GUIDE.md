# TTN Variant AI Agent - User Guide

## Table of Contents
1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Output Interpretation](#output-interpretation)
5. [Advanced Features](#advanced-features)
6. [Troubleshooting](#troubleshooting)
7. [API Reference](#api-reference)

## Introduction

The TTN Variant AI Agent is a comprehensive tool for analyzing variants in the TTN (Titin) gene. It combines:
- **AI-powered pathogenicity prediction** using Evo2
- **Literature mining** from PubMed
- **Protein domain visualization**
- **Automated HTML report generation**

### What Makes This Tool Unique?

1. **No Manual Searching**: Automatically finds relevant literature
2. **Visual Representation**: See exactly where variants are located
3. **AI Prediction**: State-of-the-art pathogenicity scoring
4. **Beautiful Reports**: Professional HTML reports ready to share

## Installation

### Prerequisites
- Python 3.8 or higher
- 16GB RAM (recommended for Evo2)
- Internet connection
- ~10GB disk space (for reference genome)

### Quick Install

```bash
# Clone repository
git clone <repository>
cd TTN_Variant_AI

# Run setup script
chmod +x setup.sh
./setup.sh

# Activate environment
source venv/bin/activate

# Test installation
python test_agent.py
```

### Manual Install

```bash
# Create virtual environment
python3 -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
pip install git+https://github.com/ArcInstitute/evo2.git

# Download reference genome
cd data
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz
cd ..
```

## Usage

### Basic Usage

Analyze a single variant:
```bash
python main.py 2-178612477-T-A
```

### Input Format

Variants can be specified in multiple formats:
- `2-178612477-T-A` (hyphen separator)
- `2_178612477_T_A` (underscore separator)
- `2:178612477:T:A` (colon separator)

Format: `chromosome-position-reference-alternate`

### Command Line Options

```bash
python main.py <variant> [OPTIONS]

Options:
  --output PATH       Output HTML file path
  --skip-evo2        Skip Evo2 prediction (faster)
  --help             Show help message
```

### Examples

**Standard Analysis**:
```bash
python main.py 2-178527121-A-G
```

**Fast Analysis (skip Evo2)**:
```bash
python main.py 2-178527121-A-G --skip-evo2
```

**Custom Output**:
```bash
python main.py 2-178527121-A-G --output my_report.html
```

### Batch Processing

Process multiple variants from a file:

```bash
python batch_analyze.py example_variants.csv
```

Input file format (CSV):
```csv
variant_id
2-178527121-A-G
2-178527148-A-T
2-178527628-C-T
```

Or simple text file (one per line):
```
2-178527121-A-G
2-178527148-A-T
2-178527628-C-T
```

## Output Interpretation

### HTML Report Sections

#### 1. Variant Summary
Basic information about the variant:
- Genomic coordinates (GRCh38)
- Reference and alternate alleles
- Gene information

#### 2. Evo2 Pathogenicity Prediction

**Delta Score Interpretation**:
- **< 0**: Variant likely disrupts protein function (PATHOGENIC)
- **≥ 0**: Variant likely maintains protein function (BENIGN)

**Confidence Levels**:
| Absolute Delta Score | Confidence |
|---------------------|------------|
| > 0.005 | High |
| 0.001 - 0.005 | Medium |
| < 0.001 | Low |

**Example Interpretations**:
```
Delta Score: -0.002667
Prediction: PATHOGENIC
Confidence: Medium (0.002667)
→ Variant likely deleterious

Delta Score: 0.000735
Prediction: BENIGN
Confidence: Low (0.000735)
→ Variant likely tolerated
```

#### 3. Protein Domain Localization

Visual representation showing:
- **Z-disk** (Red): N-terminal region, sarcomere anchoring
- **I-band** (Cyan): Elastic region
- **A-band** (Blue): Thick filament binding
- **M-band** (Green): C-terminal region

The figure shows:
1. **Domain Overview**: Variant location across all domains
2. **Transcript Views**: Position in different isoforms (N2BA, N2B, Novex-1/2/3)

#### 4. Literature Review

PubMed search results with extracted information:
- **Phenotype**: Cardiac, skeletal muscle, or both
- **Inheritance**: AD, AR, X-linked, or de novo
- **Age of Onset**: Congenital, infantile, childhood, adult, or late-onset
- **Links**: Direct links to PubMed articles

### Output Files

After running, you'll find:
```
outputs/
├── variant_report_TIMESTAMP.html    # Main HTML report
├── images/
│   └── titin_schematic_*.png       # Protein schematic
└── ttn_agent.log                    # Detailed logs
```

## Advanced Features

### Custom Configuration

Edit `config.py` to customize:

```python
# Model selection
EVO2_MODEL = "evo2_1b_base"  # Fast, less accurate
EVO2_MODEL = "evo2_7b_base"  # Slow, more accurate

# PubMed settings
PUBMED_MAX_RESULTS = 20      # Number of articles to retrieve
PUBMED_EMAIL = "your@email.com"

# Report styling
REPORT_HEADER_COLOR = "#2C3E50"
REPORT_ACCENT_COLOR = "#3498DB"
```

### Using as a Python Module

```python
from utils import (
    parse_variant,
    Evo2Predictor,
    PubMedSearcher,
    ImageGenerator,
    HTMLReportGenerator
)

# Parse variant
variant_info = parse_variant("2-178612477-T-A")

# Predict pathogenicity
predictor = Evo2Predictor()
result = predictor.predict(variant_info)

print(f"Prediction: {result['prediction']}")
print(f"Delta Score: {result['delta_score']}")
```

### Programmatic Batch Processing

```python
import pandas as pd
from utils import parse_variant, Evo2Predictor

# Load variants
variants_df = pd.read_csv("variants.csv")

# Initialize predictor
predictor = Evo2Predictor()

# Process variants
results = []
for variant_id in variants_df['variant_id']:
    variant_info = parse_variant(variant_id)
    result = predictor.predict(variant_info)
    results.append({
        'variant': variant_id,
        'prediction': result['prediction'],
        'delta_score': result['delta_score']
    })

# Save results
results_df = pd.DataFrame(results)
results_df.to_csv("predictions.csv", index=False)
```

## Troubleshooting

### Common Issues

#### "Reference genome not found"

**Solution**: Download the reference genome:
```bash
cd data
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip *.gz
```

#### "Out of memory" during Evo2 prediction

**Solution 1**: Use smaller model
```python
# In config.py
EVO2_MODEL = "evo2_1b_base"
```

**Solution 2**: Skip Evo2
```bash
python main.py <variant> --skip-evo2
```

#### "PubMed search failed"

**Solution**: Check internet and set email
```bash
export PUBMED_EMAIL="your@email.com"
```

#### "Reference mismatch at position X"

**Possible causes**:
1. Wrong reference genome version (use GRCh38)
2. Incorrect reference allele
3. Position not in reference

**Solution**: Verify variant coordinates and reference genome

### Performance Optimization

#### Slow First Run
- **Cause**: Model loading takes time
- **Solution**: Model stays in memory for subsequent runs

#### Batch Processing Tips
```bash
# Process without Evo2 for speed
python batch_analyze.py variants.csv --skip-evo2

# Use smaller model
# Edit config.py: EVO2_MODEL = "evo2_1b_base"
python batch_analyze.py variants.csv
```

### Debug Mode

Enable detailed logging:
```bash
# Check log file
tail -f outputs/ttn_agent.log

# Run with Python debugger
python -m pdb main.py 2-178612477-T-A
```

## API Reference

### Core Functions

#### `parse_variant(variant_string: str) -> Dict`
Parse variant notation into structured format.

**Parameters**:
- `variant_string`: Variant in format chr-pos-ref-alt

**Returns**:
```python
{
    'chrom': '2',
    'pos': 178612477,
    'ref': 'T',
    'alt': 'A',
    'variant_id': '2-178612477-T-A'
}
```

#### `Evo2Predictor.predict(variant_info: Dict) -> Dict`
Predict variant pathogenicity.

**Returns**:
```python
{
    'success': True,
    'prediction': 'pathogenic',
    'delta_score': -0.002667,
    'ref_score': -2.860429,
    'var_score': -2.863096,
    'confidence': 0.002667
}
```

#### `PubMedSearcher.search(variant_info: Dict) -> List[Dict]`
Search PubMed for relevant articles.

**Returns**: List of article dictionaries with fields:
- `pmid`: PubMed ID
- `title`: Article title
- `authors`: Author list
- `journal`: Journal name
- `year`: Publication year
- `phenotype`: Extracted phenotype
- `inheritance`: Inheritance pattern
- `age_onset`: Age of onset
- `pubmed_link`: URL to article

### Configuration Constants

See `config.py` for all configurable parameters:
- `EVO2_MODEL`: Model version
- `PATHOGENIC_THRESHOLD`: Classification threshold
- `PUBMED_MAX_RESULTS`: Number of articles
- `TTN_DOMAINS`: Domain coordinates
- `TTN_TRANSCRIPTS`: Transcript information

## Best Practices

### For Researchers
1. Always verify predictions with additional evidence
2. Review literature findings carefully
3. Consider transcript-specific effects
4. Compare with ClinVar when available

### For Clinical Use
⚠️ **Important**: This tool is for research only
- Do not use for clinical decision-making
- Validate findings with approved clinical tests
- Consult genetic counselors and specialists

### Data Management
- Keep organized records of analyses
- Document prediction versions used
- Save HTML reports for future reference
- Back up batch processing results

## Getting Help

### Resources
- **README.md**: General overview
- **QUICKSTART.md**: Fast setup guide
- **This guide**: Comprehensive documentation
- **test_agent.py**: Component testing

### Support
- Check logs: `outputs/ttn_agent.log`
- Run tests: `python test_agent.py`
- Review examples: `example_variants.csv`

### Reporting Issues
When reporting issues, include:
1. Variant being analyzed
2. Command used
3. Error message
4. Log file contents
5. Python/system version

## Appendix

### Supported Variant Types
Currently supports:
- ✅ Single Nucleotide Variants (SNVs)
- ❌ Insertions/Deletions (indels) - coming soon
- ❌ Complex variants - coming soon

### TTN Gene Information
- **Chromosome**: 2
- **Location**: 178,525,989-178,807,423 (GRCh38)
- **Size**: ~281 kb
- **Protein length**: Up to 34,350 amino acids
- **Function**: Muscle elasticity and contraction

### Citations
If using this tool in research, please cite:
1. Evo2 model: [Arc Institute](https://github.com/ArcInstitute/evo2)
2. Reference genome: GRCh38 (NCBI)
3. Literature source: PubMed (NCBI)

---

**Version**: 1.0.0  
**Last Updated**: 2024  
**License**: Research Use Only