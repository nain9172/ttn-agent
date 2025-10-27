# TTN Variant AI Agent

A comprehensive AI agent for analyzing TTN (Titin) gene variants, combining pathogenicity prediction (Evo2), literature mining (PubMed), and protein domain visualization.

## 🌟 Features

- **Evo2 Pathogenicity Prediction**: Uses state-of-the-art Evo2 model for variant impact assessment
- **PubMed Literature Mining**: Automatically searches and extracts relevant clinical information
- **Protein Domain Visualization**: Generates detailed schematics showing variant location across TTN isoforms
- **Comprehensive HTML Reports**: Beautiful, interactive reports with all analysis results

## 📋 Requirements

- Python 3.8+
- Reference genome (GRCh38)
- Internet connection (for PubMed searches)

## 🚀 Quick Start

### 1. Installation

```bash
# Run setup script
chmod +x setup.sh
./setup.sh

# Activate environment
source venv/bin/activate
```

### 2. Basic Usage

```bash
# Analyze a variant
python main.py 2-178612477-T-A

# View the generated report
# Opens in outputs/variant_report_*.html
```

### 3. Example Variants

```bash
# Pathogenic example
python main.py 2-178527121-A-G

# Benign example
python main.py 2-178527628-C-T
```

## 📖 Detailed Installation

### Prerequisites

```bash
# Check Python version
python3 --version  # Should be 3.8+

# Install system dependencies (Ubuntu/Debian)
sudo apt-get install python3-venv python3-pip

# macOS
brew install python3
```

### Manual Installation

```bash
# 1. Create virtual environment
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# 2. Install dependencies
pip install -r requirements.txt

# 3. Install Evo2
pip install git+https://github.com/ArcInstitute/evo2.git

# 4. Download reference genome
mkdir -p data
cd data
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz
cd ..

# 5. Configure PubMed (optional but recommended)
export PUBMED_EMAIL="your_email@example.com"
```

## 💻 Usage

### Command Line Interface

```bash
python main.py <variant> [OPTIONS]

Arguments:
  variant              Variant in format: chr-pos-ref-alt
                      Example: 2-178612477-T-A

Options:
  --output PATH       Output HTML file path
  --skip-evo2        Skip Evo2 prediction (faster)
  --help             Show help message
```

### Input Format

Variants can be specified with different separators:
- `2-178612477-T-A` (hyphen)
- `2_178612477_T_A` (underscore)
- `2:178612477:T:A` (colon)

Format: `chromosome-position-reference-alternate`

### Examples

**Standard Analysis:**
```bash
python main.py 2-178527121-A-G
```

**Fast Mode (skip Evo2):**
```bash
python main.py 2-178527121-A-G --skip-evo2
```

**Custom Output:**
```bash
python main.py 2-178527121-A-G --output my_analysis.html
```

### Batch Processing

Process multiple variants from a file:

```bash
# Create input file (CSV or text)
cat > variants.txt << EOF
2-178527121-A-G
2-178527148-A-T
2-178527628-C-T
EOF

# Run batch analysis
python batch_analyze.py variants.txt

# Results in outputs/batch_TIMESTAMP/
```

## 📊 Output

### HTML Report Sections

1. **Variant Summary**
   - Genomic coordinates (GRCh38)
   - Reference/alternate alleles
   - Gene information

2. **Evo2 Pathogenicity Prediction**
   - Prediction: PATHOGENIC or BENIGN
   - Delta score (variant - reference)
   - Confidence level
   - Interpretation guide

3. **Protein Domain Localization**
   - Visual schematic of TTN protein
   - Variant location across domains (Z-disk, I-band, A-band, M-band)
   - Multiple transcript isoforms (N2BA, N2B, Novex-1/2/3)

4. **Literature Review (PubMed)**
   - Relevant articles
   - Phenotype (cardiac/skeletal)
   - Inheritance pattern
   - Age of onset
   - Links to papers

### Understanding Predictions

**Evo2 Delta Score:**
- **< 0** → PATHOGENIC (variant disrupts protein function)
- **≥ 0** → BENIGN (variant maintains protein function)

**Confidence:**
| Absolute Delta | Confidence |
|---------------|------------|
| > 0.005       | High       |
| 0.001-0.005   | Medium     |
| < 0.001       | Low        |

### Output Files

```
outputs/
├── variant_report_20240101_120000.html  # Main report
├── images/
│   └── titin_schematic_*.png           # Protein schematic
└── ttn_agent.log                        # Detailed logs
```

## 🔧 Configuration

Edit `config.py` to customize:

```python
# Model selection
EVO2_MODEL = "evo2_1b_base"  # Fast (16GB RAM)
# EVO2_MODEL = "evo2_7b_base"  # Accurate (64GB RAM)

# PubMed settings
PUBMED_EMAIL = "your@email.com"
PUBMED_MAX_RESULTS = 20

# Report styling
REPORT_HEADER_COLOR = "#2C3E50"
REPORT_ACCENT_COLOR = "#3498DB"
```

## 🧪 Testing

Run the test suite:

```bash
python test_agent.py
```

This tests:
- Module imports
- Variant parser
- Evo2 predictor
- PubMed search
- Image generator
- Complete pipeline

## 🐛 Troubleshooting

### Common Issues

**"Reference genome not found"**
```bash
cd data
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz
```

**"Out of memory"**
```python
# Use smaller model in config.py
EVO2_MODEL = "evo2_1b_base"

# Or skip Evo2
python main.py <variant> --skip-evo2
```

**"PubMed search failed"**
```python
# Set email in config.py
PUBMED_EMAIL = "your@email.com"

# Or as environment variable
export PUBMED_EMAIL="your@email.com"
```

**"Failed to load Evo2 model"**
```bash
pip uninstall evo2
pip install git+https://github.com/ArcInstitute/evo2.git
```

## 📚 Documentation

- **QUICKSTART.md**: Fast setup guide
- **USER_GUIDE.md**: Comprehensive manual
- **This README**: Overview

## 🔬 How It Works

1. **Variant Parsing**: Validates and standardizes input
2. **Evo2 Prediction**: 
   - Extracts 8kb window around variant
   - Scores reference and variant sequences
   - Calculates delta score
3. **PubMed Search**: 
   - Searches for relevant literature
   - Extracts clinical information
4. **Visualization**: 
   - Maps variant to protein domains
   - Generates schematic for all isoforms
5. **Report Generation**: 
   - Compiles all results
   - Generates HTML report

## ⚠️ Important Notes

### Limitations

- **Research Use Only**: Not for clinical decisions
- **SNVs Only**: Currently supports single nucleotide variants
- **TTN Gene**: Designed specifically for TTN variants
- **Approximate Positions**: Protein positions are estimates

### Requirements

- 16GB RAM recommended for Evo2
- Internet connection for PubMed
- GRCh38 reference genome (~1GB)

## 🤝 Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## 📝 License

This project is for research and educational purposes.

## 🙏 Acknowledgments

- **Arc Institute**: Evo2 model
- **NCBI**: PubMed database
- **CardioDB**: TTN domain information

## 📧 Support

For issues or questions:
- Check `outputs/ttn_agent.log`
- Run `python test_agent.py`
- Review documentation

## 🔗 References

1. [Evo2 GitHub](https://github.com/ArcInstitute/evo2)
2. [TTN Gene - GeneCards](https://www.genecards.org/cgi-bin/carddisp.pl?gene=TTN)
3. [CardioDB TTN Info](https://www.cardiodb.org/titin/titin_transcripts.php)
4. [PubMed](https://pubmed.ncbi.nlm.nih.gov/)

---

**Version**: 1.0.0  
**Last Updated**: 2024  
**Status**: Production Ready