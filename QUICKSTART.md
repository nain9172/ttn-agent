# Quick Start Guide

## 🚀 Installation (5 minutes)

```bash
# 1. Run setup script
chmod +x setup.sh
./setup.sh

# 2. Activate environment
source venv/bin/activate

# 3. Test installation
python test_agent.py
```

## 📝 Basic Usage

### Analyze a single variant

```bash
python main.py 2-178612477-T-A
```

### View help

```bash
python main.py --help
```

## 🎯 Example Variants

### Pathogenic Example
```bash
python main.py 2-178527121-A-G
```
Expected: **PATHOGENIC** (delta_score < 0)

### Benign Example
```bash
python main.py 2-178527628-C-T
```
Expected: **BENIGN** (delta_score ≥ 0)

## 📊 Understanding the Output

### Console Output
```
Starting TTN Variant AI Agent for variant: 2-178612477-T-A
Step 1: Parsing variant...
Step 2: Running Evo2 prediction...
Step 3: Searching PubMed...
Step 4: Generating protein schematic...
Step 5: Generating HTML report...
✅ Report generated successfully: outputs/variant_report_TIMESTAMP.html
```

### HTML Report Sections

1. **Variant Summary**
   - Genomic coordinates
   - Reference/alternate alleles
   - Gene information

2. **Evo2 Prediction**
   - Pathogenicity prediction
   - Delta score
   - Confidence level

3. **Protein Domain Localization**
   - Visual schematic
   - Domain mapping (Z-disk, I-band, A-band, M-band)
   - Transcript isoforms

4. **Literature Review**
   - Relevant PubMed articles
   - Phenotype information
   - Inheritance patterns
   - Age of onset

## ⚙️ Configuration

### Skip Evo2 (faster, but no prediction)
```bash
python main.py 2-178527628-C-T --skip-evo2
```

### Custom output location
```bash
conda activate ai && python main.py 2-178529960-G-T --output my_report.html
```

### Set PubMed credentials
```bash
export PUBMED_EMAIL="your_email@example.com"
export PUBMED_API_KEY="your_api_key"  # Optional
```

Or edit `config.py`:
```python
PUBMED_EMAIL = "your_email@example.com"
PUBMED_API_KEY = "your_api_key"  # Optional
```
cd /home/ryan910702/ttn_agent && python3 -c "
from utils.variant_parser import parse_variant
from utils.clinvar_parser import ClinVarParser

# 測試兩個不同的變異
variants = ['2-178529960-G-T', '2-178527121-T-C']

parser = ClinVarParser()

for var_str in variants:
    print(f'\n測試變異: {var_str}')
    print('='*60)
    variant_info = parse_variant(var_str)
    result = parser.parse_variant(variant_info)
    
    if result and result.get('pmid_list'):
        print(f'  找到 {len(result[\"pmid_list\"])} 個 PMIDs')
        print(f'  影響類型: {result[\"variant_impact\"]}')
    else:
        print('  未找到 ClinVar 數據')
"
## 🔧 Troubleshooting

### Out of Memory
Use smaller model in `config.py`:
```python
EVO2_MODEL = "evo2_1b_base"  # Instead of evo2_7b_base
```

### Reference Genome Not Found
```bash
cd data
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz
```

### PubMed Connection Issues
Check your internet connection and verify email is set:
```python
# In config.py
PUBMED_EMAIL = "your_email@example.com"
```

## 📈 Performance Tips

1. **First run is slow**: Model loading takes time
2. **Subsequent runs are faster**: Model stays in memory
3. **Use smaller model**: For testing, use `evo2_1b_base`
4. **Skip Evo2**: Use `--skip-evo2` for quick literature reviews

## 🎓 Understanding Scores

### Delta Score Interpretation

| Score Range | Interpretation | Confidence |
|-------------|----------------|------------|
| < -0.005 | Likely pathogenic | High |
| -0.005 to -0.001 | Possibly pathogenic | Medium |
| -0.001 to 0 | Borderline pathogenic | Low |
| 0 to 0.001 | Borderline benign | Low |
| > 0.001 | Likely benign | Medium-High |

### Example Results

**Pathogenic Variant** (2-178527121-A-G)
```
Delta Score: -0.000238
Prediction: PATHOGENIC
```

**Benign Variant** (2-178527628-C-T)
```
Delta Score: 0.000304
Prediction: BENIGN
```

## 📚 Next Steps

1. Analyze your variants of interest
2. Compare predictions with ClinVar
3. Review literature findings
4. Visualize protein domain locations
5. Generate comprehensive reports

## 🆘 Getting Help

- Check logs: `outputs/ttn_agent.log`
- Run tests: `python test_agent.py`
- View README: `README.md`
- Report issues: [GitHub Issues]

## ⚠️ Important Reminder

This tool is for **research purposes only**. Do not use for clinical decision-making without proper validation.