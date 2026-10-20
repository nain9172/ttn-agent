"""網站的 Information 頁（site/info.html）：結構 band、isoform 的功能與臨床關聯、各種分數與推論的說明。

build_site.py 建站時寫出；報告頁的導覽列與側邊目錄都連到這一頁。內容是給臨床端看的英文說明，
Evo2 的段落與報告裡 Evo2 區塊的「?」提示共用同一段文字（utils/html_report.EVO2_HELP_HTML）。
"""

from __future__ import annotations

from utils.html_report import EVO2_HELP_HTML, FONT_STACK

BANDS = [
    ("Z-disk", "#FF6B6B", "exons 1–28",
     "Anchors the N-terminus of titin in the Z-disk. The Z1Z2 Ig domains bind telethonin (T-cap) and the "
     "Z-repeats bind &alpha;-actinin, cross-linking titin from neighbouring sarcomeres.",
     "Constitutively expressed in heart and skeletal muscle. Truncating variants here are found in dilated "
     "cardiomyopathy (DCM); missense variants in the telethonin-binding region have been reported in "
     "hypertrophic and dilated cardiomyopathy, with limited evidence."),
    ("I-band", "#4ECDC4", "exons 29–252",
     "The extensible &ldquo;molecular spring&rdquo; that sets passive tension: tandem Ig segments, the cardiac "
     "N2B element (its unique sequence is encoded by exon 49), the N2A element (binds calpain-3 and "
     "CARP/ANKRD1) and the PEVK region (rich in proline, glutamate, valine and lysine). Most "
     "tissue-specific alternative splicing of titin happens here.",
     "Many I-band exons are spliced out of cardiac titin (low heart PSI). Heterozygous truncating variants in "
     "those exons are not expected to cause DCM; when biallelic they cause recessive skeletal or congenital "
     "myopathies (e.g. centronuclear myopathy, arthrogryposis). Truncations in I-band exons that are kept in "
     "the heart carry a smaller DCM risk than A-band truncations."),
    ("A-band", "#45B7D1", "exons 253–358",
     "The inextensible part of titin that lies along the thick filament. Regular super-repeats of Ig and "
     "fibronectin-3 domains bind myosin and myosin-binding protein C and act as a ruler for thick-filament "
     "assembly.",
     "Present in every isoform. Truncating variants here (A-band TTNtv) are the most common genetic cause of "
     "DCM, inherited as autosomal dominant with incomplete penetrance; they also predispose to peripartum, "
     "alcohol- and chemotherapy-related cardiomyopathy. Missense variants in the A150 fibronectin-3 domain "
     "(exon 344) cause hereditary myopathy with early respiratory failure (HMERF, autosomal dominant)."),
    ("M-band", "#96CEB4", "exons 359–364 (Mex1–Mex6)",
     "The C-terminus overlaps in the M-band from opposite half-sarcomeres. This region includes the titin "
     "kinase domain at the A/M-band junction (a stretch sensor for signalling) and binding sites for "
     "calpain-3, obscurin and myomesin; Mex5 (exon 363, &ldquo;is7&rdquo;) is alternatively spliced.",
     "Mex6 (exon 364) variants cause tibial muscular dystrophy (TMD, Udd myopathy) when heterozygous and "
     "limb-girdle muscular dystrophy 2J / R10 when homozygous (e.g. the Finnish FINmaj variant). Biallelic "
     "truncations in Mex1–Mex3 cause Salih myopathy (early-onset myopathy with fatal cardiomyopathy, "
     "autosomal recessive). These final exons escape part of nonsense-mediated decay, so the mechanism differs "
     "from A-band TTNtv."),
]

ISOFORMS = [
    ("Meta-transcript", "Reference", "#7f8c8d",
     "An inferred transcript containing all 364 exons (LRG_391). Not a real protein; it is the reference "
     "for exon numbering used throughout these reports."),
    ("N2BA", "Cardiac", "#e74c3c",
     "Long, compliant cardiac isoform containing both the N2B and N2A elements and a longer PEVK region; "
     "includes exon 49. The N2BA : N2B ratio is high in the fetal heart and rises again in failing / dilated "
     "hearts, making the myocardium more compliant."),
    ("N2B", "Cardiac", "#e74c3c",
     "Short, stiff cardiac isoform with the N2B element (exon 49) but no N2A element. Predominant in the "
     "adult left ventricle and responsible for most of its passive stiffness."),
    ("N2A", "Skeletal muscle", "#e67e22",
     "Skeletal-muscle isoform: contains the N2A element and a long, muscle-specific PEVK region, and skips "
     "exon 49. Variants in exons used only by N2A mainly affect skeletal muscle."),
    ("Novex-1 / Novex-2", "Minor cardiac", "#95a5a6",
     "Minor cardiac isoforms that each include one extra &ldquo;novex&rdquo; I-band exon (exon 45 or 46). "
     "Their function is unclear."),
    ("Novex-3", "Short isoform", "#95a5a6",
     "A ~700 kDa isoform that ends in the proximal I-band with its own C-terminal exon (exon 48); it binds "
     "obscurin. Variants after exon 48 do not affect Novex-3."),
]

GLOSSARY = [
    ("TTNtv", "TTN truncating variant: nonsense, frameshift or canonical splice-site variant."),
    ("PSI", "Percent spliced in: the percentage of a tissue&rsquo;s TTN transcripts that include an exon. "
            "Exons with heart PSI &ge; 90% are treated as constitutively expressed in the heart."),
    ("AD / AR", "Autosomal dominant (one pathogenic allele is enough) / autosomal recessive (two pathogenic "
                "alleles, homozygous or compound heterozygous, are needed)."),
    ("DCM / HCM", "Dilated / hypertrophic cardiomyopathy."),
    ("ARVC / LVNC", "Arrhythmogenic right ventricular cardiomyopathy / left ventricular non-compaction."),
    ("TMD", "Tibial muscular dystrophy (Udd distal myopathy), autosomal dominant, Mex6."),
    ("LGMD2J / LGMDR10", "Limb-girdle muscular dystrophy caused by biallelic Mex6 variants."),
    ("HMERF", "Hereditary myopathy with early respiratory failure, autosomal dominant, exon 344."),
    ("CNM", "Centronuclear myopathy; the TTN form is autosomal recessive."),
    ("HGVS c.", "Coding-DNA notation on transcript NM_001267550.2 (the meta-transcript), e.g. c.107743A&gt;C."),
]


def info_html() -> str:
    bands = "".join(
        f"""<div class="card" style="border-top-color:{color}">
            <h3><span class="sw" style="background:{color}"></span>{name} <small>{exons}</small></h3>
            <p><b>Function.</b> {func}</p>
            <p><b>Clinical associations.</b> {clin}</p>
        </div>"""
        for name, color, exons, func, clin in BANDS
    )
    isoforms = "".join(
        f"""<tr><td><b>{name}</b></td><td><span class="tag" style="border-color:{color};color:{color}">{tissue}</span></td>
            <td>{desc}</td></tr>"""
        for name, tissue, color, desc in ISOFORMS
    )
    glossary = "".join(f"<dt>{k}</dt><dd>{v}</dd>" for k, v in GLOSSARY)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>TTN Report Information</title>
<style>
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{ font-family: {FONT_STACK}; color: #24303f; background: #f4f6fa; line-height: 1.6; }}
header {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: #fff; padding: 30px 24px 26px; }}
header .inner, main {{ max-width: 1100px; margin: 0 auto; }}
header h1 {{ font-size: 1.8em; }}
header p {{ opacity: 0.9; margin-top: 6px; }}
header a {{ color: #fff; font-weight: 600; text-decoration: none; border: 1px solid rgba(255,255,255,0.5);
            border-radius: 999px; padding: 4px 14px; display: inline-block; margin-top: 12px; font-size: 0.9em; }}
main {{ padding: 24px 24px 60px; }}
nav.jump {{ display: flex; flex-wrap: wrap; gap: 8px; margin-bottom: 22px; }}
nav.jump a {{ background: #fff; border-radius: 999px; padding: 5px 14px; text-decoration: none; color: #3d4d63;
              box-shadow: 0 2px 8px rgba(31,45,68,0.08); font-size: 0.9em; }}
section {{ background: #fff; border-radius: 12px; padding: 24px 26px; margin-bottom: 22px;
           box-shadow: 0 6px 20px rgba(31,45,68,0.07); scroll-margin-top: 12px; }}
h2 {{ font-size: 1.4em; color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 8px; margin-bottom: 16px; }}
h3 {{ font-size: 1.1em; margin-bottom: 8px; display: flex; flex-wrap: wrap; align-items: baseline; gap: 4px 8px; }}
h3 .sw {{ align-self: center; }}
h3 small {{ font-weight: 400; color: #6b7a90; }}
p {{ margin-bottom: 10px; }}
ul {{ margin: 6px 0 10px 22px; }}
.cards {{ display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 16px; }}
.card {{ border: 1px solid #e4e9f0; border-top: 5px solid; border-radius: 10px; padding: 16px 18px; }}
.sw {{ width: 14px; height: 14px; border-radius: 3px; display: inline-block; border: 1px solid rgba(0,0,0,0.3); }}
table {{ width: 100%; border-collapse: collapse; }}
td {{ padding: 10px; border-bottom: 1px solid #eef1f6; vertical-align: top; }}
.tag {{ display: inline-block; border: 1px solid; border-radius: 999px; padding: 1px 10px; font-size: 0.82em;
        font-weight: 700; white-space: nowrap; }}
.note {{ background: #f4f8fc; border: 1px solid #d6e4f0; border-radius: 8px; padding: 12px 16px; }}
.warn {{ background: #fff4e0; border: 1px solid #f0ad4e; border-radius: 8px; padding: 12px 16px; color: #6b4500; }}
code {{ background: #eef1f6; border-radius: 4px; padding: 1px 6px; font-size: 0.92em; }}
dl {{ display: grid; grid-template-columns: max-content 1fr; gap: 6px 18px; }}
dt {{ font-weight: 700; }}
@media (max-width: 760px) {{ .cards {{ grid-template-columns: 1fr; }} }}
@media (max-width: 640px) {{ dl {{ grid-template-columns: 1fr; }} section {{ padding: 18px 16px; }} }}
</style>
</head>
<body>
<header><div class="inner">
    <h1>Information</h1>
    <p>Background for reading the TTN variant reports: titin structure, isoforms, scores and what is inferred.</p>
    <a href="index.html">&larr; All reports</a>
</div></header>
<main>
<nav class="jump">
    <a href="#ids">Variant IDs</a><a href="#bands">Structural bands</a><a href="#isoforms">Isoforms</a>
    <a href="#evo2">Evo2 score</a><a href="#databases">Databases</a><a href="#inference">Inference</a>
    <a href="#glossary">Glossary</a>
</nav>

<section id="ids">
    <h2>Variant IDs and strand</h2>
    <p>TTN lies on the <b>minus strand</b> of chromosome 2. Variant IDs on this site use the standard
    <b>GRCh38 genomic plus-strand</b> notation <code>chr-position-ref-alt</code>, the same as VCF files, ClinVar
    and gnomAD. Each report also gives the change on the TTN coding strand as HGVS <code>c.</code> notation on
    NM_001267550.2, where the bases are complementary (A&harr;T, C&harr;G).</p>
    <p class="note">Example: <code>2-178527245-T-G</code> (genomic) is <code>NM_001267550.2:c.107743A&gt;C</code>
    (coding strand). When you submit variants, either orientation is accepted; it is detected from the
    reference sequence and converted automatically.</p>
</section>

<section id="bands">
    <h2>Structural bands of titin</h2>
    <p>One titin molecule spans half a sarcomere, from the Z-disk to the M-band. The reports place each
    variant in one of four regions, numbered by meta-transcript exon.</p>
    <div class="cards">{bands}</div>
</section>

<section id="isoforms">
    <h2>Isoforms</h2>
    <p>Alternative splicing, mostly in the I-band, produces tissue-specific isoforms. The schematic in each
    report draws them on one genomic axis, with cardiac isoforms framed in red and the skeletal-muscle isoform
    in orange, and marks <b>exon 49</b> (hg38 chr2:178,741,921&ndash;178,739,141, 2,781 bp encoding the N2B
    unique sequence): N2BA and N2B contain it, N2A skips it.</p>
    <table>{isoforms}</table>
</section>

<section id="evo2">
    <h2>Evo2 score</h2>
    <div class="note">{EVO2_HELP_HTML}</div>
</section>

<section id="databases">
    <h2>Database annotations</h2>
    <p>The <b>Curated Database Annotations</b> section is a lookup, not a prediction. It lists what each source
    records for this exact variant (same position, reference and alternate allele):</p>
    <ul>
        <li><b>HGMD</b>: curated disease-causing (DM) variants reported in the literature.</li>
        <li><b>PubTator</b>: variants text-mined from PubMed articles about TTN.</li>
        <li><b>PubMind</b>: literature-mined variant&ndash;disease records.</li>
        <li><b>ClinVar</b>: the variant ID in each report header links to its ClinVar record when one exists.</li>
    </ul>
</section>

<section id="inference">
    <h2>What is inferred (not observed data)</h2>
    <p class="warn">Sections with an amber <b>INFERENCE</b> badge and hatched cards marked <b>INFERRED</b> are
    generated by a language model and rules. They are hypotheses for prioritising follow-up, not patient data,
    not database records and not evidence for ACMG classification.</p>
    <p><b>Mechanistic inference</b> runs when no article mentions the variant. It predicts a disease and tissue
    by analogy: the variant&rsquo;s band, exon and type, how strongly the exon is used in heart and skeletal muscle
    (PSI), the Evo2 call, and the nearest known pathogenic HGMD variants in the same exon or band.</p>
    <p><b>Mode of inheritance</b> estimates whether the variant would act as autosomal dominant or recessive,
    using the rules below, the disease, records for the variant itself and nearby HGMD variants; a language model
    makes the final call from this evidence.</p>
    <ul>
        <li>Truncating variants in exons kept in the heart, especially the A-band &rarr; dominant DCM
            (haploinsufficiency).</li>
        <li>Truncating variants in exons spliced out of the heart &rarr; recessive myopathy (needs a second allele).</li>
        <li>Mex6 (exon 364) &rarr; TMD when heterozygous, LGMD2J when homozygous (both).</li>
        <li>Missense in exon 344 &rarr; HMERF (dominant); other missense variants are usually recessive
            (in trans with a truncating variant) or of uncertain effect.</li>
    </ul>
    <p class="note"><b>How reliable is the inheritance call?</b> Tested leave-one-out on 1,712 HGMD TTN variants whose
    mode could be derived from the record (homozygous / compound-heterozygous comments, or a disease with a known
    mode; 1,530 AD, 182 AR): the rules agreed with the label for 82.8% of variants and the language model for 80.4%
    (always answering &ldquo;dominant&rdquo; would score 89.4%). Both recognise dominant variants well (recall 87&ndash;90%)
    but recessive ones poorly (24&ndash;26%), mainly because many recessive cases carry an ordinary A-band truncating
    allele that looks identical to a dominant DCM variant. Treat an AR call as a prompt to look for a second allele,
    and an AD call as uninformative about a possible second allele.</p>
    <p>Literature fields marked <b>LLM-generated</b> (disease, tissue, age, inheritance and the AI analysis of
    each article) were extracted by the language model from the article and should be checked against it.</p>
</section>

<section id="glossary">
    <h2>Glossary</h2>
    <dl>{glossary}</dl>
</section>
</main>
</body>
</html>"""
