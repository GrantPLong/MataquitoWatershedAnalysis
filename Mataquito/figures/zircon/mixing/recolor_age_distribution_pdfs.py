# recolor_age_distribution_pdfs.py
#
# Purpose
# -------
# Recolor the parent colors in the detritalPy-mix age_distribution figures
# (output of dMix.plotMix in all_mixing_scenarios_seeded.py) to a
# colorblind-friendly pair, in response to Reviewer 2. The mixing results are
# unchanged; only the RGB color operators in the vector PDFs are replaced.
#
# Background
# ----------
# matplotlib writes each fill/stroke color into the PDF content stream as
#   "r g b rg"   (fill color)   or   "r g b RG"   (stroke color)
# with r, g, b in [0, 1]. The original figures used the detritalPy 'Default'
# palette for the two parents:
#   parent 1: darkgreen  = (0, 100, 0)/255     -> "0 0.3921568627 0"
#   parent 2: firebrick  = (178, 34, 34)/255   -> "0.6980392157 0.1333333333 0.1333333333"
# These two triplets are the only RGB colors present in the age_distribution
# PDFs (all other elements are black/white, written as grayscale operators).
# They are replaced with:
#   parent 1: darkorange = (255, 140, 0)/255   -> "1 0.5490196078 0"
#   parent 2: steelblue  = (70, 130, 180)/255  -> "0.2745098039 0.5098039216 0.7058823529"
# Transparency is stored separately (ExtGState), so it is preserved.
#
# Usage
# -----
#   pip install pikepdf
#   python recolor_age_distribution_pdfs.py
#
# Output location is set by OUT_DIR below:
#   TEST RUN  : OUT_DIR = os.path.join(FIG_DIR, 'recolored')  -> originals untouched
#   FINAL RUN : OUT_DIR = FIG_DIR                              -> originals overwritten
# If the originals are overwritten, they can be restored from git with:
#   git checkout -- Mataquito/figures/zircon/mixing/

import glob
import os
import re

import pikepdf

# Folder containing the mixing figures
FIG_DIR = '/Users/Glong1/Desktop/Andes/Chapter_1_Manuscript/MataquitoWatershedAnalysis/Mataquito/figures/zircon/mixing'

# Folder to write recolored files to (see Usage above). Default is the TEST RUN setting.
OUT_DIR = os.path.join(FIG_DIR, 'recolored')
os.makedirs(OUT_DIR, exist_ok=True)   # create the folder if it does not exist

# Old RGB triplet -> new RGB triplet.
# matplotlib wraps long lines in the content stream, so a triplet can be split
# across a newline (e.g. "0.6980392157\n0.1333333333 0.1333333333 rg"). Each
# pattern therefore allows any whitespace (\s+) between the three components,
# and requires the triplet to be followed by a fill (rg) or stroke (RG)
# operator. The lookbehind (?<![\d.]) prevents matching inside a longer number.
COLOR_MAP = [
    (re.compile(rb'(?<![\d.])0\s+0\.3921568627\s+0(?=\s+(?:rg|RG)\b)'),
     b'1 0.5490196078 0'),                                   # darkgreen -> darkorange
    (re.compile(rb'(?<![\d.])0\.6980392157\s+0\.1333333333\s+0\.1333333333(?=\s+(?:rg|RG)\b)'),
     b'0.2745098039 0.5098039216 0.7058823529'),             # firebrick -> steelblue
]


def recolor_stream(stream):
    """Replace old color triplets in one content stream; return number of replacements."""
    data = stream.read_bytes()          # decompressed stream content
    n = 0
    for pattern, new in COLOR_MAP:
        data, k = pattern.subn(new, data)   # replace all matches; k = number replaced
        n += k
    if n > 0:
        stream.write(data)              # write back (pikepdf recompresses on save)
    return n


# Process only the age_distribution figures directly in FIG_DIR
# (glob is not recursive, so files already in OUT_DIR are never re-read)
files = sorted(glob.glob(os.path.join(FIG_DIR, '*age_distribution*.pdf')))
print(f'{len(files)} age_distribution files found; writing to {OUT_DIR}')

for f in files:
    pdf = pikepdf.open(f, allow_overwriting_input=True)
    n_total = 0
    for page in pdf.pages:
        # Page content may be a single stream or an array of streams
        contents = page.Contents
        streams = list(contents) if isinstance(contents, pikepdf.Array) else [contents]
        # Form XObjects (reusable drawing objects, e.g. markers) also carry color operators
        xobjects = page.Resources.get('/XObject', {})
        streams += [x for _, x in xobjects.items() if x.get('/Subtype') == '/Form']
        for s in streams:
            n_total += recolor_stream(s)
    pdf.save(os.path.join(OUT_DIR, os.path.basename(f)))   # same file name, in OUT_DIR
    pdf.close()
    # A count of 0 means the file did not contain the expected colors -- inspect it manually
    print(f'{n_total:5d} replacements  {os.path.basename(f)}')
