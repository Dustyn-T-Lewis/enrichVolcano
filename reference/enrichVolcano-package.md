# enrichVolcano: Volcano-in-Ring and Concordance Plots from Enrichment Results

Turns differential abundance results into gene-set enrichment figures.
Reads results from common proteomics tools and study files with a sample
sheet and contrasts, runs 'fgsea' and limma's camera and fry tests, or
converts enrichment results computed elsewhere, into one validated
object. From it the package draws a differential abundance volcano
inside a ring of enrichment terms, and a scatter comparing two contrasts
term by term, and writes them to one PDF. Redundant terms can be flagged
for display with the EnrichmentMap overlap rule or 'fgsea' conditional
collapsing. No p-value changes. Seven example proteomics studies are
included.

## See also

Useful links:

- <https://github.com/Dustyn-T-Lewis/enrichVolcano>

- <https://Dustyn-T-Lewis.github.io/enrichVolcano/>

- Report bugs at
  <https://github.com/Dustyn-T-Lewis/enrichVolcano/issues>

## Author

**Maintainer**: Dustyn Lewis <dtlts14@gmail.com>

Authors:

- Dustyn Lewis <dtlts14@gmail.com>
