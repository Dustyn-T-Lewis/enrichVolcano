# Clean a pathway name for display

Strips the database prefix, expands common acronyms, title-cases, wraps,
and applies MitoCarta-hierarchy shortening for `MITOCARTA_` pathways.

## Usage

``` r
clean_label(name, width = 15)
```

## Arguments

- name:

  Character vector of raw pathway names.

- width:

  Wrap width in characters. The default suits ring arcs; at 15 or
  narrower a few short names also get hand-placed line breaks.

## Value

Character vector of cleaned, wrapped labels.

## Examples

``` r
clean_label(c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "REACTOME_TCA_CYCLE"))
#> [1] "OXPHOS"    "TCA Cycle"
```
