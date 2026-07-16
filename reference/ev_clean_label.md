# Clean a pathway name for ring display

Strips the database prefix, expands common acronyms, title-cases, wraps,
and applies MitoCarta-hierarchy shortening for `MITOCARTA_` pathways.

## Usage

``` r
ev_clean_label(name)
```

## Arguments

- name:

  Character vector of raw pathway names.

## Value

Character vector of cleaned, wrapped labels.

## Examples

``` r
ev_clean_label(c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "REACTOME_TCA_CYCLE"))
#> [1] "OXPHOS"    "TCA Cycle"
```
